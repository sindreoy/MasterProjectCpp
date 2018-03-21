//
// Created by Sindre Bakke Øyen on 18.03.2018.
//

#include "PBModel.h"

/* Constructors */
PBModel::PBModel() = default;

PBModel::PBModel(char const *f, realtype kb1, realtype kb2, realtype kc1, realtype kc2,
                 const Grid &g, const SystemProperties &s,
                 const Fluid &cont, const Fluid &disp):
        filename(f), grid(g), sysProps(s), cont(cont), disp(disp), cvode_mem(nullptr),
        kerns(Kernels(kb1, kb2, kc1, kc2, 1, g, s, cont, disp)){
    /* Get rows and columns of csv file and initialize M and N respectively */
    this->getRowsAndCols();

    /* Allocate memory for member variables */
    this->t      = gsl_vector_alloc(this->M);
    this->r      = gsl_vector_alloc(this->N);
    this->fv     = gsl_matrix_alloc(this->M, this->N);
    this->psi    = gsl_matrix_calloc(this->M, this->N);
    this->tau    = gsl_vector_alloc(this->M);
    this->NPsi   = N_VNew_Serial(this->N);

    /* Set experimental data: r, t, fv and rescale */
    this->getDistributions();
    this->rescaleInitial();

    /* Set final time of experiment and update kernels */
    realtype tf = gsl_vector_get(this->t, this->M - 1);
    this->kerns.setTf(tf);

    /* Assign nondimensional time tau = t / tf */
    gsl_vector_memcpy(this->tau , this->t);
    gsl_vector_scale(this->tau, 1/this->kerns.getTf());
    this->tRequested = gsl_vector_get(this->tau, 0);

    /* Assign nondimensional distribution psi = fv * Rm */
    gsl_vector_view fv0 = gsl_matrix_row(this->fv, 0);
    this->psiN = gsl_matrix_row(this->psi, 0);
    gsl_vector_memcpy(&(this->psiN.vector), &fv0.vector);
    gsl_matrix_scale(this->psi, sysProps.getRm());

    /* Assign initial condition to solution vector */
    NV_DATA_S(this->NPsi) = this->psiN.vector.data;

    /* Prepare memory for integration */
    this->prepareCVMemory();
}

/* Helper methods */
void PBModel::getRowsAndCols(){
    FILE *f = fopen(filename, "r");
    if (f != nullptr) {
        // Initialize variables
        size_t rows = 0, cols = 0;
        size_t i=0, j=1000, tmp = 0;
        bool flag = false;
        realtype val = 0;
        const char s[2] = ",";
        const int bufSize = 1000;
        char line[bufSize], *toFree;

        // Now loop over lines
        while(fgets(line, sizeof line, f) != nullptr) {
            rows++;
            tmp = 0;
            i = 0;
            toFree = strdup(line);
            while ((strsep(&toFree, s)) != nullptr) {
                tmp++;
                if (rows > TRASHROWS && tmp > TRASHCOLS){
                    if (toFree != nullptr) {
                        val = strtod(toFree, nullptr);
                    } else val = 1;
                    if (val != 0) flag = true;
                    else if (flag && val == 0){
                        i++;
                    }
                }
            }
            flag = false;
            if (rows > TRASHROWS && j > i) j = i;
            if (tmp > cols) {
                cols = tmp;
            }
        }
        if (j > TRUNCATETHRESHOLD){
            cols = cols - TRASHCOLS - j + TRUNCATETHRESHOLD;
        } else cols = cols - TRASHCOLS;
        rows = rows - TRASHROWS;
        this->M = rows;
        this->N = cols;
        fclose(f);
    } else {
        perror(filename);
    }
}

void PBModel::getDistributions() {
    /* Declare needed variables */
    char *hours, *minutes;
    realtype h, m;
    size_t i = 0, k = 0; // Index variables

    const char s[2] = ",";
    const size_t bufSize = 10*(this->N);
    char line[bufSize], *token, *toFree, *timeStr;

    /* Open file and start reading */
    FILE *f = fopen(filename, "r");
    if (f != nullptr) {
        while (fgets(line, sizeof line, f) != nullptr) {
            if (i == 0) {
                i++;
                continue;
            }
            toFree = strdup(line);
            k = 0;
            while((token = strsep(&toFree, s)) != nullptr){
                if (k < TRASHCOLS){
                    if (k == 2){
                        if (i==1){
                            k++;
                            continue;
                        } else {
                            // Fetch time column
                            timeStr = strsep(&token, " ");
                            timeStr = token;
                            hours = strsep(&timeStr, ":");
                            minutes = timeStr;
                            h = strtod(hours, nullptr);
                            m = strtod(minutes, nullptr);
                            gsl_vector_set(t, i-TRASHROWS, h*3600+m*60);
                            k++;
                            continue;
                        }
                    } else {
                        k++;
                        continue;
                    }
                }
                if (k-TRASHCOLS < this->N) {
                    switch (i) {
                        case 1: {
                            /* Given sizes are diameter; we need radii. Also they should be microns */
                            gsl_vector_set(r, k - TRASHCOLS, strtod(token, nullptr) / 2 * 1.e-6);
                            break;
                        }
                        default: {
                            gsl_matrix_set(fv, i - TRASHROWS, k - TRASHCOLS, strtod(token, nullptr));
                            break;
                        }
                    }
                    k++;
                } else continue;
            }
            i++;
        }
        fclose(f);
        gsl_vector_add_constant(t, -gsl_vector_get(t, 0));
        for (i = 1; i < t->size; i++){
            if (gsl_vector_get(t, i) - gsl_vector_get(t, i-1) == 0){
                gsl_vector_set(t, i, gsl_vector_get(t, i) + 30);
                // t_data[i] += 30;
            }
        }
        gsl_vector_set(t, 0, 0);
        if (gsl_vector_get(t, 1) - gsl_vector_get(t, 0) == 0){
            gsl_vector_set(t, 1, gsl_vector_get(t, 1) + 30);
        }
    } else perror(filename);
}

void PBModel::rescaleInitial(){
    /* Takes distributions, corresponding radii and phase fraction and rescales the distributions
     * fv = phi/I * f0
     * Approximate I by trapezoids: I = sum([r(i+1)-r(i)] * [f(i+1)+f(i)]) from i=0 to N-1
     */
    size_t i, j;
    realtype I, rj, rjj, fij, fijj;

    for (i = 0; i < M; i++){ /* Loop over rows */
        I  = 0;
        for ( j = 0; j < N-1; j++ ){ /* Loop over columns */
            rj  = gsl_vector_get(r, j); rjj = gsl_vector_get(r, j+1);
            fij = gsl_matrix_get(fv, i, j); fijj = gsl_matrix_get(fv, i, j+1);
            I  += ( rjj-rj ) * ( fijj+fij );
        }
        I = I/2;
        gsl_vector_view rowJ = gsl_matrix_row(fv, i);
        gsl_vector_scale(&rowJ.vector, PHI/I);
    }
}

int PBModel::prepareCVMemory(){
    int flag = 0;
    /* Call CVodeCreate to create the solver memory and specify the
     * Backward Differentiation Formula and the use of a Newton iteration */
    this->cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (checkFlag((void *)this->cvode_mem, "CVodeCreate", 0)) return(1);

    /* Call CVodeInit to initialize the integrator memory and specify the
     * user's right hand side function in y'=f(t,y), the inital time T0, and
     * the initial dependent variable vector y. */
    flag = CVodeInit(this->cvode_mem, this->getRHS, gsl_vector_get(this->tau, 0), this->NPsi);
    if (checkFlag(&flag, "CVodeInit", 1)) return(1);

    /* Call CVodeSStolerances to specify the scalar relative tolerance
     * and scalar absolute tolerance */
    flag = CVodeSStolerances(this->cvode_mem, RTOL, ATOL);
    if (checkFlag(&flag, "CVodeSStolerances", 1)) return(1);

    /* Create dense SUNMatrix for use in linear solves */
    this->A = SUNDenseMatrix(this->N, this->N);
    if(checkFlag((void *)this->A, "SUNDenseMatrix", 0)) return(1);

    /* Create dense SUNLinearSolver object for use by CVode */
    this->LS = SUNDenseLinearSolver(this->NPsi, this->A);
    if(checkFlag((void *)this->LS, "SUNDenseLinearSolver", 0)) return(1);

    /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
    flag = CVDlsSetLinearSolver(this->cvode_mem, this->LS, this->A);
    if(checkFlag(&flag, "CVDlsSetLinearSolver", 1)) return(1);
    std::cout << flag << std::endl;
    return flag;
}

int PBModel::checkFlag(void *flagvalue, const char *funcname, int opt){
    int *errflag;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL) {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return(1); }

        /* Check if flag < 0 */
    else if (opt == 1) {
        errflag = (int *) flagvalue;
        if (*errflag < 0) {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
                    funcname, *errflag);
            return(1); }}

        /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL) {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return(1); }

    return(0);
}

/* Solver methods */
int PBModel::getRHS(realtype t, N_Vector y, N_Vector ydot, void *user_data){
    // TODO: Do a bunch of math in GSL and copy the data in RHS vector to ydot.
    return 0;
}

/* getRHS will return ydot, such that N_Vector NPsi is updated. Therefore:
 * TODO: After each time iteration, copy the data of NPsi into correct row of gsl_matrix psi
 */
int PBModel::timeIterate() {
    return 0;
}

int PBModel::solvePBE(){
    return 0;
}

/* Getter methods */
gsl_matrix *PBModel::getFv() const {
    return fv;
}

gsl_vector *PBModel::getR() const {
    return r;
}

gsl_vector *PBModel::getT() const {
    return t;
}

size_t PBModel::getM() const {
    return M;
}

size_t PBModel::getN() const {
    return N;
}

const Grid &PBModel::getGrid() const {
    return grid;
}

const Kernels &PBModel::getKerns() const {
    return kerns;
}

const SystemProperties &PBModel::getSysProps() const {
    return sysProps;
}

const Fluid &PBModel::getCont() const {
    return cont;
}

const Fluid &PBModel::getDisp() const {
    return disp;
}

/* Printer methods */
void PBModel::printExperimentalDistribution(){
    size_t i, j;
    std::cout << "The droplet size density distribution:" << std::endl;
    for (i=0;i<M;i++){
        for(j=0;j<N;j++){
            std::cout << std::setw(4) << std::setprecision(4) << gsl_matrix_get(fv, i, j) << "\t";
        }
        std::cout << std::endl;
    }
}

void PBModel::printCurrentPsi(){
    realtype *data = NV_DATA_S(this->NPsi);
    size_t i = 0;
    std::cout << "Psi for the current time iteration is: " << std::endl;
    for (i = 0; i < this->N; i++){
        std::cout << std::setw(8) << std::setprecision(3) << data[i];
    }
    std::cout << std::endl;
}

void PBModel::printPsi(){
    size_t i = 0, j = 0;
    std::cout << "Nondimensionalized droplet size density distribution:" << std::endl;
    for (i = 0; i < this->M; i++){
        for (j = 0; j < this->N; j++){
            std::cout << std::setw(4) << std::setprecision(3) << gsl_matrix_get(psi, i, j) << "\t";
        }
        std::cout << std::endl;
    }
}

void PBModel::printDimensions(){
    std::cout << "Rows: " << this->M << ", Columns: " << this->N << std::endl;
}

void PBModel::printTime() {
    size_t i = 0;
    std::cout << "Time of measurement: " << std::endl;
    for (i = 0; i < this->M; i++){
        std::cout << std::setw(5) << std::setprecision(4) << gsl_vector_get(this->t, i);
    }
    std::cout << std::endl;
}

void PBModel::printTau(){
    size_t i = 0;
    std::cout << "Nondimensionalized time vector: " << std::endl;
    for (i = 0; i < this->M; i++){
        std::cout << std::setw(6) << std::setprecision(2) << gsl_vector_get(this->tau, i);
    }
    std::cout << std::endl;
}

void PBModel::printSizeClasses() {
    size_t i = 0;
    std::cout << "Measured size classes: " << std::endl;
    for (i = 0; i < this->N; i++){
        std::cout << std::setw(10) << std::setprecision(7) << gsl_vector_get(this->r, i);
    }
    std::cout << std::endl;
}

/* Destructors */
PBModel::~PBModel(){
    gsl_matrix_free(this->fv);
    gsl_matrix_free(this->psi);
    gsl_vector_free(this->r);
    gsl_vector_free(this->t);
    gsl_vector_free(this->tau);
    /* NPsi->data points to a row in psi. psi is freed, so we cannot free NPsi yet.
     * Point NPsi->data to nullptr before freeing, so we don't encounter memory issues.
     */
    NV_DATA_S(this->NPsi) = nullptr;
    N_VDestroy_Serial(this->NPsi);
}
