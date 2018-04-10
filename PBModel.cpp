//
// Created by Sindre Bakke Øyen on 18.03.2018.
//

#include <gsl/gsl_vector_double.h>
#include "PBModel.h"

/* Constructors */
PBModel::PBModel() = default;

PBModel::PBModel(char const *f, realtype kb1, realtype kb2, realtype kc1, realtype kc2,
                 const Grid &g, const SystemProperties &s,
                 const Fluid &cont, const Fluid &disp):
        filename(f), grid(g), sysProps(s), cont(cont), disp(disp), cvode_mem(nullptr),
        kerns(Kernels(kb1, kb2, kc1, kc2, 1, g, s, cont, disp)){
    int flag = 0;
    /* Get rows and columns of csv file and initialize M and N respectively */
    this->getRowsAndCols();

    /* Allocate memory for member variables */
    this->t      = gsl_vector_alloc(this->M);
    this->r      = gsl_vector_alloc(this->N);
    this->fv     = gsl_matrix_alloc(this->M, this->N);
    this->tau    = gsl_vector_alloc(this->M);
    this->fvSim  = gsl_matrix_calloc(this->M, this->N);

    /* These must be on another domain (xi, not r) */
    this->psi    = gsl_matrix_calloc(this->M, grid.getN());
    this->NPsi   = N_VNew_Serial(grid.getN());

    /* Set experimental data: r, t, fv and rescale */
    this->getDistributions();
    this->rescaleInitial();

    /* Set final time of experiment and update kernels */
    realtype tf = gsl_vector_get(this->t, this->M - 1);
    this->kerns.setTf(tf);

    /* Assign nondimensional time tau = t / tf */
    gsl_vector_memcpy(this->tau , this->t);
    gsl_vector_scale(this->tau, 1/this->kerns.getTf());
    this->tRequested = gsl_vector_get(this->tau, 1);

    /* Assign nondimensional initial distribution */
    this->psiN   = gsl_matrix_row(this->psi, 0);
    flag = this->preparePsi();
    if (flag != 0) perror("Failed to interpolate fv onto psi");

    /* Prepare memory for integration */
    flag = this->prepareCVMemory();
    if (flag == 1) perror("Failed to prepare ODE memory");
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

int PBModel::preparePsi(){
    realtype xN, yN;
    realtype *psiData = NV_DATA_S(this->NPsi);
    /* Allocate memory for Steffen spline on experimental psi */
    gsl_vector_view fv0 = gsl_matrix_row(this->fv, 0);
    gsl_vector *psi0 = gsl_vector_alloc(this->N);
    gsl_vector_memcpy(psi0, &fv0.vector);
    gsl_vector_scale(psi0, this->sysProps.getRm());

    /* Create temporary experimental radius / Rm */
    gsl_vector *tmp = gsl_vector_alloc(this->N);
    gsl_vector_memcpy(tmp, this->r);
    gsl_vector_scale(tmp, 1/sysProps.getRm());

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen, this->N);
    gsl_spline_init(spline, &tmp->data[0], &psi0->data[0], this->N);

    for (size_t i = 0; i < grid.getN(); i++){
        xN = gsl_vector_get(grid.getXi(), i);
        /* We are not allowed to "interpolate" outside of experimental radius */
        if (xN > gsl_vector_get(tmp, 0) && xN < gsl_vector_get(tmp, this->N-1)) {
            yN = gsl_spline_eval(spline, xN, acc);
        } else yN = 0; /* If we are outside of experimental radius, 0 our distribution */
        gsl_vector_set(&psiN.vector, i, yN);
        /* Assign initial condition to solution vector */
        psiData[i] = yN;
    }

    gsl_vector_free(psi0);
    gsl_vector_free(tmp);
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

    return 0;
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
    flag = CVodeInit(this->cvode_mem, dydt, gsl_vector_get(this->tau, 0), this->NPsi);
    if (checkFlag(&flag, "CVodeInit", 1)) return(1);

    /* Call CVodeSStolerances to specify the scalar relative tolerance
     * and scalar absolute tolerance */
    flag = CVodeSStolerances(this->cvode_mem, RTOL, ATOL);
    if (checkFlag(&flag, "CVodeSStolerances", 1)) return(1);

    /* Set the pointer to user-defined data */
    flag = CVodeSetUserData(cvode_mem, this);
    if(checkFlag(&flag, "CVodeSetUserData", 1)) return(1);

    /* Create dense SUNMatrix for use in linear solves */
    this->A = SUNDenseMatrix(this->grid.getN(), this->grid.getN());
    if(checkFlag((void *)this->A, "SUNDenseMatrix", 0)) return(1);

    /* Create dense SUNLinearSolver object for use by CVode */
    this->LS = SUNDenseLinearSolver(this->NPsi, this->A);
    if(checkFlag((void *)this->LS, "SUNDenseLinearSolver", 0)) return(1);

    /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
    flag = CVDlsSetLinearSolver(this->cvode_mem, this->LS, this->A);
    if(checkFlag(&flag, "CVDlsSetLinearSolver", 1)) return(1);
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
int PBModel::getRHS(N_Vector y, N_Vector ydot){
    size_t i = 0, j = 0;
    double integralBB; /* Stores integral value for BB */
    double integralBC; /* Stores integral value for BC */
    double integralDC; /* Stores integral value for DC */
    realtype *ydata = NV_DATA_S(y);
    realtype *ydotdata = NV_DATA_S(ydot);

    /* Fetch different grids */
    size_t Np          = this->grid.getN();
    gsl_vector *w      = grid.getW();
    gsl_vector *xi     = grid.getXi();
    gsl_matrix *xipBB  = grid.getXipBB();
    gsl_matrix *xipBC  = grid.getXipBC();
    gsl_matrix *xippBC = grid.getXippBC();

    /* Fetch different kernels */
    gsl_matrix *KBB    = kerns.getKBB();
    gsl_vector *KDB    = kerns.getKDB();
    gsl_matrix *KBC    = kerns.getKBC();
    gsl_matrix *KDC    = kerns.getKDC();

    /* Point psiN.vector to ydata */
    this->psiN.vector.data = ydata;

    /* Allocate memory for our interpolated distributions */
    gsl_matrix *psipBB  = gsl_matrix_alloc(Np, Np);
    gsl_matrix *psipBC  = gsl_matrix_alloc(Np, Np);
    gsl_matrix *psippBC = gsl_matrix_alloc(Np, Np);

    /* Interpolate from (x,y)-pairs to (xx,yy)-pairs */
    interpolatePsi(xi, &this->psiN.vector, xipBB, psipBB);
    interpolatePsi(xi, &this->psiN.vector, xipBC, psipBC);
    interpolatePsi(xi, &this->psiN.vector, xippBC, psippBC);

    /* Allocate memory for integrands (and integrals?) */
    /* IBB, IBC, IDC (and BB, DB, BC, DC?) */
    gsl_matrix *IBB = gsl_matrix_calloc(Np, Np);
    gsl_matrix *IBC = gsl_matrix_calloc(Np, Np);
    gsl_matrix *IDC = gsl_matrix_calloc(Np, Np);
    gsl_vector *B   = gsl_vector_calloc(Np);
    gsl_vector *C   = gsl_vector_calloc(Np);
    /* B = IBB*w - DB, C = IBC*w - IDC*w */
    /* ydot[i] = B[i] + C[i] */
    /* Loop over rows and columns and evaluate RHS */
    for (i = 1; i < Np; i++){
        realtype kdbi = gsl_vector_get(KDB, i);
        realtype xii = gsl_vector_get(xi, i);
        realtype psii = gsl_vector_get(&this->psiN.vector, i);
        for (j = 0; j < Np; j++){
            /* Fetch indexed variables for easier typing */
            realtype xipbbij = gsl_matrix_get(xipBB, i, j);
            realtype xipbcij = gsl_matrix_get(xipBC, i, j);
            realtype xippbcij = gsl_matrix_get(xippBC, i, j);
            realtype xij = gsl_vector_get(xi, j);
            realtype kbbij = gsl_matrix_get(KBB, i, j);
            realtype kbcij = gsl_matrix_get(KBC, i, j);
            realtype kdcij = gsl_matrix_get(KDC, i, j);
            realtype psipbbij = gsl_matrix_get(psipBB, i, j);
            realtype psipbcij = gsl_matrix_get(psipBC, i, j);
            realtype psippbcij = gsl_matrix_get(psippBC, i, j);
            realtype psij = gsl_vector_get(&this->psiN.vector, j);

            /* Set integrand birth breakage */
            realtype denomBB = SUNRpowerI(xipbbij, 3);
            if (denomBB == 0) { gsl_matrix_set(IBB, i, j, 0); }
            else {
                gsl_matrix_set(IBB, i, j,
                               kbbij * psipbbij
                               / denomBB
                               * (1 - xii)
                );
            }
            /* Set integrand birth coalescence */
            realtype denomBC1 = SUNRpowerI(xipbcij, 3);
            realtype denomBC2 = SUNRpowerI(xippbcij, 3);
            if (denomBC1 == 0 || denomBC2 == 0){ gsl_matrix_set(IBC, i, j, 0); }
            else {
                gsl_matrix_set(IBC, i, j,
                               kbcij
                               * psipbcij / denomBC1
                               * psippbcij / denomBC2
                               * (realtype) SUNRpowerI(xii / xippbcij, 2)
                               * xii / (realtype) SUNRpowerR(2.0, 1.0 / 3.0)
                );
            }
            /* Set integrand death coalescence */
            realtype denomDB = SUNRpowerI(xij, 3);
            if (denomDB == 0) { gsl_matrix_set(IDC, i, j, 0); }
            else {
                gsl_matrix_set(IDC, i, j,
                               kdcij * psij / denomDB
                );
            }
        }
        /* Now do inner product of integrand and weights to evaluate integrals */
        /* Populate B by BB - DB
         * BB[i] =  xii^3 * IBB[i, :] * w = xii^3 * ddot(IBB[i, :], w)
         */
        gsl_vector_view IBBrow = gsl_matrix_row(IBB, i);
        gsl_vector_view IBCrow = gsl_matrix_row(IBC, i);
        gsl_vector_view IDCrow = gsl_matrix_row(IDC, i);

        /* Breakage */
        gsl_blas_ddot(w, &IBBrow.vector, &integralBB);
        /* B[i] = Birth breakage[i] - Death breakage[i] */
        gsl_vector_set(B, i,
                       (realtype)SUNRpowerI(xii, 3)
                       *(realtype)integralBB
                       -kdbi*psii
        );

        /* Coalescence */
        gsl_blas_ddot(w, &IBCrow.vector, &integralBC);
        gsl_blas_ddot(w, &IDCrow.vector, &integralDC);
        /* C[i] = Birth coalescence[i] - Death coalescence[i] */
        gsl_vector_set(C, i,
                       (realtype)SUNRpowerI(xii, 3)
                       *(realtype)integralBC
                       -psii*(realtype)integralDC
        );
        ydotdata[i] = (realtype) (gsl_vector_get(B, i) + gsl_vector_get(C, i));
    }
    ydotdata[0] = 0;

    /* Free allocated memory that is only used in current scope */
    gsl_matrix_free(psipBB);
    gsl_matrix_free(psipBC);
    gsl_matrix_free(psippBC);
    gsl_matrix_free(IBB);
    gsl_matrix_free(IBC);
    gsl_matrix_free(IDC);
    gsl_vector_free(B);
    gsl_vector_free(C);
    return 0;
}

int PBModel::interpolatePsi(const gsl_vector *x, const gsl_vector *y, const gsl_matrix *xx, gsl_matrix *yy){
    /* x and y is original data, xx and yy is interpolated data */
    size_t i=0, j=0;
    realtype xN, yN;
    realtype x0   = gsl_vector_get(x, 0);
    realtype xend = gsl_vector_get(x, x->size - 1);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();

    gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen, x->size);
    gsl_spline_init(spline, x->data, y->data, x->size);
    /* Interpolate onto xx domain and get yy values */
    for (i = 0; i < xx->size1; i++) {     /* Loop over rows       */
        for (j = 0; j < xx->size2; j++) { /* Loop over columns    */
            xN = gsl_matrix_get(xx, i, j);
            if (xN > x0 && xN < xend) {
                yN = gsl_spline_eval(spline, xN, acc);
            } else yN = 0; /* If we are outside of experimental radius, 0 our distribution */
            gsl_matrix_set(yy, i, j, yN);
        }
    }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
    return 0;
}

int PBModel::interpolateFv(const gsl_vector *x, const gsl_vector *y, const gsl_vector *xx, gsl_vector *yy){
    /* x and y is original data. xx and yy is interpolated data */
    size_t i=0;
    realtype xN, yN;
    realtype x0   = gsl_vector_get(x, 0);
    realtype xend = gsl_vector_get(x, x->size -1);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();

    gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen, x->size);
    gsl_spline_init(spline, x->data, y->data, x->size);
    /* Interpolate onto xx domain and get yy values */
    for (i = 0; i < xx->size; i++){
        xN = gsl_vector_get(xx, i);
        if (xN > x0 && xN < xend){
            yN = gsl_spline_eval(spline, xN, acc);
        } else yN = 0;
        gsl_vector_set(yy, i, yN);
    }
    gsl_spline_free (spline);
    gsl_interp_accel_free(acc);
    return 0;
}

int PBModel::timeIterate() {
    int flag = CVode(this->cvode_mem, this->tRequested, this->NPsi, &(this->tout), CV_NORMAL);
    if(checkFlag(&flag, "CVode", 1)) return 1;
//    std::cout << "time requested: " << this->tRequested << ", time produced: " << this->tout << std::endl;
    return 0;
}

int PBModel::solvePBE(){
    size_t i = 0;
    int flag = 0;
    for (i = 1; i < this->M; i++){
        /* Request new return time for ODE solver */
        this->tRequested = gsl_vector_get(this->tau, i);

        /* Take one time iteration at solving the ODE */
        flag = this->timeIterate();
        if (flag == 1) { break; }

        /* Copy data of current row into psi (psiN->data points to solution data) */
        gsl_vector_view row = gsl_matrix_row(this->psi, i);
        gsl_vector_memcpy(&row.vector, &(this->psiN.vector));
    }

    /** Interpolate solution back onto experimental radial domain   **/
    /** i.e. psi(xi, tau) --> fv(r, t)                              **/
    /* Matrix fvSim will hold fv on experimental domain             */
    /* Create temporary matrix to hold fv on discretized domain     */
    gsl_matrix *tmpPsi   = gsl_matrix_alloc(this->M, this->grid.getN());
    /* Create temporary vector to hold simulated radii on [0, Rm]   */
    gsl_vector *tmpR     = gsl_vector_alloc(this->grid.getN());
    /* Copy original data */
    gsl_vector_memcpy(tmpR, this->grid.getXi());
    gsl_matrix_memcpy(tmpPsi, this->psi);
    /* Scale original data */
    gsl_vector_scale(tmpR, sysProps.getRm());
    gsl_matrix_scale(tmpPsi, 1/sysProps.getRm());

    for (i = 0; i < this->M; i++) {
        gsl_vector_view tmpPsii   = gsl_matrix_row(tmpPsi, i);
        gsl_vector_view fvSimi = gsl_matrix_row(this->fvSim, i);
        interpolateFv(tmpR, &tmpPsii.vector, this->r, &fvSimi.vector);
    }
//    this->printFvSimulated();
    gsl_matrix_free(tmpPsi);
    gsl_vector_free(tmpR);
    return 0;
}

realtype PBModel::getResidualij(size_t i, size_t j){
    /* fvSim was set in solvePBE method and should by now
     * hold simulated fv on experimental radial domain
     */
    realtype fvSimVal = gsl_matrix_get(this->fvSim, i, j);
    realtype fvExpVal = gsl_matrix_get(this->fv, i, j);
    return (fvSimVal - fvExpVal);
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
            std::cout << std::setw(8) << std::setprecision(3) << gsl_matrix_get(fv, i, j) << "\t";
        }
        std::cout << std::endl;
    }
}

void PBModel::printSizeClasses() {
    size_t i = 0;
    std::cout << "Measured size classes: " << std::endl;
    for (i = 0; i < this->N; i++){
        std::cout << std::setw(10) << std::setprecision(7) << gsl_vector_get(this->r, i);
    }
    std::cout << std::endl;
}

void PBModel::printCurrentPsi(){
    realtype *data = NV_DATA_S(this->NPsi);
    size_t i = 0;
    std::cout << "Psi for the current time iteration is: " << std::endl;
    for (i = 0; i < grid.getN(); i++){
        std::cout << std::setw(8) << std::setprecision(2) << data[i];
    }
    std::cout << std::endl;
}

void PBModel::printPsi(){
    size_t i = 0, j = 0;
    std::cout << "Nondimensionalized droplet size density distribution:" << std::endl;
    for (i = 0; i < this->M; i++){
        for (j = 0; j < grid.getN(); j++){
            std::cout << std::setw(12) << std::setprecision(3) << gsl_matrix_get(psi, i, j);
        }
        std::cout << std::endl;
    }
}

void PBModel::printFvSimulated(){
    size_t i = 0, j = 0;
    std::cout << "Droplet size density distribution fvSim(r, t):" << std::endl;
    for (i = 0; i < this->M; i++){
        for (j = 0; j < this->N; j++){
            std::cout << std::setw(12) << std::setprecision(3) << gsl_matrix_get(this->fvSim, i, j);
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


/* Exporter methods */
int PBModel::exportFvSimulatedWithExperimental(const std::string &fileName) {
    if (fileExists(fileName)) {return 1;} /* fileExists declared inline in header */
    size_t i = 0, j = 0;
    std::ofstream outfile;

    outfile.open(fileName);
    outfile << "#r,fv\n";
    for (i = 0; i < this->N; i++) {
        outfile << gsl_vector_get(this->r, i) << ",";
        for (j = 0; j < this->M; j++) {
            outfile << gsl_matrix_get(this->fvSim, j, i) << ",";
        }
        for (j = 0; j < this->M; j++){
            outfile << gsl_matrix_get(this->fv, j, i) << ",";
        }
        outfile << std::endl;
    }
    outfile.close();
    return 0;
}

int PBModel::exportPsiWithExperimental(const std::string &fileName) {
    if (fileExists(fileName)) { return 1; } /* fileExists declared inline in header */
    // TODO: Create this method by interpolating back onto measured radii and export both as in exportFvSimulatedWithExperimental
    return 0;
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
