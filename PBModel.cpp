//
// Created by Sindre Bakke Øyen on 18.03.2018.
//

#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>
#include "PBModel.h"

/* Constructors */
PBModel::PBModel() = default;

PBModel::PBModel(char const *f, realtype kb1, realtype kb2, realtype kc1, realtype kc2,
                 const Grid &g, const SystemProperties &s,
                 const Fluid &cont, const Fluid &disp,
                 size_t decision):
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
    if (decision) {
        std::ifstream fin("../results/raw_logNormal4.txt");
        std::string line;
        getline(fin, line);
        size_t i = 0;
        while (getline(fin, line)){
            i++;
        }
        gsl_vector_free(this->r);
        this->r = gsl_vector_calloc(i);
        gsl_matrix_free(this->fv);
        this->fv = gsl_matrix_calloc(this->M, i);
        gsl_matrix_free(this->fv);
        this->fvSim = gsl_matrix_calloc(this->M, i);
        fin.clear();
        fin.seekg(0, fin.beg);
        getline(fin, line);
        realtype val1 = 0, val2 = 0;
        size_t j = 0;
        while (std::getline(fin, line)){
            std::stringstream linestream(line);
            linestream >> val1 >> val2;
            gsl_vector_set(r, j, val1);
            gsl_matrix_set(this->fv, 0, j, val2);
            j++;
        }
        fin.close();
        this->N = i;
        gsl_vector_view fvj0 = gsl_matrix_row(this->fv, 0);
        for (j = 1; j < this->M; j++){
            gsl_vector_view fvjj = gsl_matrix_row(this->fv, j);
            gsl_vector_memcpy(&fvjj.vector, &fvj0.vector);
        }
    }

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

    /* Prepare CVode memory */
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

int PBModel::releaseCVMemory(){
    CVodeFree(&this->cvode_mem);
    SUNMatDestroy(this->A);
    SUNLinSolFree(this->LS);
    return 0;
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

bool PBModel::checkMassBalance() {
    /* Only to be called after solvePBE method */
    size_t i = 0;
    realtype phaseFraction = 0;
    realtype phasef = 0;
    for (i = 0; i < this->M; i++){
        gsl_vector_view psii = gsl_matrix_row(this->psi, i);
        gsl_vector_view fvi  = gsl_matrix_row(this->fvSim, i);
        gsl_blas_ddot(this->grid.getW(), &psii.vector, &phaseFraction);
        gsl_blas_ddot(this->r, &fvi.vector, &phasef);
        if ((realtype) SUNRabs(phaseFraction - PHI)/PHI * 100 > 5){
            /* More than relative 5% deviation. Mass not conserved */
            return false;
        }
    }
    /* At no time the mass was not conserved --> mass was conserved, return true */
    return true;
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

    /* Open two files to write B and C to file */
    std::ofstream bin("../results/breakage.dat", std::fstream::app);
    std::ofstream cin("../results/coalescence.dat", std::fstream::app);
    bin << this->tRequested * gsl_vector_get(this->t, this->M-1) << "\t";
    cin << this->tRequested / gsl_vector_get(this->t, this->M-1) << "\t";
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

        /* Write coalescence and breakage to file */
        bin << gsl_vector_get(B, i) << "\t";
        cin << gsl_vector_get(C, i) << "\t";
    }
    ydotdata[0] = 0;
    bin << std::endl;
    cin << std::endl;
    bin.close();
    cin.close();
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
    int flag = 0;
    /* Prepare psi for CVode */
    flag = this->preparePsi();
    if (flag != 0) perror("Failed to interpolate fv onto psi");
    /* Prepare memory for integration */
    flag = this->releaseCVMemory();
    if (flag == 1) perror("Failed to release ODE memory");
    flag = this->prepareCVMemory();
    if (flag == 1) perror("Failed to prepare ODE memory");

    size_t i = 0;
    /* Write breakage and coalescence contributions to file */
    std::ofstream bin("../results/breakage.dat");
    std::ofstream cin("../results/coalescence.dat");
    bin << "xi\n";
    cin << "xi\n";
    for (i = 0; i < this->grid.getN(); i++){
        bin << gsl_vector_get(this->grid.getXi(), i) << "\t";
        cin << gsl_vector_get(this->grid.getXi(), i) << "\t";
    }
    bin << std::endl;
    cin << std::endl;
    bin.close();
    cin.close();

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

double PBModel::getModeledMean(size_t t){
    double s, rdfvsim, dr;
    s = 0;
    size_t i = 0;
    for (i = 1; i < this->N; i++){
        dr = gsl_vector_get(this->r, i) - gsl_vector_get(this->r, i-1);
        rdfvsim = gsl_vector_get(this->r, i)
                  * (gsl_matrix_get(this->fvSim, t, i-1) + gsl_matrix_get(this->fvSim, t, i));
        s += dr * rdfvsim / 2;
    }
    return (s / PHI * 1.e6);
}

double PBModel::getExperimentalMean(size_t t) {
    double s, rdfvsim, dr;
    s = 0;
    size_t i = 0;
    for (i = 1; i < this->N; i++){
        dr = gsl_vector_get(this->r, i) - gsl_vector_get(this->r, i-1);
        rdfvsim = gsl_vector_get(this->r, i)
                  * (gsl_matrix_get(this->fv, t, i-1) + gsl_matrix_get(this->fv, t, i));
        s += dr * rdfvsim / 2;
    }
    return (s / PHI * 1.e6);
}

realtype PBModel::getResidualMean(size_t t){
    /* t is time instant */
    return (this->getModeledMean(t) - this->getExperimentalMean(t));
}

double PBModel::getWeightedResidual(size_t i, size_t j, double m, double s){
    double x = gsl_vector_get(this->r, j);
    double weight = this->getWeight(x, m, s);
    return (this->getResidualij(i, j) * weight);

    return 1.0;
}

double PBModel::getWeight(double x, double m, double s) {
    return (1.0 / (1.0+exp((m-x)/s)));
}


/* Levenberg-Marquardt parameter estimation */
int PBModel::levenbergMarquardtCostFunction(const gsl_vector *x, gsl_vector *f) {
    /* Evaluates the cost function at x
     * x is the vector of parameters kb1, kb2, kc1, kc2 */
    realtype kb1 = gsl_vector_get(x, 0) / 1.e5;
    realtype kb2 = gsl_vector_get(x, 1) / 1.e4;
    realtype kc1 = gsl_vector_get(x, 2) / 1.e4;
    realtype kc2 = gsl_vector_get(x, 3) / 1.e-2;
    this->kerns.setNewKs(kb1, kb2, kc1, kc2);
    this->solvePBE();
    size_t times = f->size / this->N;
    size_t i = 0, j = 0;
    for (i = 0; i < times-1; i++){
        for (j = 0; j < this->N; j++){
            size_t idx = i * this->N + j;
            gsl_vector_set(f, idx, this->getResidualij(i, j));
        }
    }
    for (j = 0; j < this->N; j++){
        size_t idx = (times-1)*this->N + j;
        gsl_vector_set(f, idx, this->getResidualij(M-1, j));
    }
    return GSL_SUCCESS;
}

int PBModel::costFunctionMean(const gsl_vector *x, gsl_vector *f){
    double kb1 = gsl_vector_get(x, 0) / 1.e4;
    double kb2 = gsl_vector_get(x, 1) / 1.e3;
    double kc1 = gsl_vector_get(x, 2) / 1.e5;
    double kc2 = gsl_vector_get(x, 3) / 1.e-1;
    this->kerns.setNewKs(kb1, kb2, kc1, kc2);
    this->solvePBE();
    size_t i = 0;
    size_t nRes = f->size; /* Number of residual means */
    for (i = 0; i < nRes; i++){
        gsl_vector_set(f, i, this->getResidualMean(i));
    }
    return GSL_SUCCESS;
}

int PBModel::paramesterEstimationSSE() {
    const size_t Ntmin = 60;    /* Minimum number of distributions chosen   */
    const size_t Ntmax = 81;    /* Maximum number of distributions chosen   */
    size_t Nt = Ntmin;          /* Number of distributions chosen           */
    const size_t p = 4;         /* Number of parameters                     */
    do {
        size_t N = Nt * this->N;      /* Number of residuals */
        size_t n = N;

        const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
        gsl_multifit_nlinear_workspace *w;
        gsl_multifit_nlinear_fdf fdf;
        gsl_multifit_nlinear_parameters fdf_params =
                gsl_multifit_nlinear_default_parameters();
        fdf_params.h_df = 1.e-2;

        std::cout << "Number of distributions: " << Nt << std::endl;

        gsl_vector *f;  /* Function */
        gsl_matrix *J;  /* Jacobian */
        gsl_matrix *covar = gsl_matrix_alloc(p, p);

        PBModel *d = this;
        /* starting values */
        double x1_scaling = 1.e5, x2_scaling = 1.e4, x3_scaling = 1.e4, x4_scaling = 1.e-2;
        double x_init[4] = {this->kerns.getKb1() * 1.e5, this->kerns.getKb2() * 1.e4,
                            this->kerns.getKc1() * 1.e4, this->kerns.getKc2() * 1.e-2};
        gsl_vector_view x = gsl_vector_view_array(x_init, p);
        double chisq, chisq0;
        int status, info;

        const double xtol = 1e-8;
        const double gtol = 1e-8;
        const double ftol = 1.e-4;

        /* define the function to be minimized */
        fdf.f = gatewayCostSSE;
        fdf.df = NULL;      /* set to NULL for finite-difference Jacobian */
        fdf.fvv = NULL;     /* not using geodesic acceleration */
        fdf.n = n;
        fdf.p = p;
        fdf.params = d;

        /* allocate workspace with default parameters */
        w = gsl_multifit_nlinear_alloc(T, &fdf_params, n, p);

        /* initialize solver with starting point and weights */
        gsl_multifit_nlinear_init(&x.vector, &fdf, w);

        /* compute initial cost function */
        f = gsl_multifit_nlinear_residual(w);
        gsl_blas_ddot(f, f, &chisq0);

        /* solve the system with a maximum of 200 iterations */
        status = gsl_multifit_nlinear_driver(200, xtol, gtol, ftol,
                                             paramEstimationCallbackSSE, NULL, &info, w);

        /* compute covariance of best fit parameters */
        J = gsl_multifit_nlinear_jac(w);
        gsl_multifit_nlinear_covar(J, 0.0, covar);

        /* compute final cost */
        gsl_blas_ddot(f, f, &chisq);

#define FIT(i) gsl_vector_get(w->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

        time_t rawtime;
        struct tm *timeinfo;
        char buffer[80];
        time(&rawtime);
        timeinfo = localtime(&rawtime);
        strftime(buffer, sizeof(buffer), "%d-%m-%Y-%I:%M:%S", timeinfo);
        std::string str(buffer);
        std::stringstream ss;
        ss << Nt;
        std::string outputFilename = "../results/parameterEstimation/" + ss.str() + "_dists_" + str + ".dat";
        std::ofstream outfile;
        outfile.open(outputFilename);

        fprintf(stderr, "summary from method '%s/%s'\n",
                gsl_multifit_nlinear_name(w),
                gsl_multifit_nlinear_trs_name(w));
        fprintf(stderr, "number of iterations: %zu\n",
                gsl_multifit_nlinear_niter(w));
        fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
        fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
        fprintf(stderr, "reason for stopping: %s\n",
                (info == 1) ? "small step size" : "small gradient");
        fprintf(stderr, "initial |f(x)| = %f\n", sqrt(chisq0));
        fprintf(stderr, "final   |f(x)| = %f\n", sqrt(chisq));

        {
            double dof = n - p;
            double c = GSL_MAX_DBL(1, sqrt(chisq / dof));

            fprintf(stderr, "chisq/dof = %g\n", chisq / dof);

            fprintf(stderr, "kb1      = %.3g +/- %.3g\n", FIT(0), c * ERR(0));
            fprintf(stderr, "kb2      = %.3g +/- %.3g\n", FIT(1), c * ERR(1));
            fprintf(stderr, "kc1      = %.3g +/- %.3g\n", FIT(2), c * ERR(2));
            fprintf(stderr, "kc2      = %.3g +/- %.3g\n", FIT(3), c * ERR(3));

            outfile << "#kb1,kb2,kc1,kc2,kb10,kb20,kc10,kc20,chisq/dof,#dists,initial,final,iter\n";
            outfile << FIT(0) / x1_scaling << "," << FIT(1) / x2_scaling
                    << "," << FIT(2) / x3_scaling << "," << FIT(3) / x4_scaling
                    << "," << x_init[0] / x1_scaling << "," << x_init[1] / x2_scaling
                    << "," << x_init[2] / x3_scaling << "," << x_init[3] / x4_scaling
                    << "," << chisq / dof << "," << Nt
                    << "," << sqrt(chisq0) << "," << sqrt(chisq)
                    << "," << gsl_multifit_nlinear_niter(w) << "\n";
            outfile << c * ERR(0) / x1_scaling << "," << c * ERR(1) / x2_scaling
                    << "," << c * ERR(2) / x3_scaling << "," << c * ERR(3) / x4_scaling << "\n";
            outfile << fdf_params.h_df << "\n";
        }
        outfile.close();
        fprintf(stderr, "status = %s\n", gsl_strerror(status));

        gsl_multifit_nlinear_free(w);
        gsl_matrix_free(covar);
        Nt++;
    } while (Nt < Ntmax);
    return 0;
}

int PBModel::parameterEstimationMean() {
    const size_t N = this->M;
    const size_t p = 4;
    const size_t n = N;

    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_fdf fdf;
    gsl_multifit_nlinear_parameters fdf_params =
            gsl_multifit_nlinear_default_parameters();
    fdf_params.h_df = 1.e-2;

    gsl_vector *f;  /* Function */
    gsl_matrix *J;  /* Jacobian */
    gsl_matrix *covar = gsl_matrix_alloc(p, p);

    PBModel *d = this;
    /* starting values */
    double x1_scaling = 1.e4, x2_scaling = 1.e3, x3_scaling = 1.e5, x4_scaling = 1.e-1;
    double x_init[4] = {this->kerns.getKb1() * x1_scaling, this->kerns.getKb2() * x2_scaling,
                        this->kerns.getKc1() * x3_scaling, this->kerns.getKc2() * x4_scaling};
    gsl_vector_view x = gsl_vector_view_array(x_init, p);
    double chisq, chisq0;
    int status, info;

    const double xtol = 1e-8;
    const double gtol = 1e-8;
    const double ftol = 1.e-4;

    /* define the function to be minimized */
    fdf.f = gatewayCostMean;
    fdf.df = NULL;      /* set to NULL for finite-difference Jacobian */
    fdf.fvv = NULL;     /* not using geodesic acceleration */
    fdf.n = n;
    fdf.p = p;
    fdf.params = d;

    /* allocate workspace with default parameters */
    w = gsl_multifit_nlinear_alloc(T, &fdf_params, n, p);

    /* initialize solver with starting point and weights */
    gsl_multifit_nlinear_init(&x.vector, &fdf, w);

    /* compute initial cost function */
    f = gsl_multifit_nlinear_residual(w);
    gsl_blas_ddot(f, f, &chisq0);

    /* solve the system with a maximum of 200 iterations */
    status = gsl_multifit_nlinear_driver(200, xtol, gtol, ftol,
                                         paramEstimationCallbackMean, NULL, &info, w);

    /* compute covariance of best fit parameters */
    J = gsl_multifit_nlinear_jac(w);
    gsl_multifit_nlinear_covar(J, 0.0, covar);

    /* compute final cost */
    gsl_blas_ddot(f, f, &chisq);

#define FIT(i) gsl_vector_get(w->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

    time_t rawtime;
    struct tm *timeinfo;
    char buffer[80];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer, sizeof(buffer), "%d-%m-%Y-%I:%M:%S", timeinfo);
    std::string str(buffer);
    std::string outputFilename = "../results/parameterEstimation/means_all_dists" + str + ".dat";
    std::ofstream outfile;
    outfile.open(outputFilename);

    fprintf(stderr, "summary from method '%s/%s'\n",
            gsl_multifit_nlinear_name(w),
            gsl_multifit_nlinear_trs_name(w));
    fprintf(stderr, "number of iterations: %zu\n",
            gsl_multifit_nlinear_niter(w));
    fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
    fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
    fprintf(stderr, "reason for stopping: %s\n",
            (info == 1) ? "small step size" : "small gradient");
    fprintf(stderr, "initial |f(x)| = %f\n", sqrt(chisq0));
    fprintf(stderr, "final   |f(x)| = %f\n", sqrt(chisq));

    {
        double dof = n - p;
        double c = GSL_MAX_DBL(1, sqrt(chisq / dof));

        fprintf(stderr, "chisq/dof = %g\n", chisq / dof);

        fprintf(stderr, "kb1      = %.3g +/- %.3g\n", FIT(0), c * ERR(0));
        fprintf(stderr, "kb2      = %.3g +/- %.3g\n", FIT(1), c * ERR(1));
        fprintf(stderr, "kc1      = %.3g +/- %.3g\n", FIT(2), c * ERR(2));
        fprintf(stderr, "kc2      = %.3g +/- %.3g\n", FIT(3), c * ERR(3));

        outfile << "#kb1,kb2,kc1,kc2,kb10,kb20,kc10,kc20,chisq/dof,initial,final,iter\n";
        outfile << FIT(0) / x1_scaling << "," << FIT(1) / x2_scaling
                << "," << FIT(2) / x3_scaling << "," << FIT(3) / x4_scaling
                << "," << x_init[0] / x1_scaling << "," << x_init[1] / x2_scaling
                << "," << x_init[2] / x3_scaling << "," << x_init[3] / x4_scaling
                << "," << chisq / dof << ","
                << "," << sqrt(chisq0) << "," << sqrt(chisq)
                << "," << gsl_multifit_nlinear_niter(w) << "\n";
        outfile << c * ERR(0) / x1_scaling << "," << c * ERR(1) / x2_scaling
                << "," << c * ERR(2) / x3_scaling << "," << c * ERR(3) / x4_scaling << "\n";
        outfile << fdf_params.h_df << "\n";
    }
    outfile.close();
    fprintf(stderr, "status = %s\n", gsl_strerror(status));

    gsl_multifit_nlinear_free(w);
    gsl_matrix_free(covar);
    return 0;
}

///* Fletcher-Reeves constrained optimization (parameter estimation) */
//double PBModel::fletcherReevesCostFunction(const gsl_vector *v) {
//    realtype kb1 = gsl_vector_get(v, 0);
//    realtype kb2 = gsl_vector_get(v, 1);
//    realtype kc1 = gsl_vector_get(v, 2);
//    realtype kc2 = gsl_vector_get(v, 3);
//    this->kerns.setNewKs(kb1, kb2, kc1, kc2);
//    this->solvePBE();
//    double result = 0;
//    size_t i = 0, j = 0;
//    for (i = 0; i < this->M; i++){
//        for (j = 0; j < this->N; j++){
//            result += pow(this->getResidualij(i, j), 2);
//        }
//    }
//    return result;
//}
//void PBModel::fletcherReevesParamEstimation(){
//    size_t iter = 0;
//    int status;
//
//    const gsl_multimin_fdfminimizer_type *T;
//    gsl_multimin_fdfminimizer *s;
//
//    PBModel *m = this;
//
//    gsl_vector *x;
//    gsl_multimin_function_fdf func;
//    func.n = 4;
//    func.f = fletcherReevesGatewayCost;
//    func.df = NULL;
//    func.fdf = NULL;
//    func.params = m;
//
//    /* Starting point */
//    x = gsl_vector_alloc(4);
//    double x_init[4] = { this->kerns.getKb1(), this->kerns.getKb2(),
//                         this->kerns.getKc1(), this->kerns.getKc2() };
//    gsl_vector_set(x, 0, this->kerns.getKb1());
//    gsl_vector_set(x, 1, this->kerns.getKb2());
//    gsl_vector_set(x, 2, this->kerns.getKc1());
//    gsl_vector_set(x, 3, this->kerns.getKc2());
//
//    T = gsl_multimin_fdfminimizer_conjugate_fr;
//    s = gsl_multimin_fdfminimizer_alloc (T, 4);
//}

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
int PBModel::exportFvSimulatedWithExperimental() {
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer,sizeof(buffer),"%d-%m-%Y-%I:%M:%S",timeinfo);
    std::string str(buffer);
    std::string outputFilename = "../results/solutionFiles/pbe-" + str + ".dat";
    std::ofstream outfile;

    size_t i = 0, j = 0;
    outfile.open(outputFilename);
    outfile << "#r,#fv\n";
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

int PBModel::exportFv(){
    /* Create new matrix and vector to dimensionalize results */
    gsl_matrix *tmpfv = gsl_matrix_alloc(this->M, this->grid.getN());
    gsl_matrix_memcpy(tmpfv, psi);
    gsl_matrix_scale(tmpfv, 1/this->sysProps.getRm());

    gsl_vector *tmpr = gsl_vector_alloc(this->grid.getN());
    gsl_vector_memcpy(tmpr, this->grid.getXi());
    gsl_vector_scale(tmpr, this->sysProps.getRm());

    /* Export to file */
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer,sizeof(buffer),"%d-%m-%Y-%I:%M:%S",timeinfo);
    std::string str(buffer);
    std::string outputFilename = "../results/solutionFiles/pbe-" + str + ".dat";
    std::ofstream outfile;

    size_t i = 0, j = 0;
    outfile.open(outputFilename);
    outfile << "#r,#fv\n";
    for (i = 0; i < this->grid.getN(); i++) {
        outfile << gsl_vector_get(tmpr, i) << ",";
        for (j = 0; j < this->M; j++) {
            outfile << gsl_matrix_get(tmpfv, j, i) << ",";
        }
        outfile << std::endl;
    }
    outfile.close();

    /* Free temporary variables */
    gsl_matrix_free(tmpfv);
    gsl_vector_free(tmpr);
    return 0;
}

int PBModel::exportPsi() {
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer,sizeof(buffer),"%d-%m-%Y-%I:%M:%S",timeinfo);
    std::string str(buffer);
    std::string outputFilename = "../results/solutionFiles/pbe-" + str + ".dat";
    std::ofstream outfile;

    size_t i = 0, j = 0;
    outfile.open(outputFilename);
    outfile << "#xi,#psi\n";
    for (i = 0; i < this->grid.getN(); i++) {
        outfile << gsl_vector_get(this->grid.getXi(), i) << ",";
        for (j = 0; j < this->M; j++) {
            outfile << gsl_matrix_get(this->psi, j, i) << ",";
        }
        outfile << std::endl;
    }
    outfile.close();
    return 0;
}

int PBModel::exportMeans(){
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer,sizeof(buffer),"%d-%m-%Y-%I:%M:%S",timeinfo);
    std::string str(buffer);
    std::string outputFilename = "../results/solutionFiles/means-" + str + ".dat";
    std::ofstream outfile(outputFilename);
    if (!outfile.good()) return 1;

    outfile << "#modeled,#experimental" << std::endl;
    size_t t;
    for (t = 0; t < this->M; t++){
        outfile << this->getModeledMean(t) << "," << this->getExperimentalMean(t) << std::endl;
    }
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
    this->releaseCVMemory();
    NV_DATA_S(this->NPsi) = nullptr;
    N_VDestroy_Serial(this->NPsi);
}


