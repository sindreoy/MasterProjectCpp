//
// Created by Sindre Bakke Øyen on 18.03.2018.
//

#ifndef MASTERPROJECTCPP_MODEL_H
#define MASTERPROJECTCPP_MODEL_H
/****************************************************************************************/
/* Preamble                                                                             */
/****************************************************************************************/
/* Built-in header files */
#include <cmath>
#include <ostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>

/* External library header files */
#include <sundials/sundials_math.h>     /* Math functions, power etc                        */
#include <sundials/sundials_types.h>    /* Data types such as realtype                      */
#include <cvode/cvode.h>                /* prototypes for CVODE fcts., consts.              */
#include <nvector/nvector_serial.h>     /* access to serial N_Vector                        */
#include <sunmatrix/sunmatrix_dense.h>  /* access to band SUNMatrix                         */
#include <sunlinsol/sunlinsol_dense.h>  /* access to band SUNLinearSolver                   */
#include <cvode/cvode_direct.h>         /* access to CVDls interface                        */

#include <gsl/gsl_spline.h>             /* Spline interpolation from GSL                    */
#include <gsl/gsl_interp.h>             /* Interpolation header from GSL                    */
#include <gsl/gsl_multifit_nlinear.h>   /* Non-linear multifit regression                   */
#include <gsl/gsl_multimin.h>           /* Multidimensional minimization                    */
#include <gsl/gsl_matrix.h>

/* User-defined header files */
#include "Kernels.h"                    /* Contains kernels, override to use other kernels  */
#include "Grid.h"                       /* Contains Gaussian quadrature rule                */
#include "SystemProperties.h"           /* Contains data from experimental setup + fluid    */

/* Define constants for program to run */
#define RTOL                RCONST(1.0e-4)              /* Relative integration tolerance              */
#define ATOL                RCONST(1.0e-8)              /* Absolute integration tolerance              */
#define PHI                 RCONST(0.7e-2)              /* Phase fraction of oil in water              */
#define TRASHROWS           RCONST(2)                   /* Rows in csv not containing relevant data    */
#define TRASHCOLS           RCONST(9)                   /* Columns in csv not containing relevant data */
#define TRUNCATETHRESHOLD   RCONST(6)                   /* Truncate 0's in csv if many consecutive 0's */

/****************************************************************************************/
/* Class declaration                                                                    */
/****************************************************************************************/
class PBModel {
private:
    char const *filename;
    size_t M, N;                /* Only used for experimental data */

    /* Experimental data */
    gsl_matrix *fv;             /* size MxN         */
    gsl_vector *r, *t;          /* size N, size M   */

    /* Modeled data */
    realtype tout, tRequested;
    gsl_matrix *psi;            /* size Mxgrid.getN()                                       */
    gsl_vector_view psiN;       /* size grid.getN()                                         */
    gsl_vector *tau;            /* size M                                                   */
    gsl_matrix *fvSim;          /* Holds modeled fv on experimental radial domain (size MxN)*/

    /* Sundials variables for evaluating ODE */
    SUNMatrix A;
    N_Vector NPsi;              /* size grid.getN() */
    SUNLinearSolver LS;
    void *cvode_mem;

    /* Classes to help evaluate model */
    const Grid grid;
    Kernels kerns;
    const SystemProperties sysProps;
    const Fluid cont, disp;
public:
    /* Constructors */
    PBModel();
    PBModel(char const *f, realtype kb1, realtype kb2, realtype kc1, realtype kc2,
            const Grid &g, const SystemProperties &s, const Fluid &cont, const Fluid &disp,
            size_t decision);

    /* Helper methods */
    void getRowsAndCols();
    void getDistributions();
    void rescaleInitial();
    int preparePsi();
    int prepareCVMemory(); /* Allocates memory for ODE solver and prepares it for solution */
    int releaseCVMemory();
    int checkFlag(void *flagvalue, const char *funcname, int opt);
    bool checkMassBalance();


    /* Solver methods */
    int getRHS(N_Vector y, N_Vector ydot);
    static int interpolatePsi(const gsl_vector *x, const gsl_vector *y, const gsl_matrix *xx, gsl_matrix *yy);
    static int interpolateFv(const gsl_vector *x, const gsl_vector *y, const gsl_vector *xx, gsl_vector *yy);
    int timeIterate();
    int solvePBE();
    realtype getResidualij(size_t i, size_t j);
    double getModeledMean(size_t t);
    double getExperimentalMean(size_t t);
    realtype getResidualMean(size_t t);
    double getWeightedResidual(size_t i, size_t j, double m, double s);
    double getWeight(double x, double m, double s);


    /* Non-linear least squares parameter estimation */
    int levenbergMarquardtCostFunction(const gsl_vector *x, gsl_vector *f);
    int levenbergMarquardtParamEstimation();


    /* Fletcher-Reeves constrained optimization (parameter estimation) */
    double fletcherReevesCostFunction(const gsl_vector *v);
    void fletcherReevesParamEstimation();


    /* Setter methods */


    /* Getter methods */
    gsl_matrix *getFv() const;
    gsl_vector *getR() const;
    gsl_vector *getT() const;
    size_t getM() const;
    size_t getN() const;

    const Grid &getGrid() const;
    const Kernels &getKerns() const;
    const SystemProperties &getSysProps() const;
    const Fluid &getCont() const;
    const Fluid &getDisp() const;


    /* Print methods */
    void printExperimentalDistribution();
    void printSizeClasses();
    void printCurrentPsi();
    void printPsi();
    void printFvSimulated();
    void printTime();
    void printTau();
    void printDimensions();


    /* Exporter methods */
    int exportFvSimulatedWithExperimental();
    int exportFv();
    int exportPsi();
    int exportMeans();


    inline static bool fileExists(const std::string &fileName){
        struct stat buf;
        return (stat(fileName.c_str(), &buf) != -1);
    }

    /* Destructors */
    ~PBModel();
};
/* CVode trick */
inline int dydt(realtype t, N_Vector y, N_Vector ydot, void *user_data){
    PBModel *obj = static_cast<PBModel *> (user_data);
    int err = obj->getRHS(y, ydot);
    return err;
}

/* Levenberg-Marquardt trick */
inline int levenbergMarquardtGatewayCost(const gsl_vector *x, void *data, gsl_vector *f){
    PBModel *obj = static_cast<PBModel *> (data);
    int err = obj->levenbergMarquardtCostFunction(x, f);
    return err;
}
inline void levenbergMarquardtCallback(const size_t iter, void *params,
                                       const gsl_multifit_nlinear_workspace *w){
    gsl_vector *f = gsl_multifit_nlinear_residual(w);
    gsl_vector *x = gsl_multifit_nlinear_position(w);
    gsl_matrix *J = gsl_multifit_nlinear_jac(w);
    double rcond;

    /* compute reciprocal condition number of J(x) */
    gsl_multifit_nlinear_rcond(&rcond, w);

    fprintf(stderr, "iter %2zu: kb1 = %.10g, kb2 = %.10g, kc1 = %.10g, kc2, = %.10g,"
                    " cond(J) = %8.4f, |f(x)| = %.4f\n",
            iter,
            gsl_vector_get(x, 0),
            gsl_vector_get(x, 1),
            gsl_vector_get(x, 2),
            gsl_vector_get(x, 3),
            1.0 / rcond,
            gsl_blas_dnrm2(f));
//    for (size_t i = 0; i < J->size1; i++){
//        for (size_t j = 0; j < J->size2; j++){
//            fprintf(stderr, "%.2g\t", gsl_matrix_get(J, i, j));
//        }
//        fprintf(stderr, "\n");
//    }
}

/* Fletcher-Reeves trick */
inline double fletcherReevesGatewayCost(const gsl_vector *v, void *params){
    PBModel *obj = static_cast<PBModel *> (params);
    double err = obj->fletcherReevesCostFunction(v);
    return err;
}
#endif //MASTERPROJECTCPP_MODEL_H
