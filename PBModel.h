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

/* External library header files */
#include <sundials/sundials_math.h>     /* Math functions, power etc                    */
#include <sundials/sundials_types.h>    /* Datatypes such as realtype                   */
#include <cvode/cvode.h>                /* prototypes for CVODE fcts., consts.          */
#include <nvector/nvector_serial.h>     /* access to serial N_Vector                    */
#include <sunmatrix/sunmatrix_dense.h>  /* access to band SUNMatrix                     */
#include <sunlinsol/sunlinsol_dense.h>  /* access to band SUNLinearSolver               */
#include <cvode/cvode_direct.h>         /* access to CVDls interface                    */

/* User-defined header files */
#include "Kernels.h"                    /* Contains kernels, override to use other kerns*/
#include "Grid.h"                       /* Contains Gaussian quadrature rule            */
#include "SystemProperties.h"           /* Contains data from experimental setup + fluid*/
#include "Constants.h"                  /* Contains model fitted params and constants   */

/* Define constants for program to run */
#define RTOL RCONST(1.0e-4)             /* Relative integration tolerance               */
#define ATOL RCONST(1.0e-8)             /* Absolute integration tolerance               */
#define PHI  RCONST(0.7e-2)             /* Phase fraction of oil in water               */
#define TRASHROWS RCONST(2)             /* Rows in csv not containing relevant data     */
#define TRASHCOLS RCONST(9)             /* Columns in csv not containing relevant data  */
#define TRUNCATETHRESHOLD RCONST(6)     /* Truncate 0's in csv if many consecutive 0's  */

/****************************************************************************************/
/* Class declaration                                                                    */
/****************************************************************************************/
class PBModel {
private:
    char const *filename;
    size_t M, N;
    /* Experimental data */
    gsl_matrix *fv;
    gsl_vector *r, *t;

    /* Modeled data */
    /* TODO: Create gsl_vector_view psiN = gsl_matrix_row(gslPsi, timeIndex)
     *       Then NPsi can be that specific time instant --> Solve that time instant
     *       getRHS only has to return NPsi, so no problems with column-major <--> row-major :D
     */
    gsl_matrix *gslPsi;
    N_Vector NPsi, time;
    SUNMatrix A, sunPsi;
    SUNLinearSolver LS;
    void *cvode_mem;

    Grid grid;
    Kernels kerns;
    Constants consts;
    SystemProperties sysProps;
public:
    PBModel();
    PBModel(char const *f, const Grid &g, const Kernels &k,
          const Constants &c, const SystemProperties &s);

    void getRowsAndCols(size_t &rows, size_t &cols);
    void getDistributions(size_t &rows, size_t &cols);
    void rescaleInitial();

    // TODO: Create conversions between different datatypes in Sundials and GSL
//    void gslVec2NVec(N_Vector &nv, const gsl_matrix &gv);
//    void NVec2GslVec(gsl_vector &gv, const N_Vector &nv);
//    void getRHS();
//    void timeIterate();

    gsl_matrix *getFv() const;
    gsl_vector *getR() const;
    gsl_vector *getT() const;
    size_t getM() const;
    size_t getN() const;

    void printDistribution();
    void printTime();
    void printSizeClasses();
    void printDimensions();
    ~PBModel();
};

#endif //MASTERPROJECTCPP_MODEL_H
