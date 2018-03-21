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
#include <cstdio>
#include <cstdlib>
#include <cstring>

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
    /* TODO: Create gsl_vector_view psiN = gsl_matrix_row(psi, timeIndex)
     *       Then NPsi can be that specific time instant --> Solve that time instant
     *       getRHS only has to return NPsi, so no problems with column-major <--> row-major :D
     *       Set NPsi->data to nullptr before freeing object
     */
    realtype tProduced, tRequested;
    gsl_matrix *psi;
    gsl_vector_view psiN;
    gsl_vector *tau;
    N_Vector NPsi;
    SUNMatrix A;
    SUNLinearSolver LS;
    void *cvode_mem;

    Grid grid;
    Kernels kerns;
    SystemProperties sysProps;
    Fluid cont, disp;
public:
    /* Constructors */
    PBModel();
    PBModel(char const *f, realtype kb1, realtype kb2, realtype kc1, realtype kc2,
            const Grid &g, const SystemProperties &s, const Fluid &cont, const Fluid &disp);

    /* Helper methods */
    void getRowsAndCols();
    void getDistributions();
    void rescaleInitial();
    int prepareCVMemory(); /* Allocates memory for ODE solver and prepares it for solution */
    int checkFlag(void *flagvalue, const char *funcname, int opt);

    /* Solver methods */
    static int getRHS(realtype t, N_Vector y, N_Vector ydot, void *user_data);
    static int timeIterate();
    static int solvePBE();

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
    void printCurrentPsi();
    void printPsi();
    void printTime();
    void printTau();
    void printSizeClasses();
    void printDimensions();

    /* Destructors */
    ~PBModel();
};

#endif //MASTERPROJECTCPP_MODEL_H
