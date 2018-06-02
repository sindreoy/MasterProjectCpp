#include <iostream>
#include <chrono>
#include <ctime>
#include <gsl/gsl_vector_double.h>

/* User-defined header files */
#include "PBModel.h"

int main() {
    /*****************************************************************************************/
    /* Description of program                                                                */
    /*****************************************************************************************/
    /* This is the main program that solves the Population Balance Equation (PBE)
     * We have some classes to help us solve the model:
     *  - Grid              :: Contains all variables needed for Gaussian quadrature rule
     *  - Fluid             :: Contains density, surface tension and viscosity for a fluid
     *  - SystemProperties  :: Contains variables such as maximum radius, volume of tanks etc
     *  - Constants         :: Contains parameters such as k1,k2,k3,k4,kb1,kb2,kc1,kc2
     *  - Kernels           :: Contains kernels for breakage (KBB,kDB) and coalescence (KBC,KDC)
     *  - PBModel           :: Solves the entire model by the use of an ODE solver:
     *                         Utilizes all above classes
     */

    /*****************************************************************************************/
    /* Declaration of variables                                                              */
    /*****************************************************************************************/
    const size_t Np = 200;                      /* Number of grid points    */
    char const *fileName = "crudeB.csv";        /* Experimental data        */
    realtype kb1 = 8.e-6,
            kb2 = 2.e-4,
            kc1 = 1.e-4,
            kc2 = 3.e2;
//    /* Guess for parameter estimating residuals */
    kb1 = 2.99e-5, kb2 = 5.556421e-4, kc1 = 5.57417e-4, kc2 = 227.596;
//    kb1*=1; kb2*=1; kc1*= 1; kc2*= 1;

    /* Guess for parameter estimating weighted residuals */
//    kb1 = 7.1905e-05, kb2 = 8.01876e-4, kc1 = 7.46537e-4, kc2 = 151.12;

    /* Guess for parameter estimating means */
//    kb1 = 3.59381e-4, kb2 = 5.99484e-3, kc1 = 1.29155e-4, kc2 = 227.592;




    /* Solve with parameters from optimizing residuals */
////    kb1 = 3.99121e-05, kb2 = 0.000505857, kc1 = 0.00116999, kc2 = 312.743;  /* Coarse   */
//    kb1 = 2.92088e-05, kb2 = 0.000601003, kc1 = 0.000573217, kc2 = 251.809; /* Refined  */

    /* Solve with parameters from optimizing weighted residuals */
    kb1 = 4.22981e-05, kb2 = 0.000631455, kc1 = 0.000881273, kc2 = 262.607; /* Refined, coarse are same as SSE */

    /* Solve with parameters from optimizing means */
////    kb1 = 0.00035419, kb2=0.00611499, kc1=0.000122041, kc2=232.589;     /* Coarse   */
//    kb1=0.000248857,kb2=	0.00496479,kc1=	6.58189e-05,kc2=31.9293;    /* Refined  */
    /*****************************************************************************************/
    /* Instantiation and solution                                                            */
    /*****************************************************************************************/
    /* Helper classes */
    Grid g = Grid(Np, 0, 1, 0, 0, 2);
    Fluid disp = Fluid(0.837e3, 22.0e-3, 16.88e-3); /* Oil      */
    Fluid cont = Fluid(1.0e3, 1, 1);                /* Water    */
    SystemProperties s = SystemProperties(500.0e-6, 725.0e-6, 0.366, disp);
    PBModel m = PBModel(fileName, kb1, kb2, kc1, kc2, g, s, cont, disp, 0);
//    m.paramesterEstimationSSE();
//    m.parameterEstimationMean();
    m.solvePBE();
    m.exportFvSimulatedWithExperimental();
//    m.exportMeans();
    return 0;
}