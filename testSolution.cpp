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
    kb1 = 4.6416e-05, kb2 = 4.6416e-04, kc1 = 1.e-3, kc2 = 227.592;
    kb1*=1; kb2*=1; kc1*= 1; kc2*= 1;
    /*****************************************************************************************/
    /* Instantiation and solution                                                            */
    /*****************************************************************************************/
    /* Helper classes */
    Grid g = Grid(Np, 0, 1, 0, 0, 2);
    Fluid disp = Fluid(0.837e3, 22.0e-3, 16.88e-3); /* Oil      */
    Fluid cont = Fluid(1.0e3, 1, 1);                /* Water    */
    SystemProperties s = SystemProperties(500.0e-6, 725.0e-6, 0.366, disp);
    PBModel m = PBModel(fileName, kb1, kb2, kc1, kc2, g, s, cont, disp, 0);
    m.levenbergMarquardtParamEstimation();
//    m.solvePBE();
//    m.exportFvSimulatedWithExperimental();
    return 0;
}