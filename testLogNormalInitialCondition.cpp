//
// Created by Sindre Bakke Øyen on 18.04.2018.
//

#include <iostream>
#include <chrono>
#include <ctime>

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
    const size_t Np = 400;                      /* Number of grid points    */
    char const *fileName = "crudeB.csv";        /* Experimental data        */
//    realtype kb1 = 2.6426477281e-07,          /* Model fitted parameters  */
//            kb2 = 0,
//            kc1 = 1.0138643992e+00*kb1,
//            kc2 = 1.6791362386e-01;
    realtype kb1 = 4.e-6*0,
            kb2 = 1.e-3,
            kc1 = 1.e-4,
            kc2 = 3.e2;
    /*****************************************************************************************/
    /* Instantiation and solution                                                            */
    /*****************************************************************************************/
    /* Helper classes */
    Grid g = Grid(Np, 0, 1, 0, 0, 2);
    Fluid disp = Fluid(0.837e3, 22.0e-3, 16.88e-3); /* Oil      */
    Fluid cont = Fluid(1.0e3, 1, 1);                /* Water    */
    SystemProperties s = SystemProperties(500.0e-6, 725.0e-6, 0.366, disp);
    PBModel m = PBModel(fileName, kb1, kb2, kc1, kc2, g, s, cont, disp, 1);
    m.solvePBE();
    m.exportFv();
    return 0;
}