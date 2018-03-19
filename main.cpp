#include <iostream>

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
     *
     */
    /*****************************************************************************************/
    /* Declaration of variables                                                              */
    /*****************************************************************************************/
    const size_t Np = 10;
    char const *fileName = "crudeB.csv";

    /*****************************************************************************************/
    /* Instantiation and solution                                                            */
    /*****************************************************************************************/
    Grid g = Grid(Np, 0, 1, 0, 0, 2);
    Fluid disp = Fluid(0.837e3, 22.0e-3, 16.88e-3);
    Fluid cont = Fluid(1.0e3, 1, 1);
    SystemProperties s = SystemProperties(500.0e-6, 725.0e-6, 0.366, 2965, disp);
    Constants c = Constants(1.8190e-12,1.819e-9,1,1.0e3,s, cont, disp);
    Kernels k = Kernels(c,g,s);
    PBModel m = PBModel(fileName, g, k, c, s);
    m.printDimensions();
    m.printDistribution();
    m.printTime();
    m.printSizeClasses();
    // TODO: Create the model to solve
    // TODO: Remove contents of cmake-build-debug directory from github
    return 0;
}