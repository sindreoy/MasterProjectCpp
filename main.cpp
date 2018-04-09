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
     */

    /*****************************************************************************************/
    /* Declaration of variables                                                              */
    /*****************************************************************************************/
    const size_t Np = 200;                      /* Number of grid points    */
    char const *fileName      = "crudeB.csv";   /* Experimental data        */
    const std::string outName = "output1.dat";   /* Output filename          */
    realtype kb1 = 1.8190e-12,                  /* Model fitted parameters  */
             kb2 = 1.8190e-9,
             kc1 = 1,
             kc2 = 1e3;

    /*****************************************************************************************/
    /* Instantiation and solution                                                            */
    /*****************************************************************************************/
    /* Helper classes */
    Grid g = Grid(Np, 0, 1, 0, 0, 2);
    Fluid disp = Fluid(0.837e3, 22.0e-3, 16.88e-3); /* Oil      */
    Fluid cont = Fluid(1.0e3, 1, 1);                /* Water    */
    SystemProperties s = SystemProperties(500.0e-6, 725.0e-6, 0.366, disp);
    PBModel m = PBModel(fileName, kb1, kb2, kc1, kc2, g, s, cont, disp);
    m.solvePBE();
    m.exportFvSimulated(outName);

    /** Evaluate residuals for t = tf (end of time horizon) **/
    std::cout << "\n\n\n\n";
    size_t M = m.getM(), N = m.getN(), i = 0;
    realtype currentRes = 0;
    for (i = 0; i < N; i++){
        currentRes = m.getResidualij(M-1, i);
        std::cout << currentRes << std::endl;
    }
    /*****************************************************************************************/
    /* Input output handling                                                                 */
    /*****************************************************************************************/
    /* We are processing a input file of format
     *                                   4 variables
                      0.000000000000000e+00 kb1
                      0.000000000000000e+00 kb2
                      0.000000000000000e+00 kc1
                      0.000000000000000e+00 kc2
                                          1 functions
                                          1 ASV_1:response_fn_1
                                          4 derivative_variables
                                          1 DVV_1:kb1
                                          2 DVV_2:kb2
                                          3 DVV_3:kc1
                                          4 DVV_4:kc2
                                          0 analysis_components
                                          1 eval_id
     * Where the first one in 1 ASV_1:response_fn_1 means we require f-value only
     * 1) Define enum var_t { param1, param2, ..., paramN }
     * 2) Create map<str, var_t> to store key-value pairs.
     * 3) Fetch parameter values from input file
     * 4) Fetch active set vector (ASV);
     *      value 1: function value, value 2: derivative value, value 4: hessian value
     * 5) Compute output response
     * 6) Write response to outfile. Handle ASV[0,1,2,...,Nparams] & (1, 2 and 4)
     * 7) fout.flush and fout.close
     */


    return 0;
}