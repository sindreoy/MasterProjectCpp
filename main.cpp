/* Built-in header files */
#include <iostream>
#include <vector>

/* User-defined header files */
#include "PBModel.h"

int main(int argc, char **argv) {
    std::ifstream fin(argv[1]);
    if (!fin) {
        std::cerr << "\nError: failure opening " << argv[1] << std::endl;
        exit(-1);
    }
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
    /* Input handling from Dakota                                                            */
    /*****************************************************************************************/
    /* We are processing a input file of format
     *                                    4 variables
     *                 0.000000000000000e+00 kb1
     *                 0.000000000000000e+00 kb2
     *                 0.000000000000000e+00 kc1
     *                 0.000000000000000e+00 kc2
     *                                     80 functions
     *                                     1 ASV_1:least_sq_term_1
     *                                     1 ASV_2:least_sq_term_2
     *                                     1 ASV_3:least_sq_term_3
     *                                     ...
     *                                     1 ASV_80:least_sq_term_80
     *                                     4 derivative_variables
     *                                     1 DVV_1:kb1
     *                                     2 DVV_2:kb2
     *                                     3 DVV_3:kc1
     *                                     4 DVV_4:kc2
     *                                     0 analysis_components
     *                                     1 eval_id
     */
    size_t i, j, k, num_vars, num_fns, num_deriv_vars;  /* num means number of                  */
    std::string vars_text, fns_text, dvv_text;          /* Description text (2nd column above)  */

    // Get the parameter std::vector and ignore the labels
    fin >> num_vars >> vars_text;
    std::vector<double> x(num_vars);
    for (i=0; i<num_vars; i++) {
        fin >> x[i];
        fin.ignore(256, '\n');
    }

    // Get the ASV std::vector and ignore the labels
    /* Possible ASV values:
     * 1 (function value)
     * 2 (gradient)
     * 3 (function value and gradient)
     * 4 (hessian)
     * 5 (function value and hessian)
     * 6 (gradient and hessian)
     * 7 (function value, gradient and hessian)
     */
    fin >> num_fns >> fns_text;
    std::vector<int> ASV(num_fns);
    for (i=0; i<num_fns; i++) {
        fin >> ASV[i];
        fin.ignore(256, '\n');
    }

    // Get the DVV std::vector and ignore the labels
    fin >> num_deriv_vars >> dvv_text;
    std::vector<int> DVV(num_deriv_vars);
    for (i=0; i<num_deriv_vars; i++) {
        fin >> DVV[i];
        fin.ignore(256, '\n');
    }

    /*****************************************************************************************/
    /* Declaration of variables                                                              */
    /*****************************************************************************************/
    const size_t Np = 200;                      /* Number of grid points    */
    char const *fileName      = "crudeB.csv";   /* Experimental data        */
    const std::string outName = "output1.dat";  /* Output filename          */
    realtype kb1 = x[0]     ,                   /* Model fitted parameters  */
             kb2 = x[1]     ,
             kc1 = x[2]*kb1 ,
             kc2 = x[3];
    /* x[2] is actually the ratio kc1/kb1. kc1 = kb1*x[2] */

    /*****************************************************************************************/
    /* Instantiation and solution                                                            */
    /*****************************************************************************************/
    /* Helper classes */
    Grid g = Grid(Np, 0, 1, 0, 0, 2);
    Fluid disp = Fluid(0.837e3, 22.0e-3, 16.88e-3); /* Oil      */
    Fluid cont = Fluid(1.0e3, 1, 1);                /* Water    */
    SystemProperties s = SystemProperties(500.0e-6, 725.0e-6, 0.366, disp);
    /* Solution class */
    PBModel m = PBModel(fileName, kb1, kb2, kc1, kc2, g, s, cont, disp);
    m.solvePBE();
//    m.exportFvSimulatedWithExperimental(outName);

    /*****************************************************************************************/
    /* Output handling to Dakota                                                             */
    /*****************************************************************************************/
    std::ofstream fout(argv[2]);
    if (!fout) {
        std::cerr << "\nError: failure creating " << argv[2] << std::endl;
        exit(-1);
    }
    fout.precision(15); // 16 total digits
    fout.setf(std::ios::scientific);
    fout.setf(std::ios::right);

    /** Evaluate residuals for t = tf (end of time horizon) **/
    size_t M = m.getM(), N = m.getN();
    realtype currentRes = 0;
    for (i = 0; i < N; i++){
        if (ASV[i] & 1) {
            currentRes = m.getResidualij(M - 1, i);
            fout << "                     " << currentRes << " f" << i+1 << std::endl;
        }
    }
    fout.flush();
    fout.close();
    return 0;
}