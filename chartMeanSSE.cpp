//
// Created by Sindre Bakke Øyen on 30.04.2018.
//

#include <iostream>
#include <chrono>
#include <ctime>
#include <vector>

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
    realtype kb1 = 7.e-6,   /* From testLogNormalInitialCondition.cpp */
            kb2 = 2.e-4,    /* From testLogNormalInitialCondition.cpp */
            kc1 = 1.e-4,    /* From testLogNormalInitialCondition.cpp */
            kc2 = 3.e2;     /* From testLogNormalInitialCondition.cpp */
    kb1 = 8.e-6;
    kb2 = kb2;
    kc1 = 5.e-4;
    kc2 = 4.e2;

    size_t Npts = 50;
    std::vector<gsl_matrix *> residualVector(Npts*Npts, nullptr); /* 19x19 matrix of matrices */
    std::vector<PBModel *> m(Npts*Npts, nullptr);
    std::vector<realtype> k1vec(Npts, 0);
    std::vector<realtype> k2vec(Npts, 0);
    realtype k1lb = log10(5.e-5), k1ub = log10(5.e-3);
    realtype k2lb = log10(40), k2ub = log10(4000);
    size_t i, j, q;
    for (i = 0; i < Npts; i++){
        k1vec[i] = pow(10.0, (realtype) i / (Npts-1) * (k1ub-k1lb) + k1lb);
        k2vec[i] = pow(10.0, (realtype) i / (Npts-1) * (k2ub-k2lb) + k2lb);
    }
    realtype k1, k2;

    /*****************************************************************************************/
    /* Instantiation and solution                                                            */
    /*****************************************************************************************/
    /* Helper classes */
    Grid g = Grid(Np, 0, 1, 0, 0, 2);
    Fluid disp = Fluid(0.837e3, 22.0e-3, 16.88e-3); /* Oil      */
    Fluid cont = Fluid(1.0e3, 1, 1);                /* Water    */
    SystemProperties s = SystemProperties(500.0e-6, 725.0e-6, 0.366, disp);
    realtype summation = 0;
    size_t index_ij = 0, N = 0, M = 0;
    for (i = 0; i < Npts; i++){
        k1 = k1vec[i];
        for (j = 0; j < Npts; j++){
            index_ij = i*Npts + j;
            k2 = k2vec[j];
            m[index_ij] = new PBModel(fileName, kb1, kb2, k1, k2, g, s, cont, disp, 0);
            m[index_ij]->solvePBE();
            M = m[index_ij]->getM();
            N = m[index_ij]->getN();
            residualVector[index_ij] = gsl_matrix_alloc(M, 1);
            for (q = 0; q < M; q++) {
                gsl_matrix_set(
                        residualVector[index_ij], q, 0,
                        m[index_ij]->getResidualMean(q)
                );
            }
            delete m[index_ij];
        }
    }

    /* Write to file */
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer,sizeof(buffer),"%d-%m-%Y-%I:%M:%S",timeinfo);
    std::string str(buffer);

//    std::ofstream fparams("../results/paramsBreakage.dat");
    std::ofstream fparams("../results/parameterFiles/paramsCoalescence" + str + ".dat");
    std::ofstream fresiduals("../results/residualFiles/residualsMeans" + str + ".dat");
//    std::ofstream fresiduals("../results/residualsNoDynamics.dat");
    fparams << "#k1,#k2\n";
    for (i = 0; i < Npts; i++) {
        fparams << k1vec[i] << "," << k2vec[i] << std::endl;
    }
    fresiduals << "#residuals\n";
    for (i = 0; i < Npts; i++){
        for (j = 0; j < Npts; j++) {
            index_ij = i * Npts + j;
            gsl_matrix *currentMat = residualVector[index_ij];
            for (q = 0; q < M; q++) {
                fresiduals << gsl_matrix_get(currentMat, q, 0) << ",";
            }
            fresiduals << std::endl;
        }
    }
    fparams.close();
    fresiduals.close();
    return 0;
}