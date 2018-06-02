//
// Created by Sindre Bakke Øyen on 10.05.2018.
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
    const size_t Np = 180;                      /* Number of grid points    */
    char const *fileName = "crudeB.csv";        /* Experimental data        */
    realtype kb1 = 7.e-6,   /* From testLogNormalInitialCondition.cpp */
            kb2 = 2.e-4,    /* From testLogNormalInitialCondition.cpp */
            kc1 = 1.e-4,    /* From testLogNormalInitialCondition.cpp */
            kc2 = 3.e2;     /* From testLogNormalInitialCondition.cpp */
    kb1 = 8.e-6;
    kb2 = kb2;
    kc1 = 5.e-4;
    kc2 = 4.e2;

    size_t Npts = 15;
    std::vector<gsl_matrix *> SSEvector(Npts*Npts, nullptr);    /* 10x10 matrix of matrices */
    std::vector<gsl_matrix *> wSSEvector(Npts*Npts, nullptr);   /* 10x10 matrix of matrices */
    std::vector<gsl_matrix *> meanSSEvector(Npts*Npts, nullptr);   /* 10x10 matrix of matrices */
    std::vector<PBModel *> m(Npts*Npts*Npts*Npts, nullptr);
    std::vector<realtype> kb1vec(Npts, 0);
    std::vector<realtype> kb2vec(Npts, 0);
    std::vector<realtype> kc1vec(Npts, 0);
    std::vector<realtype> kc2vec(Npts, 0);

    /* Set search space */
//    realtype kb1lb = log10(1.e-7), kb1ub = log10(1.e-4);
//    realtype kb2lb = log10(1.e-5), kb2ub = log10(1.e0);
//    realtype kc1lb = log10(1.e-7), kc1ub = log10(1.e-4);
//    realtype kc2lb = log10(1.e-2), kc2ub = log10(4.e3);
    realtype kb1lb = log10(2.15e-6), kb1ub = log10(1.e-3);
    realtype kb2lb = log10(1.29e-4), kb2ub = log10(2.15e-2);
    realtype kc1lb = log10(1.67e-5), kc1ub = log10(1.e-3);
    realtype kc2lb = log10(12.95), kc2ub = log10(4.e3);

    size_t i, j, q, r, v, w;
    for (i = 0; i < Npts; i++){
        kb1vec[i] = pow(10.0, (realtype) i / (Npts-1) * (kb1ub-kb1lb) + kb1lb);
        kb2vec[i] = pow(10.0, (realtype) i / (Npts-1) * (kb2ub-kb2lb) + kb2lb);
        kc1vec[i] = pow(10.0, (realtype) i / (Npts-1) * (kc1ub-kc1lb) + kc1lb);
        kc2vec[i] = pow(10.0, (realtype) i / (Npts-1) * (kc2ub-kc2lb) + kc2lb);
    }

    /*****************************************************************************************/
    /* Instantiation and solution                                                            */
    /*****************************************************************************************/
    /* Helper classes */
    Grid g = Grid(Np, 0, 1, 0, 0, 2);
    Fluid disp = Fluid(0.837e3, 22.0e-3, 16.88e-3); /* Oil      */
    Fluid cont = Fluid(1.0e3, 1, 1);                /* Water    */
    SystemProperties s = SystemProperties(500.0e-6, 725.0e-6, 0.366, disp);

    /* Looping on all parameters */
    size_t index_ijqr = 0, index_ij = 0, N = 0, M = 0;
    realtype res = 0, wres = 0, meanres = 0, currentRes = 0, mid = 6.e-6, slack = 1.e-6;
    for (i = 0; i < Npts; i++){
        kb1 = kb1vec[i];
        for (j = 0; j < Npts; j++){
            kb2 = kb2vec[j];
            index_ij = Npts*i + j;
            SSEvector[index_ij] = gsl_matrix_calloc(Npts, Npts);
            wSSEvector[index_ij] = gsl_matrix_calloc(Npts, Npts);
            meanSSEvector[index_ij] = gsl_matrix_calloc(Npts, Npts);
            for (q = 0; q < Npts; q++) {
                kc1 = kc1vec[q];
                for (r = 0; r < Npts; r++) {
                    kc2 = kc2vec[r];
                    index_ijqr = (size_t) (i * pow(Npts, 3) + j * pow(Npts, 2) + q * pow(Npts, 1) + r * pow(Npts, 0));
                    std::cout << index_ijqr << std::endl;
                    m[index_ijqr] = new PBModel(fileName, kb1, kb2, kc1, kc2, g, s, cont, disp, 0);
                    m[index_ijqr]->solvePBE();
                    M = m[index_ijqr]->getM();
                    N = m[index_ijqr]->getN();
                    res = 0;
                    wres = 0;
                    meanres = 0;
                    for (v = 0; v < M; v++){
                        for (w = 0; w < N; w++){
                            realtype dropSize = gsl_vector_get(m[index_ijqr]->getR(), w);
                            currentRes = pow(m[index_ijqr]->getResidualij(v, w), 2);
                            res += currentRes;
                            wres += currentRes * m[index_ijqr]->getWeight(dropSize, mid, slack);
                        }
                        meanres += pow(m[index_ijqr]->getResidualMean(v), 2);
                    }
                    if (m[index_ijqr]->checkMassBalance()){
                        gsl_matrix_set(SSEvector[index_ij], q, r, res);
                        gsl_matrix_set(wSSEvector[index_ij], q, r, wres);
                        gsl_matrix_set(meanSSEvector[index_ij], q, r, meanres);
                    }
                    else {
                        gsl_matrix_set(SSEvector[index_ij], q, r, NAN);
                        gsl_matrix_set(wSSEvector[index_ij], q, r, NAN);
                        gsl_matrix_set(meanSSEvector[index_ij], q, r, NAN);
                    }
                    delete m[index_ijqr];
                }
            }
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

    std::ofstream fbparams("../results/parameterFiles/paramsBreakage" + str + ".dat");
    std::ofstream fcparams("../results/parameterFiles/paramsCoalescence" + str + ".dat");
    std::ofstream fresiduals("../results/residualFiles/SSEallTimes" + str + ".dat");
    std::ofstream fwresiduals("../results/residualFiles/wSSEallTimes" + str + ".dat");
    std::ofstream fmeanresiduals("../results/residualFiles/meanSSEallTimes" + str + ".dat");
    fbparams << "#kb1,#kb2\n";
    for (i = 0; i < Npts; i++) {
        fbparams << kb1vec[i] << "," << kb2vec[i] << std::endl;
    }
    fcparams << "#kc1,#kc2\n";
    for (i = 0; i < Npts; i++) {
        fcparams << kc1vec[i] << "," << kc2vec[i] << std::endl;
    }
    fresiduals << "#SSE\n";
    for (i = 0; i < Npts; i++){
        for (j = 0; j < Npts; j++) {
            index_ij = i * Npts + j;
            gsl_matrix *currentMat = SSEvector[index_ij];
            for (q = 0; q < Npts; q++) {
                for (r = 0; r < Npts; r++) {
                    fresiduals << gsl_matrix_get(currentMat, q, r) << ",";
                }
                fresiduals << std::endl;
            }
            gsl_matrix_free(currentMat);
        }
    }
    fwresiduals << "#wSSE\n";
    for (i = 0; i < Npts; i++){
        for (j = 0; j < Npts; j++) {
            index_ij = i * Npts + j;
            gsl_matrix *currentMat = wSSEvector[index_ij];
            for (q = 0; q < Npts; q++) {
                for (r = 0; r < Npts; r++) {
                    fwresiduals << gsl_matrix_get(currentMat, q, r) << ",";
                }
                fwresiduals << std::endl;
            }
            gsl_matrix_free(currentMat);
        }
    }
    fmeanresiduals << "#meanSSE\n";
    for (i = 0; i < Npts; i++){
        for (j = 0; j < Npts; j++) {
            index_ij = i * Npts + j;
            gsl_matrix *currentMat = meanSSEvector[index_ij];
            for (q = 0; q < Npts; q++) {
                for (r = 0; r < Npts; r++) {
                    fmeanresiduals << gsl_matrix_get(currentMat, q, r) << ",";
                }
                fmeanresiduals << std::endl;
            }
            gsl_matrix_free(currentMat);
        }
    }

    fbparams.close();
    fcparams.close();
    fresiduals.close();
    fwresiduals.close();
    fmeanresiduals.close();
    return 0;
}