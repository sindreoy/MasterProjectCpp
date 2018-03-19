//
// Created by Sindre Bakke Øyen on 07.03.2018.
//

#include <gsl/gsl_matrix.h>
#include "Kernels.h"

Kernels::Kernels() = default;

Kernels::Kernels(const Kernels &k){
    size_t N;
    N = k.getKBB()->size1;

    /* Allocate memory for kernels */
    this->KBB = gsl_matrix_alloc(N, N);
    this->KBC = gsl_matrix_alloc(N, N);
    this->KDC = gsl_matrix_alloc(N, N);
    this->KDB = gsl_vector_alloc(N);
    /* Copy the values (not pointers) so we don't get memory leak */
    gsl_matrix_memcpy(this->KBB, k.getKBB());
    gsl_matrix_memcpy(this->KBC, k.getKBC());
    gsl_matrix_memcpy(this->KDC, k.getKDC());
    gsl_vector_memcpy(this->KDB, k.getKDB());
}

Kernels::Kernels(const Constants &consts, const Grid &grid, const SystemProperties &sysProps):
    KBB(gsl_matrix_calloc(grid.getN(), grid.getN())), KBC(gsl_matrix_calloc(grid.getN(), grid.getN())),
    KDC(gsl_matrix_calloc(grid.getN(), grid.getN())), KDB(gsl_vector_calloc(grid.getN())){

    this->setBreakageKernels(consts, grid);
    this->setCoalescenceKernels(consts, grid);
}

void Kernels::setBreakageKernels(const Constants &c, const Grid &g) {
    size_t i, j, N;
    realtype xi_i, xipBBij;
    gsl_vector *xi;
    gsl_matrix *xipBB;

    N = g.getN();
    xi = g.getXi();
    xipBB = g.getXipBB();

    for (i = 0; i < N; i++){
        xi_i = gsl_vector_get(xi, i);

        /* Death breakage */
        gsl_vector_set(KDB, i,
                       1/SUNRpowerR(xi_i, (realtype) 2.0/3.0)
                          *SUNRexp(-c.getK2()/SUNRpowerR(xi_i, (realtype) 5.0/3.0)));
        for (j = 0; j < N; j++){
            xipBBij = gsl_matrix_get(xipBB, i, j);

            /* Birth breakage */
            gsl_matrix_set(KBB, i, j,
                           (2*1/(realtype)SUNRpowerR(xipBBij, (realtype) 2.0/3.0))
                           *SUNRexp(-c.getK2()/(realtype)SUNRpowerR(xipBBij, (realtype) 5.0/3.0))
                           *(2.4/(realtype)SUNRpowerI(xipBBij, 3))
                           *SUNRexp(-4.5*(realtype)SUNRpowerI(
                                   2*(realtype)SUNRpowerI(xi_i, 3) - (realtype)SUNRpowerI(xipBBij, 3), 2
                           ) /(realtype)SUNRpowerI(xipBBij, 6)
                           )
                           *3*(realtype)SUNRpowerI(xi_i, 2)
            );
        }
    }
    /* Set first value to 0, to avoid NaN */
    gsl_vector_set(KDB, 0, 0);
    gsl_matrix_set(KBB, 0, 0, 0);
}

void Kernels::setCoalescenceKernels(const Constants &c, const Grid &g) {
    size_t i, j, N;
    realtype xi_i, xi_j, xipBCij, xippBCij, k4;

    gsl_vector *xi      = g.getXi();
    gsl_matrix *xipBC   = g.getXipBC();
    gsl_matrix *xippBC  = g.getXippBC();

    N  = g.getN();
    k4 = c.getK4();
    for (i = 0; i < N; i++){
        xi_i = gsl_vector_get(xi, i);
        for (j = 0; j < N; j++){
            xi_j     = gsl_vector_get(xi, j);
            xipBCij  = gsl_matrix_get(xipBC, i, j);
            xippBCij = gsl_matrix_get(xippBC, i, j);
            /* Birth coalescence */
            gsl_matrix_set(KBC, i, j,
                           (realtype)SUNRpowerI(xipBCij+xippBCij, 2)
                           *(realtype)SUNRpowerR(
                                   (realtype)SUNRpowerR(xipBCij, (realtype) 2.0/3.0)
                                   +(realtype)SUNRpowerR(xippBCij, (realtype) 2.0/3.0),
                                   (realtype) 1.0/2.0
                           )*(realtype)SUNRexp(
                                   -k4
                                   *(realtype)SUNRpowerR(1/xipBCij+1/xippBCij, (realtype)-5.0/6.0)
                           )
            );

            /* Death coalescence */
            gsl_matrix_set(KDC, i, j,
                           (realtype)SUNRpowerI(xi_j+xi_i, 2)
                           *(realtype)SUNRpowerR(
                                   (realtype)SUNRpowerR(xi_j, (realtype)2.0/3.0)
                                   +(realtype)SUNRpowerR(xi_i, (realtype)2.0/3.0),
                                   (realtype)1.0/2.0
                           )*(realtype)SUNRexp(
                                   -k4
                                   *(realtype)SUNRpowerR(1/xi_j + 1/xi_i, (realtype)-5.0/6.0)
                           )
            );
        }
    }
}

gsl_matrix *Kernels::getKBB() const {
    return KBB;
}

gsl_matrix *Kernels::getKBC() const {
    return KBC;
}

gsl_matrix *Kernels::getKDC() const {
    return KDC;
}

gsl_vector *Kernels::getKDB() const {
    return KDB;
}

std::ostream& operator<<(std::ostream &os, const Kernels &kernels) {
    size_t i, j, N;
    gsl_matrix *KBB = kernels.getKBB();
    gsl_vector *KDB = kernels.getKDB();
    gsl_matrix *KBC = kernels.getKBC();
    gsl_matrix *KDC = kernels.getKDC();

    N = KDB->size;
    os << "KBB (Kernel birth breakage):\n";
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            os << std::setw(8) << std::setprecision(3) << gsl_matrix_get(KBB, i, j) << "\t";
        }
        os << "\n";
    }
    os << "\n\nKDB (Kernel death breakage):\n";
    for (i = 0; i < N; i++){
            os << std::setw(8) << std::setprecision(3) << gsl_vector_get(KDB, i) << "\t";
    }
    os << "\n\nKBC (Kernel birth coalescence):\n";
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            os << std::setw(8) << std::setprecision(3) << gsl_matrix_get(KBC, i, j) << "\t";
        }
        os << "\n";
    }
    os << "\n\nKDC (Kernel death coalescence):\n";
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            os << std::setw(8) << std::setprecision(3) << gsl_matrix_get(KDC, i, j) << "\t";
        }
        os << "\n";
    }
    return os;
}

Kernels::~Kernels() {
    std::cout << "Destroying kernels" << std::endl;
    gsl_vector_free(this->KDB);
    gsl_matrix_free(this->KBB);
    gsl_matrix_free(this->KBC);
    gsl_matrix_free(this->KDC);
}
