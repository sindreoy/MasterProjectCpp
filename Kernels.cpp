//
// Created by Sindre Bakke Øyen on 07.03.2018.
//

#include <gsl/gsl_matrix.h>
#include "Kernels.h"

/* Constructors */
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

Kernels::Kernels(realtype kb1, realtype kb2, realtype kc1, realtype kc2, realtype tf,
                 const Grid &grid, const SystemProperties &sysProps,
                 const Fluid &cont, const Fluid &disp) :
    kb1(kb1), kb2(kb2), kc1(kc1), kc2(kc2), tf(tf), grid(grid),
    KBB(gsl_matrix_calloc(grid.getN(), grid.getN())),
    KBC(gsl_matrix_calloc(grid.getN(), grid.getN())),
    KDC(gsl_matrix_calloc(grid.getN(), grid.getN())),
    KDB(gsl_vector_calloc(grid.getN())){
    this->initializeKs(cont, disp, sysProps);
    this->setBreakageKernels();
    this->setCoalescenceKernels();
}

void Kernels::setBreakageKernels() {
    size_t i, j, N;
    realtype xi_i, xipBBij;
    gsl_vector *xi;
    gsl_matrix *xipBB;

    N = grid.getN();
    xi = grid.getXi();
    xipBB = grid.getXipBB();

    for (i = 1; i < N; i++){
        xi_i = gsl_vector_get(xi, i);

        /* Death breakage */
        gsl_vector_set(this->KDB, i,
                       this->k1
                       *1/SUNRpowerR(xi_i, (realtype) 2.0/3.0)
                       *SUNRexp(-this->k2/SUNRpowerR(xi_i, (realtype) 5.0/3.0)));
        for (j = 1; j < N; j++){
            xipBBij = gsl_matrix_get(xipBB, i, j);

            /* Birth breakage */
            gsl_matrix_set(this->KBB, i, j,
                           this->k1
                           *(2*1/(realtype)SUNRpowerR(xipBBij, (realtype) 2.0/3.0))
                           *SUNRexp(-this->k2/(realtype)SUNRpowerR(xipBBij, (realtype) 5.0/3.0))
                           *(2.4/(realtype)SUNRpowerI(xipBBij, 3))
                           *SUNRexp(-4.5*(realtype)SUNRpowerI(
                                   2*(realtype)SUNRpowerI(xi_i, 3) - (realtype)SUNRpowerI(xipBBij, 3), 2
                           ) /(realtype)SUNRpowerI(xipBBij, 6)
                           )
                           *3*(realtype)SUNRpowerI(xi_i, 2)
            );
        }
    }
//    /* Set first value to 0, to avoid NaN */
//    gsl_vector_set(KDB, 0, 0);
//    gsl_matrix_set(KBB, 0, 0, 0);
}

void Kernels::setCoalescenceKernels() {
    size_t i, j, N;
    realtype xi_i, xi_j, xipBCij, xippBCij;

    gsl_vector *xi      = grid.getXi();
    gsl_matrix *xipBC   = grid.getXipBC();
    gsl_matrix *xippBC  = grid.getXippBC();

    N  = grid.getN();
    for (i = 1; i < N; i++){
        xi_i = gsl_vector_get(xi, i);
        for (j = 1; j < N; j++){
            xi_j     = gsl_vector_get(xi, j);
            xipBCij  = gsl_matrix_get(xipBC, i, j);
            xippBCij = gsl_matrix_get(xippBC, i, j);
            /* Birth coalescence */
            gsl_matrix_set(this->KBC, i, j,
                           this->k3
                           *(realtype)SUNRpowerI(xipBCij+xippBCij, 2)
                           *(realtype)SUNRpowerR(
                                   (realtype)SUNRpowerR(xipBCij, (realtype) 2.0/3.0)
                                   +(realtype)SUNRpowerR(xippBCij, (realtype) 2.0/3.0),
                                   (realtype) 1.0/2.0
                           )*(realtype)SUNRexp(
                                   -this->k4
                                   *(realtype)SUNRpowerR(1/xipBCij+1/xippBCij, (realtype)-5.0/6.0)
                           )
            );

            /* Death coalescence */
            gsl_matrix_set(this->KDC, i, j,
                           this->k3
                           *(realtype)SUNRpowerI(xi_j+xi_i, 2)
                           *(realtype)SUNRpowerR(
                                   (realtype)SUNRpowerR(xi_j, (realtype)2.0/3.0)
                                   +(realtype)SUNRpowerR(xi_i, (realtype)2.0/3.0),
                                   (realtype)1.0/2.0
                           )*(realtype)SUNRexp(
                                   -this->k4
                                   *(realtype)SUNRpowerR(1/xi_j + 1/xi_i, (realtype)-5.0/6.0)
                           )
            );
        }
    }
}

/* Setter methods */
void Kernels::initializeKs(const Fluid &cont, const Fluid &disp, const SystemProperties &s){
    realtype Rm = s.getRm();
    realtype eps = s.getEps();
    realtype rhoc = cont.getRho();
    realtype rhod = disp.getRho();
    realtype sigma = disp.getSigma();
    realtype Vm = s.getVm();

    realtype R23 = SUNRpowerR(Rm, 2.0/3);
    realtype R53 = SUNRpowerR(Rm, 5.0/3);
    realtype R73 = SUNRpowerR(Rm, 7.0/3);
    realtype R56 = SUNRpowerR(Rm, 5.0/6);
    realtype e13 = SUNRpowerR(eps, 1.0/3);
    realtype e23 = SUNRpowerR(eps, 2.0/3);
    realtype t13 = SUNRpowerR(2.0, 1.0/3);
    realtype t23 = SUNRpowerR(2.0, 2.0/3);
    realtype t53 = SUNRpowerR(2.0, 5.0/3);
    realtype rho12 = SUNRpowerR(rhoc, 1.0/2);
    realtype sigma12 = SUNRpowerR(sigma, 1.0/2);

    this->k1 = tf*kb1*e13/(t23*R23)*SUNRsqrt(rhod/rhoc);
    this->k2 = kb2*sigma / (rhod*t53*e23*R53);
    this->k3 = tf/Vm*R73*4*t13*kc1*e13;
    this->k4 = kc2*R56*rho12*e13/(2*sigma12);
}

void Kernels::setTf(realtype tf) {
    this->k1 *= tf / this->tf;
    this->k3 *= tf / this->tf;
    Kernels::tf = tf;
    this->setBreakageKernels();
    this->setCoalescenceKernels();
}

void Kernels::setKb1(realtype kb1) {
    Kernels::kb1 = kb1;
}

void Kernels::setKb2(realtype kb2) {
    Kernels::kb2 = kb2;
}

void Kernels::setKc1(realtype kc1) {
    Kernels::kc1 = kc1;
}

void Kernels::setKc2(realtype kc2) {
    Kernels::kc2 = kc2;
}

void Kernels::setNewK1(realtype kb1) {
    this->k1 *= kb1 / this->kb1;
    this->setKb1(kb1);
}

void Kernels::setNewK2(realtype kb2) {
    this->k2 *= kb2 / this->kb2;
    this->setKb2(kb2);
}

void Kernels::setNewK3(realtype kc1) {
    this->k3 *= kc1 / this->kc1;
    this->setKc1(kc1);
}

void Kernels::setNewK4(realtype kc2) {
    this->k4 *= kc2 / this->kc2;
    this->setKc2(kc2);
}

void Kernels::setNewKs(realtype kb1, realtype kb2, realtype kc1, realtype kc2){
    this->setNewK1(kb1);
    this->setNewK2(kb2);
    this->setNewK3(kc1);
    this->setNewK4(kc2);
    /* Also get new kernels because of new k's */
    this->setBreakageKernels();
    this->setCoalescenceKernels();
}

/* Getter methods */
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

realtype Kernels::getK1() const {
    return k1;
}

realtype Kernels::getK2() const {
    return k2;
}

realtype Kernels::getK3() const {
    return k3;
}

realtype Kernels::getK4() const {
    return k4;
}

realtype Kernels::getKb1() const {
    return kb1;
}

realtype Kernels::getKb2() const {
    return kb2;
}

realtype Kernels::getKc1() const {
    return kc1;
}

realtype Kernels::getKc2() const {
    return kc2;
}

realtype Kernels::getTf() const {
    return tf;
}

const Grid &Kernels::getGrid() const {
    return grid;
}

/* Relational operators */
Kernels& Kernels::operator=(const Kernels &rhs){
    this->k1 = rhs.getK1();
    this->k2 = rhs.getK2();
    this->k3 = rhs.getK3();
    this->k4 = rhs.getK4();
    this->kb1 = rhs.getKb1();
    this->kb2 = rhs.getKb2();
    this->kc1 = rhs.getKc1();
    this->kc2 = rhs.getKc2();
    this->tf = rhs.getTf();
    size_t N = rhs.getGrid().getN();
    this->KBB = gsl_matrix_calloc(N, N);
    this->KBC = gsl_matrix_calloc(N, N);
    this->KDC = gsl_matrix_calloc(N, N);
    this->KDB = gsl_vector_calloc(N);
    return *this;
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
    gsl_vector_free(this->KDB);
    gsl_matrix_free(this->KBB);
    gsl_matrix_free(this->KBC);
    gsl_matrix_free(this->KDC);
}
