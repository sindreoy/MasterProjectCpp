//
// Created by Sindre Bakke Øyen on 07.03.2018.
//

#ifndef MASTERPROJECTCPP_KERNELS_H
#define MASTERPROJECTCPP_KERNELS_H
/* User-defined header files */
#include "Grid.h"
#include "SystemProperties.h"

class Kernels {
private:
    gsl_matrix *KBB, *KBC, *KDC; /* KDC is symmetric */
    gsl_vector *KDB;

    realtype k1, k2, k3, k4, kb1, kb2, kc1, kc2, tf;
    Grid grid;
public:
    /* Constructors */
    Kernels();
    Kernels(const Kernels &k);
    Kernels(realtype kb1, realtype kb2, realtype kc1, realtype kc2, realtype tf,
            const Grid &grid, const SystemProperties &sysProps,
            const Fluid &cont, const Fluid &disp);

    /* Setter methods */
    void initializeKs(const Fluid &cont, const Fluid &disp, const SystemProperties &s);

    void setTf(realtype tf);

    void setKb1(realtype kb1);
    void setKb2(realtype kb2);
    void setKc1(realtype kc1);
    void setKc2(realtype kc2);

    void setNewK1(realtype kb1);
    void setNewK2(realtype kb2);
    void setNewK3(realtype kc1);
    void setNewK4(realtype kc2);

    void setNewKs(realtype kb1, realtype kb2, realtype kc1, realtype kc2);

    // TODO: (Optional) Avoid double for loops and use elementwise operations
    void setBreakageKernels();
    void setCoalescenceKernels();

    /* Getter methods */
    gsl_matrix *getKBB() const;
    gsl_matrix *getKBC() const;
    gsl_matrix *getKDC() const;
    gsl_vector *getKDB() const;
    realtype getK1() const;
    realtype getK2() const;
    realtype getK3() const;
    realtype getK4() const;
    realtype getKb1() const;
    realtype getKb2() const;
    realtype getKc1() const;
    realtype getKc2() const;
    realtype getTf() const;
    const Grid &getGrid() const;

    /* Relational operators */
    Kernels &operator=(const Kernels &rhs);

    /* Friend methods */
    friend std::ostream& operator<<(std::ostream &os, const Kernels &kernels);

    /* Destructors */
    ~Kernels();
};


#endif //MASTERPROJECTCPP_KERNELS_H
