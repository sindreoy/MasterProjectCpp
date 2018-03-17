//
// Created by Sindre Bakke Øyen on 07.03.2018.
//

#ifndef MASTERPROJECTCPP_KERNELS_H
#define MASTERPROJECTCPP_KERNELS_H
/* User-defined header files */
#include "Grid.h"
#include "SystemProperties.h"
#include "Constants.h"

class Kernels {
private:
    gsl_matrix *KBB, *KBC, *KDC; /* KDC is symmetric */
    gsl_vector *KDB;
public:
    Kernels();
    Kernels(const Constants &consts, const Grid &grid, const SystemProperties &sysProps);
    // TODO: (Optional) Avoid double for loops and use elementwise operations
    void setBreakageKernels(const Constants &c, const Grid &g);
    void setCoalescenceKernels(const Constants &c, const Grid &g);

    gsl_matrix *getKBB() const;
    gsl_matrix *getKBC() const;
    gsl_matrix *getKDC() const;
    gsl_vector *getKDB() const;

    friend std::ostream& operator<<(std::ostream &os, const Kernels &kernels);
    ~Kernels();
};


#endif //MASTERPROJECTCPP_KERNELS_H
