//
// Created by Sindre Bakke Øyen on 05.03.2018.
//

#ifndef MASTERPROJECTCPP_SYSTEMPROPERTIES_H
#define MASTERPROJECTCPP_SYSTEMPROPERTIES_H
/* Built-in header files */
#include <cmath>
#include <ostream>

/* External library header files */
#include <sundials/sundials_math.h>     /* Math functions, power etc            */

/* User-defined header files */
#include "Fluid.h"

class SystemProperties {
private:
    realtype Rm, Vl, Vm, P, tf, eps;
public:
    SystemProperties();

    SystemProperties(
            realtype Rm, realtype Vl, realtype P, realtype tf, const Fluid &disp);

    realtype getRm() const;
    realtype getVl() const;
    realtype getVm() const;
    realtype getP() const;
    realtype getTf() const;
    realtype getEps() const;

    friend std::ostream &operator<<(std::ostream &os, const SystemProperties &properties);

    ~SystemProperties();
};


#endif //MASTERPROJECTCPP_SYSTEMPROPERTIES_H
