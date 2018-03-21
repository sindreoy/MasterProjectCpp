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
    realtype Rm, Vl, Vm, P, eps;
public:
    /* Constructors */
    SystemProperties();
    SystemProperties(const SystemProperties &s);
    SystemProperties(
            realtype Rm, realtype Vl, realtype P, const Fluid &disp);

    /* Getter methods */
    realtype getRm() const;
    realtype getVl() const;
    realtype getVm() const;
    realtype getP() const;
    realtype getEps() const;

    friend std::ostream &operator<<(std::ostream &os, const SystemProperties &properties);

    /* Destructors */
    ~SystemProperties();
};


#endif //MASTERPROJECTCPP_SYSTEMPROPERTIES_H
