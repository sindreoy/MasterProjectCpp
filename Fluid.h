//
// Created by Sindre Bakke Øyen on 05.03.2018.
//

#ifndef MASTERPROJECTCPP_FLUID_H
#define MASTERPROJECTCPP_FLUID_H
/* Built-in header files */
#include <iostream>                 /* Used for input/output to console */
#include <ostream>                  /* Used for overloading print operator */

/* External library header files */
#include <sundials/sundials_types.h>

class Fluid {
private:
    realtype rho, sigma, nu;
public:
    Fluid();
    Fluid(realtype rho, realtype sigma, realtype nu);

    realtype getNu() const;

    realtype getSigma() const;

    realtype getRho() const;

    friend std::ostream &operator<<(std::ostream &os, const Fluid &fluid);

    ~Fluid();
};


#endif //MASTERPROJECTCPP_FLUID_H
