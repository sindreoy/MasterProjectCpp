//
// Created by Sindre Bakke Øyen on 05.03.2018.
//

#include "Fluid.h"

Fluid::Fluid() : rho(0), sigma(0), nu(0) {}

Fluid::Fluid(const Fluid &f) : rho(f.getRho()), sigma(f.getSigma()), nu(f.getNu()){}

Fluid::Fluid(realtype rho, realtype sigma, realtype nu) : rho(rho), sigma(sigma), nu(nu) {}

realtype Fluid::getRho() const {
    return rho;
}

realtype Fluid::getSigma() const {
    return sigma;
}

realtype Fluid::getNu() const {
    return nu;
}

std::ostream &operator<<(std::ostream &os, const Fluid &fluid) {
    os << "rho: " << fluid.rho << " sigma: " << fluid.sigma << " nu: " << fluid.nu;
    return os;
}

Fluid::~Fluid(){
    std::cout << "Destroying fluid" << std::endl;
}
