//
// Created by Sindre Bakke Øyen on 05.03.2018.
//

#include "SystemProperties.h"

/* Constructors */
SystemProperties::SystemProperties() = default;

SystemProperties::SystemProperties(const SystemProperties &s) :
        Rm(s.getRm()), Vl(s.getVl()), Vm(s.getVm()), P(s.getP()), eps(s.getEps()) {}

SystemProperties::SystemProperties(
        realtype Rm, realtype Vl, realtype P, const Fluid &disp) :
        Rm(Rm), Vl(Vl), P(P)
{
    this->Vm = 4.0/3 * M_PI * SUNRpowerI(this->Rm, 3);
    this->eps = this->P / (disp.getRho() * this->Vl);
}

/* Getter methods */
realtype SystemProperties::getRm() const {
    return Rm;
}

realtype SystemProperties::getVl() const {
    return Vl;
}

realtype SystemProperties::getVm() const {
    return Vm;
}

realtype SystemProperties::getP() const {
    return P;
}

realtype SystemProperties::getEps() const {
    return eps;
}

/* Setter methods */

/* Friend methods */
std::ostream &operator<<(std::ostream &os, const SystemProperties &properties) {
    os << "Rm: " << properties.Rm << " Vl: " << properties.Vl << " Vm: " << properties.Vm
       << " P: " << properties.P << " eps: " << properties.eps;
    return os;
}

/* Destructors */
SystemProperties::~SystemProperties(){}
