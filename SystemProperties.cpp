//
// Created by Sindre Bakke Øyen on 05.03.2018.
//

#include "SystemProperties.h"

SystemProperties::SystemProperties() = default;

SystemProperties::SystemProperties(const SystemProperties &s) :
        Rm(s.getRm()), Vl(s.getVl()), Vm(s.getVm()), P(s.getP()),
        tf(s.getTf()), eps(s.getEps()) {}

SystemProperties::SystemProperties(
        realtype Rm, realtype Vl, realtype P, realtype tf, const Fluid &disp) :
        Rm(Rm), Vl(Vl), P(P), tf(tf)
{
    this->Vm = 4.0/3 * M_PI * SUNRpowerI(this->Rm, 3);
    this->eps = this->P / (disp.getRho() * this->Vl);
}

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

realtype SystemProperties::getTf() const {
    return tf;
}

realtype SystemProperties::getEps() const {
    return eps;
}

std::ostream &operator<<(std::ostream &os, const SystemProperties &properties) {
    os << "Rm: " << properties.Rm << " Vl: " << properties.Vl << " Vm: " << properties.Vm << " P: " << properties.P
       << " tf: " << properties.tf << " eps: " << properties.eps;
    return os;
}

SystemProperties::~SystemProperties(){
    std::cout << "Destroying systemproperties" << std::endl;
}
