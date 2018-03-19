//
// Created by Sindre Bakke Øyen on 12.03.2018.
//

#include "Constants.h"

Constants::Constants() = default;

Constants::Constants(const Constants &c) : k1(c.getK1()), k2(c.getK2()),
                                           k3(c.getK3()), k4(c.getK4()),
                                           kb1(c.getKb1()), kb2(c.getKb2()),
                                           kc1(c.getKc1()), kc2(c.getKc2()) {}

Constants::Constants(
        realtype kb1, realtype kb2, realtype kc1, realtype kc2,
        SystemProperties &s, const Fluid &cont, const Fluid &disp)
        : kb1(kb1), kb2(kb2), kc1(kc1), kc2(kc2)
{
    realtype Rm = s.getRm();
    realtype eps = s.getEps();
    realtype rhoc = cont.getRho();
    realtype rhod = disp.getRho();
    realtype sigma = disp.getSigma();
    realtype tf = s.getTf();
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

void Constants::setKb1(realtype kb1) {
    Constants::kb1 = kb1;
}

void Constants::setKb2(realtype kb2) {
    Constants::kb2 = kb2;
}

void Constants::setKc1(realtype kc1) {
    Constants::kc1 = kc1;
}

void Constants::setKc2(realtype kc2) {
    Constants::kc2 = kc2;
}

void Constants::setNewK1(realtype kb1) {
    this->k1 *= kb1 / this->kb1;
    this->setKb1(kb1);
}

void Constants::setNewK2(realtype kb2) {
    this->k2 *= kb2 / this->kb2;
    this->setKb2(kb2);
}

void Constants::setNewK3(realtype kc1) {
    this->k3 *= kc1 / this->kc1;
    this->setKc1(kc1);
}

void Constants::setNewK4(realtype kc2) {
    this->k4 *= kc2 / this->kc2;
    this->setKc2(kc2);
}

void Constants::setNewKs(realtype kb1, realtype kb2, realtype kc1, realtype kc2){
    this->setNewK1(kb1);
    this->setNewK2(kb2);
    this->setNewK3(kc1);
    this->setNewK4(kc2);
}

realtype Constants::getKb1() const {
    return kb1;
}

realtype Constants::getKb2() const {
    return kb2;
}

realtype Constants::getKc1() const {
    return kc1;
}

realtype Constants::getKc2() const {
    return kc2;
}

realtype Constants::getK1() const {
    return k1;
}

realtype Constants::getK2() const {
    return k2;
}

realtype Constants::getK3() const {
    return k3;
}

realtype Constants::getK4() const {
    return k4;
}

std::ostream &operator<<(std::ostream &os, const Constants &constants) {
    os << "kb1: " << constants.kb1 << " kb2: " << constants.kb2 << " kc1: " << constants.kc1 << " kc2: "
       << constants.kc2 << " k1: " << constants.k1 << " k2: " << constants.k2 << " k3: " << constants.k3 << " k4: "
       << constants.k4;
    return os;
}

Constants::~Constants() {
    std::cout << "Destroying constants" << std::endl;
}
