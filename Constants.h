//
// Created by Sindre Bakke Øyen on 12.03.2018.
//

#ifndef MASTERPROJECTCPP_CONSTANTS_H
#define MASTERPROJECTCPP_CONSTANTS_H
/* User-defined header files */
#include <iostream>
#include <ostream>
#include "SystemProperties.h"

class Constants {
private:
    realtype kb1, kb2, kc1, kc2, k1, k2, k3, k4;
public:
    Constants();
    Constants(realtype kb1, realtype kb2, realtype kc1, realtype kc2,
              SystemProperties &s, const Fluid &cont, const Fluid &disp);

    void setKb1(realtype kb1);
    void setKb2(realtype kb2);
    void setKc1(realtype kc1);
    void setKc2(realtype kc2);

    void setNewK1(realtype kb1);
    void setNewK2(realtype kb2);
    void setNewK3(realtype kc1);
    void setNewK4(realtype kc2);
    void setNewKs(realtype kb1, realtype kb2, realtype kc1, realtype kc2);

    realtype getKb1() const;
    realtype getKb2() const;
    realtype getKc1() const;
    realtype getKc2() const;

    realtype getK1() const;
    realtype getK2() const;
    realtype getK3() const;
    realtype getK4() const;

    friend std::ostream &operator<<(std::ostream &os, const Constants &constants);

    ~Constants();
};


#endif //MASTERPROJECTCPP_CONSTANTS_H
