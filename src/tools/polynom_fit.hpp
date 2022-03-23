#ifndef __POLYNOM_FIT__HPP
#define __POLYNOM_FIT__HPP

#include "../../dep/Eigen/Dense"
#include "../../dep/Eigen/QR"
#include "../types.hpp"
#include <iostream>
#include <cmath>
#include <vector>

namespace dpd{


class PolynomFit{
private:
    Rvec t, v;
    Eigen::MatrixXd T, V;
    Eigen::VectorXd result;
    int order;




public:
    PolynomFit(){}
    PolynomFit(int _order, Rvec &_t, Rvec &_v);
    ~PolynomFit(){}

    void fit();
    Rvec getResults();

};
};

#endif
