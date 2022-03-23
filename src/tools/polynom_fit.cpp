#include "polynom_fit.hpp"

using namespace dpd;

PolynomFit::PolynomFit(int _order, Rvec &_t, Rvec &_v){
    t=_t;
    v=_v;
    order=_order;
    T=Eigen::MatrixXd(t.size(), order+1);
    V=Eigen::VectorXd::Map(&v.front(), v.size());

}

void PolynomFit::fit(){
    for(size_t i=0;i<t.size();++i){
        for(size_t j=0;j<order+1;++j){
            T(i,j)=pow(t.at(i), j);
        }
    }
    
    result=T.householderQr().solve(V);
    return;
}

Rvec PolynomFit::getResults(){
    Rvec coeff;
    coeff.resize(order+1);
    for(int i=0;i<order+1;i++){
        coeff[i]=result[i];
    }
    return coeff;
}

