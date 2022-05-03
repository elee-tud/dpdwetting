#ifndef __BRIDGEPOLCONC__HPP
#define __BRIDGEPOLCONC__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class BridgePolymerConc:public Property{
private:
    real dz, dx;
    real dcrit;
    real zcenter;
    real drtpl;
    real dxclz;
    int numz, numx;
    real surfacet, surfaceb;
    int numr;
    Ivec liquididx;
    Svec liquidgrps;
    int nliqptcls;
    Rvec topdensity;
    Rvec botdensity;
    Rvec interface;


    


public:
    BridgePolymerConc(){}
    BridgePolymerConc(InitialSet initset);
    ~BridgePolymerConc();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
};
};



#endif
