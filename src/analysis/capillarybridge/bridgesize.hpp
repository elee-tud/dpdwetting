#ifndef __BRIDGESIZE__HPP
#define __BRIDGESIZE__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class BridgeSize:public Property{
private:
    real dz;
    real dcrit;
    real zcenter;
    int numz, numx;
    real surfacet, surfaceb;
    real dr;
    int numr;
    Ivec liquididx;
    Svec liquidgrps;
    int nliqptcls;
    Rvec2D topdensity;
    Rvec2D botdensity;
    Rvec2D interface;
    real pilheight;


    


public:
    BridgeSize(){}
    BridgeSize(InitialSet initset);
    ~BridgeSize();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
};
};



#endif
