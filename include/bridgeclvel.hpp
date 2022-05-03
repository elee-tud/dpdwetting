#ifndef __BRIDGECLVEL__HPP
#define __BRIDGECLVEL__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class BridgeContLineVelocity:public Property{
private:
    real dz;
    real dcrit;
    real zcenter;
    real drtpl;
    int numz, numx;
    real surfacet, surfaceb;
    real dr;
    int numr;
    Ivec liquididx;
    Svec liquidgrps;
    int nliqptcls;
    Rvec topdensity;
    Rvec botdensity;
    Rvec interface;

    real height, shearrate, shearvel;

    Ivec prev_in_cl, curr_in_cl;

    


public:
    BridgeContLineVelocity(){}
    BridgeContLineVelocity(InitialSet initset);
    ~BridgeContLineVelocity();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
};
};



#endif
