#ifndef __BRIDGESLIPVEL__HPP
#define __BRIDGESLIPVEL__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class BridgeSlipVelocity:public Property{
private:
    real zsep;
    real surfacet, surfaceb;
    Ivec liquididx;
    Svec liquidgrps;
    int nliqptcls;
    real dxfc;

    Ivec nsptclsz;
    Ivec npptclsz;
    Ivec nlptclsz;
    Rvec xcenter, xmin, xmax;
    

    


public:
    BridgeSlipVelocity(){}
    BridgeSlipVelocity(InitialSet initset);
    ~BridgeSlipVelocity();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    virtual void normalizeResults();
};
};



#endif
