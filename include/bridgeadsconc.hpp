#ifndef __BRIDGEADSCONC__HPP
#define __BRIDGEADSCONC__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class BridgeAdsorptionConc:public Property{
private:
    real zcenter;
    real dr;
    int numx;
    real surfacet, surfaceb;
    Ivec liquididx;
    Svec liquidgrps;
    int nliqptcls;
    Ivec npolbot, nsolbot, npoltop, nsoltop;


    


public:
    BridgeAdsorptionConc(){}
    BridgeAdsorptionConc(InitialSet initset);
    ~BridgeAdsorptionConc();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    virtual void normalizeResults();
};
};



#endif
