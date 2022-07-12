#ifndef __BRIDGEALDENS__HPP
#define __BRIDGEALDENS__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class BridgeAdsorblayerDensity:public Property{
private:
    real dz, zcenter;
    int numx;
    real surfacet, surfaceb;
    Ivec liquididx;
    Svec liquidgrps;
    int nliqptcls;


    


public:
    BridgeAdsorblayerDensity(){}
    BridgeAdsorblayerDensity(InitialSet initset);
    ~BridgeAdsorblayerDensity();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    virtual void normalizeResults();
};
};



#endif
