#ifndef __BRIDGEVELX__HPP
#define __BRIDGEVELX__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class BridgeVelocityX:public Property{
private:
    real surfacet, surfaceb;
    real dz;
    Ivec liquididx;
    Svec liquidgrps;
    int nliqptcls;

    Ivec nsptclsx;
    Ivec npptclsx;
    

    


public:
    BridgeVelocityX(){}
    BridgeVelocityX(InitialSet initset);
    ~BridgeVelocityX();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    virtual void normalizeResults();
};
};



#endif
