#ifndef __BRIDGEZD__HPP
#define __BRIDGEZD__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class BridgeZdensity:public Property{
private:
    real dx, xmin, xmax;
    real surfacet, surfaceb;
    Ivec liquididx;
    Svec liquidgrps;
    int nliqptcls;


    


public:
    BridgeZdensity(){}
    BridgeZdensity(InitialSet initset);
    ~BridgeZdensity();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    virtual void normalizeResults();
};
};



#endif
