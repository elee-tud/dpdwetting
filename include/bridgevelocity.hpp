#ifndef __BRIDGEVELOCITY__HPP
#define __BRIDGEVELOCITY__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class BridgeVelocity:public Property{
private:
    real dz;
    real surfacet, surfaceb;
    Ivec npolatz, nsolatz;

    


public:
    BridgeVelocity(){}
    BridgeVelocity(InitialSet initset);
    ~BridgeVelocity();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    virtual void normalizeResults();
};
};



#endif
