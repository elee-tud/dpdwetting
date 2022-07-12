#ifndef __BRIDGEVELXZ__HPP
#define __BRIDGEVELXZ__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class BridgeVelocityXZ:public Property{
private:
    real surfacet, surfaceb;
    real dx, dz;
    int nx, nz;
    Ivec liquididx;
    Svec liquidgrps;
    int nliqptcls;

    Ivec nsptclsx;
    Ivec npptclsx;


    

    


public:
    BridgeVelocityXZ(){}
    BridgeVelocityXZ(InitialSet initset);
    ~BridgeVelocityXZ();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    virtual void normalizeResults();
    virtual void writeOutput();
};
};



#endif
