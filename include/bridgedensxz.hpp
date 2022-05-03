#ifndef __BRIDGEDENSXZ__HPP
#define __BRIDGEDENSXZ__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class BridgeDensityXZ:public Property{
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
    BridgeDensityXZ(){}
    BridgeDensityXZ(InitialSet initset);
    ~BridgeDensityXZ();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    virtual void normalizeResults();
    virtual void writeOutput();
};
};



#endif
