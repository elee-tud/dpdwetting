#ifndef __BRIDGEINTERFACE__HPP
#define __BRIDGEINTERFACE__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class BridgeInterface:public Property{
private:
    real dz, dy, dr, dv;
    real dcrit;
    int numz, numy, numx;
    real surfacet, surfaceb;
    Ivec liquididx;
    Svec liquidgrps;
    int nliqptcls;
    Rvec3D density;
    Rvec3D interface;
    real pilheight;
    std::ofstream outstream;


    


public:
    BridgeInterface(){}
    BridgeInterface(InitialSet initset);
    ~BridgeInterface();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    virtual void writeOutput(){}
    virtual void writeResultsForNow(int step);
};
};



#endif
