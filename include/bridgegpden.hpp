#ifndef __BRIDGEGPDEN__HPP
#define __BRIDGEGPDEN__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class BridgeGrooveParticleDensity:public Property{
private:
    real dz;
    real dcrit;
    int numz, numx;
    real surfacet, surfaceb;
    real groovet, grooveb;
    real dr;
    int numr;
    Ivec liquididx;
    Svec liquidgrps;
    int nliqptcls;
    Rvec topdensity;
    Rvec botdensity;
    Rvec interface;
    real pilheight;
    real pilwidth;
    real pilgap;
    real pilperiod;

    real pilpos0;
    real srate;
    real platevel;
    

    


public:
    BridgeGrooveParticleDensity(){}
    BridgeGrooveParticleDensity(InitialSet initset);
    ~BridgeGrooveParticleDensity();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
};
};



#endif
