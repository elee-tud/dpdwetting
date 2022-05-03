#ifndef __BRIDGEJUMPFREQ__HPP
#define __BRIDGEJUMPFREQ__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
#define STATIC 0
#define PARDIFF 1
#define PERDIFF 2
#define OUTINTER 3
namespace dpd{

class BridgeJumpingFrequency:public Property{
private:
    R3vec ref;
    int refstep;
    int refdt;
    real dsz;
    real dl, dlsqr;
    real wdisp, shearrate;
    Ivec srefidx_advtop, srefidx_rectop;
    Ivec srefidx_advbot, srefidx_recbot;
    Ivec prefidx_advtop, prefidx_rectop;
    Ivec prefidx_advbot, prefidx_recbot;
    int nsref_advtop, nsref_rectop;
    int nsref_advbot, nsref_recbot;
    int npref_advtop, npref_rectop;
    int npref_advbot, npref_recbot;
    Ivec2D scummul_rectop, pcummul_rectop;
    Ivec2D scummul_advtop, pcummul_advtop;
    Ivec2D scummul_recbot, pcummul_recbot;
    Ivec2D scummul_advbot, pcummul_advbot;
    long totnsref_adv, totnpref_adv;
    long totnsref_rec, totnpref_rec;
    real surfacet, surfaceb;
    Ivec liquididx;
    Svec liquidgrps;
    int nliqptcls;
    real dcl;


    int numx;
    real dx, dcrit;
    Rvec topdensity, botdensity;
    real advinter_top, advinter_bot;
    real recinter_top, recinter_bot;



    


public:
    BridgeJumpingFrequency(){}
    BridgeJumpingFrequency(InitialSet initset);
    ~BridgeJumpingFrequency();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    void normalizeResults();
    void calculateInterfacePosition();
};
};



#endif
