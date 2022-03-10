#ifndef __JUMPFREQ__HPP
#define __JUMPFREQ__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class JumpingFrequency:public Property{
private:
    R3vec ref;
    int refstep;
    int refdt;
    real dsz, dszsqr;
    real dl, dlsqr;
    Ivec srefidx;
    Ivec prefidx;
    int nsref, npref;
    Ivec2D scummul, pcummul;
    long totnsref, totnpref;
    real surfacet, surfaceb;
    Ivec liquididx;
    Svec liquidgrps;
    int nliqptcls;



    


public:
    JumpingFrequency(){}
    JumpingFrequency(InitialSet initset);
    ~JumpingFrequency();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    void normalizeResults();
};
};



#endif
