#ifndef __POLYMERDEPOSITORIENT__HPP
#define __POLYMERDEPOSITORIENT__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class PolymerDepositOrient:public Property{
private:
    MoleculeGroup polymers;
    int num_pol;
    Real3D ref0, ref1;
    real dz;
    ParticleGroup liquid;
    Ivec liquididx;
    Svec liquidgrps;
    int nliqptcls;
    real surfaceb;




    


public:
    PolymerDepositOrient(){}
    PolymerDepositOrient(InitialSet initset);
    ~PolymerDepositOrient();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
};
};



#endif
