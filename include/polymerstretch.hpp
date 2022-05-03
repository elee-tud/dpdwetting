#ifndef __POLYMERSTRETCH__HPP
#define __POLYMERSTRETCH__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class PolymerStretch:public Property{
private:
    MoleculeGroup polymers;
    int num_pol;
    int num_at;
    int first_idx;
    ParticleList ref;



    


public:
    PolymerStretch(){}
    PolymerStretch(InitialSet initset);
    ~PolymerStretch();

    virtual void getSpecificParameters(){}
    virtual void initializeVariables();
    virtual void calculateStep(int step);
};
};



#endif
