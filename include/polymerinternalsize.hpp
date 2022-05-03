#ifndef __POLYMERINTERNALSIZE__HPP
#define __POLYMERINTERNALSIZE__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class PolymerInternalSize:public Property{
private:
    MoleculeGroup polymers;
    int num_pol;
    int num_at;
    int first_idx;
    ParticleList ref;
    std::string molname;



    


public:
    PolymerInternalSize(){}
    PolymerInternalSize(InitialSet initset);
    ~PolymerInternalSize();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
};
};



#endif
