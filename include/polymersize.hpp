#ifndef __POLYMERSIZE__HPP
#define __POLYMERSIZE__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class PolymerSize:public Property{
private:
    MoleculeGroup polymers;
    int num_pol;
    int num_at;
    int first_idx;
    ParticleList ref;
    std::string molname;



    


public:
    PolymerSize(){}
    PolymerSize(InitialSet initset);
    ~PolymerSize();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
};
};



#endif
