#ifndef __POLYMERORIENTATION__HPP
#define __POLYMERORIENTATION__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class PolymerOrientation:public Property{
private:
    MoleculeGroup polymers;
    int num_pol;
    int num_at;
    int first_idx;
    ParticleList ref;
    std::string molname;



    


public:
    PolymerOrientation(){}
    PolymerOrientation(InitialSet initset);
    ~PolymerOrientation();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
};
};



#endif
