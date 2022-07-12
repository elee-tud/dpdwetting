#ifndef __POLYMERBONDLENGTH__HPP
#define __POLYMERBONDLENGTH__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class PolymerBondlength:public Property{
private:
    MoleculeGroup polymers;
    int num_pol;
    int num_at;
    int first_idx;
    ParticleList ref;
    real maxbl;
    std::string molname;



    


public:
    PolymerBondlength(){}
    PolymerBondlength(InitialSet initset);
    ~PolymerBondlength();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
};
};



#endif
