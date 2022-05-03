#ifndef __POLYMEREEVECRELAX__HPP
#define __POLYMEREEVECRELAX__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class PolymerEevecRelax:public Property{
private:
    MoleculeGroup polymers;
    int num_pol;
    ParticleList ref;
    R3vec2D etoe;
    int n_mod;
    int mon_first, mon_end;
    std::string molname;



    


public:
    PolymerEevecRelax(){}
    PolymerEevecRelax(InitialSet initset);
    ~PolymerEevecRelax();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    void calculateDynamicProperty();
};
};



#endif
