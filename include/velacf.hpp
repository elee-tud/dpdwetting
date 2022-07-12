#ifndef __VELACF__HPP
#define __VELACF__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
#define MAXSOL 1000
namespace dpd{

class VelocityAutoCorrelation:public Property{
private:
    MoleculeGroup polymers;
    int num_pol;
    int num_sol;
    ParticleList ref;
    R3vec2D pol_vel;
    R3vec2D sol_vel;
    Ivec solidx;



    


public:
    VelocityAutoCorrelation(){}
    VelocityAutoCorrelation(InitialSet initset);
    ~VelocityAutoCorrelation();

    virtual void getSpecificParameters(){}
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    void calculateDynamicProperty();
};
};



#endif
