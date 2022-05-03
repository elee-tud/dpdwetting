#ifndef __DROPVELOCITY__HPP
#define __DROPVELOCITY__HPP
#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{
class DropVelocity:public Property{
private:
    real dz;
    std::vector<ParticleGroup*> particles_atz;
    int numz;
    real surfacet, surfaceb;
    R3vec com;
    R3vec ref0, ref1;


public:
    DropVelocity(){}
    DropVelocity(InitialSet initset);
    ~DropVelocity(){}

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);


};
};
#endif
