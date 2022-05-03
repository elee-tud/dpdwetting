#ifndef __DROPSPHERICALSTRESS__HPP
#define __DROPSPHERICALSTRESS__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
#include "../../force.hpp"
#include "../../localdensity.hpp"
namespace dpd{
class DropSphericalStress:public Property{
private:
    Real3D com;
    Real3D ref0;
    Real3D ref1;
    ParticleGroup ptcls;
    Decomposition *decomp;
    Configuration *config;
    Interactions *inter;
    PeriodicBoundary pbc;
    Force* force;
    LocalDensity *localdensity;
    real** press;

    

public:
    DropSphericalStress(){}
    DropSphericalStress(InitialSet initset);
    ~DropSphericalStress(){}
    
    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    virtual void normalizeResults();

};
};

#endif
                        
