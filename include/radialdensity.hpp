#ifndef __RADIALDENSITY__HPP
#define __RADIALDENSITY__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class RadialDensity:public Property{
private:

    Real3D com;
    Real3D ref0;
    Real3D ref1;
    ParticleGroup ptcls;

public:
    RadialDensity(){}
    RadialDensity(InitialSet initset);
    ~RadialDensity(){}


    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    virtual void normalizeResults();

};
};
#endif

