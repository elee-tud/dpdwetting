#ifndef __SURFACECOVERAGE__HPP
#define __SURFACECOVERAGE__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class SurfaceCoverage:public Property{
private:
    MoleculeGroup polymers;
    int num_pol;
    int num_at;
    int first_idx;
    ParticleList ref;

    real dz;
    real dcrit;
    ParticleGroup particles_bottom;
    Real3D com;
    real dr;
    Rvec rdensity;
    Real3D ref0, ref1;
    int numr;
    Ivec liquididx;
    Svec liquidgrps;
    int nliqptcls;


    real acrit, surfaceb, surfacet;



    


public:
    SurfaceCoverage(){}
    SurfaceCoverage(InitialSet initset);
    ~SurfaceCoverage();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
};
};



#endif
