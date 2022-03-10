#ifndef __DROPSIZE__HPP
#define __DROPSIZE__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class DropSize:public Property{
private:
    bool with_wall;
    real dz;
    real dcrit;
    std::vector<ParticleGroup*> particles_atz;
    int numz;
    real surfacet, surfaceb;
    R3vec com;
    real dr;
    int numr;
    Ivec liquididx;
    Svec liquidgrps;
    int nliqptcls;
    Rvec rdensity;
    R3vec ref0, ref1;
    Rvec radius;
    Rvec hdensity;
    real hradius;
    real pilheight;


    


public:
    DropSize(){}
    DropSize(InitialSet initset);
    ~DropSize();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
};
};



#endif
