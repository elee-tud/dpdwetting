#ifndef __POLYMERADSORPTION__HPP
#define __POLYMERADSORPTION__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class PolymerAdsorption:public Property{
private:
    real acrit;
    real dcrit;
    real surfacet, surfaceb;
    MoleculeGroup polymers;
    int nbond_types;
    int num_pol;
    int num_at;
    int first_idx;
    Ivec ads;
    real bondl;



    


public:
    PolymerAdsorption(){}
    PolymerAdsorption(InitialSet initset);
    ~PolymerAdsorption();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);

    real findStretchingFactor(Rvec zcoords);
};
};



#endif
