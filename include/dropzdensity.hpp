#ifndef __DROPZDENSITY__HPP
#define __DROPZDENSITY__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class DropZDensity:public Property{
private:
    real drcenter;
    real drsqr;
    Ivec liquididx;
    Svec liquidgrps;
    real surfaceb, surfacet;
    Rvec center;
    int nliqptcls;


public:
    DropZDensity(){}
    DropZDensity(InitialSet initset);
    ~DropZDensity(){}


    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    virtual void normalizeResults();

};
};
#endif

