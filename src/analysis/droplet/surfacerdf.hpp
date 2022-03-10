#ifndef __SURFACERDF__HPP
#define __SURFACERDF__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class SurfaceRdf:public Property{
private:
    real dsz;
    Ivec liquididx;
    Svec liquidgrps;
    real surfaceb, surfacet;
    int nliqptcls;
    long nref_p, nref_s, nref;
    real surfdens;
    Rvec boxcx, boxcy;


public:
    SurfaceRdf(){}
    SurfaceRdf(InitialSet initset);
    ~SurfaceRdf(){}


    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    virtual void normalizeResults();

};
};
#endif

