#ifndef __BRIDGECLINE__HPP
#define __BRIDGECLINE__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
#include "../../tools/ellipse_fit.hpp"
#include "../../tools/polynom_fit.hpp"
namespace dpd{

class BridgeCline:public Property{
private:
    real dz;
    real dcrit;
    real zcenter;
    int numx;
    real surfacet, surfaceb;
    real dr;
    int numr;
    Ivec liquididx;
    Svec liquidgrps;
    int nliqptcls;
    Rvec topdensity;
    Rvec botdensity;
    real pilheight;
    Rvec2D ellipsepts;


    


    


public:
    BridgeCline(){}
    BridgeCline(InitialSet initset);
    ~BridgeCline();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
};
};



#endif
