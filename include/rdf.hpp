#ifndef __RDF__HPP
#define __RDF__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
#define MAXRDFAVG 100
namespace dpd{

class RadialDistributionFunction:public Property{
private:
    real maxr;
    int navg;



    


public:
    RadialDistributionFunction(){}
    RadialDistributionFunction(InitialSet initset);
    ~RadialDistributionFunction(){}

    virtual void getSpecificParameters(){}
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    virtual void normalizeResults();
};
};



#endif
