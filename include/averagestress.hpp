#ifndef __AVERAGESTRESS__HPP
#define __AVERAGESTRESS__HPP

#include "../property.hpp"

namespace dpd{
class AverageStress:public Property{
private:


public:
    AverageStress(){}
    AverageStress(InitialSet initset);
    ~AverageStress(){}

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);

};
};

#endif
