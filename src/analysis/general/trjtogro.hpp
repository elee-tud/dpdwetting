#ifndef __TRJTOGRO__HPP
#define __TRJTOGRO__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class TrjToGro:public Property{
private:
    bool dumpfrozen;
    std::ofstream trajstream;
    int nufbeads;


    


public:
    TrjToGro(){}
    TrjToGro(InitialSet initset);
    ~TrjToGro();

    virtual void getSpecificParameters(){}
    virtual void initializeVariables(){}
    virtual void calculateStep(int step){}
    virtual void reduceResults(){}
    virtual void generateOutFileName(){}
    virtual void noramlizeResults(){}
    virtual void dumpGro(int step, real t);
    virtual void writeOutput();
};
};



#endif
