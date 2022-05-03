#ifndef __POLYMERSMSF__HPP
#define __POLYMERSMSF__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
#define NUMAVG 1
namespace dpd{

class PolymerStructFactor:public Property{
private:
    MoleculeGroup polymers;
    int num_pol;
    int num_at;
    int first_idx;
    ParticleList ref;
    Rvec q;



    


public:
    PolymerStructFactor(){}
    PolymerStructFactor(InitialSet initset);
    ~PolymerStructFactor();

    virtual void getSpecificParameters(){}
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    virtual void normalizeResults();
    virtual void writeOutput();
};
};



#endif
