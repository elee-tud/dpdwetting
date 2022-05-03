#ifndef __GRAVITY__HPP
#define __GRAVITY__HPP

#include "../externalforce.hpp"
namespace dpd{

class Gravity:public ExternalForce{
private:
    real field;
    int direct;

public:
    Gravity(){}
    Gravity(Topology* topol, Configuration* config, Decomposition* decomp);
    ~Gravity(){}

    void calculateSingleForce(Particle* ptcl);
};
};
#endif

