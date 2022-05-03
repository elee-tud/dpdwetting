#ifndef __SPRINGPULLING__HPP
#define __SPRINGPULLING__HPP

#include "../externalforce.hpp"
namespace dpd{

class SpringPulling:public ExternalForce{
private:
    Real3D center;
    real springk;
    int direct;
    Ivec dir;
    int numdir;


public:
    SpringPulling(){}
    SpringPulling(Topology* topol, Configuration* config, Decomposition* decomp);
    ~SpringPulling(){}

    void calculateSingleForce(Particle* ptcl);
};
};
#endif

