#ifndef __BONDLENHARMONIC__HPP
#define __BONDLENHARMONIC__HPP

#include "../bonded.hpp"

namespace dpd{

class BondLenHarmonic:public Bonded{
protected:
    Rvec k;
    Rvec leq;
public:
    BondLenHarmonic(){}
    BondLenHarmonic(Topology* topol, Configuration* config, Decomposition* decomp);
    ~BondLenHarmonic(){}

    Ivec getParticleIndexForCommunication(Particle* ptcl);
    virtual int calculateParticleForce(Particle* ptcl);
    virtual int calculateParticleForce(Particle* ptcl, Real3D com, real** press, real dr);
};
};

#endif
