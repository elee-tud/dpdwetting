#ifndef __SLIPSPRINGHARMONIC__HPP
#define __SLIPSPRINGHARMONIC__HPP

#include "../bonded.hpp"
namespace dpd{

class SlipSpringHarmonic:public Bonded{
private:
    real ssk;
    real ssl;
    Control* control;

public:
    SlipSpringHarmonic(){}
    SlipSpringHarmonic(Topology* topol, Configuration* config, Decomposition* decomp);
    ~SlipSpringHarmonic(){}

    virtual Ivec getParticleIndexForCommunication(Particle* ptcl);
    virtual int calculateParticleForce(Particle* ptcl);
    virtual int calculateParticleForce(Particle* ptcl, Real3D com, real** press, real dr);

};
};


#endif
