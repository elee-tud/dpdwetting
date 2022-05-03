#ifndef __INTEGRATOR__HPP
#define __INTEGRATOR__HPP

#include "initialization.hpp"

namespace dpd{

class Integrator{
protected:
    Initialization* init;
    Decomposition* decomp;
    SetMPI* mpi;
    Control* control;
    real dt;
    real hdtsqr;
    CellList cells;
    ParticleList particles;
    Ivec mybeads;

public:
    Integrator(){}
    Integrator(Initialization* init);
    ~Integrator(){}

    virtual void updatePosition()=0;
    virtual void updatePosition(bool save_prev)=0;
    virtual void updateVelocity()=0;
    virtual void updateVelocity(bool save_prev)=0;
    virtual real getMaxForce()=0; 
};
};
#endif
