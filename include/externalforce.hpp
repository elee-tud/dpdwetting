#ifndef __EXTERNALFORCE__HPP
#define __EXTERNALFORCE__HPP
#include "decomposition.hpp"

namespace dpd{

class ExternalForce{
protected:
    Topology* topol;
    Configuration* config;
    Decomposition* decomp;
    Real3D box;
    PeriodicBoundary pbc;
    ParticleList particles;
    CellList cells;
    Control *control;
    SetMPI* mpi;
    Ivec mybeads;

public:
    ExternalForce(){}
    ExternalForce(Topology* topol, Configuration* config, Decomposition* decomp);
    ~ExternalForce(){}

    void calculateForce();
    virtual void calculateSingleForce(Particle* ptcl)=0;
};

typedef std::vector<ExternalForce*> ExternalList;
};

#endif
