#ifndef __BONDED__HPP
#define __BONDED__HPP

#include "decomposition.hpp"
#include "stress.hpp"

namespace dpd{

class Bonded{
protected:
    Topology* topol;
    Configuration* config;
    Decomposition* decomp;
    Real3D box;
    PeriodicBoundary pbc;
    ParticleList particles;
    CellList cells;
    Ivec truecells;
    Ivec ghostcells;
    Control* control;
    SetMPI* mpi;
    Error err;
    Ivec need_whole;
    Ivec need_from;
    Ivec where;
    int need_num;
    bool exist;

public:
    Bonded(){}
    Bonded(Topology* topol, Configuration* config, Decomposition* decomp);
    ~Bonded(){}

    void calculateForce();
    void calculateForce(Real3D com, real** press, real dr);
    void communicateBondedParticles();
    void reduceForce();
    virtual Ivec getParticleIndexForCommunication(Particle* ptcl)=0;
    virtual int calculateParticleForce(Particle* ptcl)=0;
    virtual int calculateParticleForce(Particle* ptcl, Real3D com, real** press, real dr)=0;

};

typedef std::vector<Bonded*> BondedList;
};
#endif

