#ifndef __MPICLASSES__HPP
#define __MPICLASSES__HPP

#include "setmpi.hpp"
#include "molecule.hpp"
#include "particle.hpp"

namespace dpd{

class MPIClasses{
private:
    SetMPI* mpi;
    bool isMaster;
    int rank;
    int size;
    Ivec mybeads;

public:
    MPIClasses(){}
    MPIClasses(SetMPI* mpi);
    ~MPIClasses(){}

    void bcastMolecule(int source, Molecule *mol);
    void bcastParticleCoord(Particle* ptcl, int source);
    void bcastParticleVeloc(Particle* ptcl, int source);
    void bcastParticleForce(Particle* ptcl, int source);
    void sendParticleCoord(Particle* ptcl, int dest, int tag);
    void recvParticleCoord(Particle* ptcl, int source, int tag);
    void sendParticleVeloc(Particle* ptcl, int dest, int tag);
    void recvParticleVeloc(Particle* ptcl, int source, int tag);
    void sendParticleForce(Particle* ptcl, int dest, int tag);
    void recvParticleForce(Particle* ptcl, int source, int tag);
    Real3D recvParticleGhostForce( Particle *ptcl, int source, int tag);
};
};
#endif
