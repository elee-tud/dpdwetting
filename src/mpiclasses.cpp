#include "mpiclasses.hpp"
using namespace dpd;
MPIClasses::MPIClasses(SetMPI* mpi):mpi(mpi){
    rank=mpi->rank();
    size=mpi->size();
    if(rank==MASTER)
        isMaster=true;
    else
        isMaster=false;
}

void MPIClasses::bcastMolecule(int source, Molecule *mol){
    int namesize=mol->name.size();
    int atomsize=mol->atoms.size();
    int bondsize=mol->bonds.size();
    int atypesize=mol->atomtypes.size();
    int btypesize=mol->bondtypes.size();
    MPI_Bcast(&namesize, 1, MPI_INT, source, MPI_COMM_WORLD);
    MPI_Bcast(&atomsize, 1, MPI_INT, source, MPI_COMM_WORLD);
    MPI_Bcast(&bondsize, 1, MPI_INT, source, MPI_COMM_WORLD);
    MPI_Bcast(&atypesize, 1, MPI_INT, source, MPI_COMM_WORLD);
    MPI_Bcast(&btypesize, 1, MPI_INT, source, MPI_COMM_WORLD);
    MPI_Bcast(&mol->is_frozen, 1, MPI_CXX_BOOL, source, MPI_COMM_WORLD);


    if(!mpi->isMaster()){
        mol->name.resize(namesize);
        mol->atoms.resize(atomsize);
        mol->masses.resize(atomsize);
        mol->bonds.resize(bondsize);
        mol->atomtypes.resize(atypesize);
        mol->bondtypes.resize(btypesize);
        for(int i=0;i<bondsize;i++)
            mol->bonds[i].resize(2);
    }

    MPI_Bcast(&mol->num_mols, 1, MPI_INT, source, MPI_COMM_WORLD);
    MPI_Bcast(&mol->num_atoms, 1, MPI_INT, source, MPI_COMM_WORLD);
    MPI_Bcast(&mol->num_bonds, 1, MPI_INT, source, MPI_COMM_WORLD);

    MPI_Bcast(const_cast<char*>(mol->name.data()), namesize, MPI_CHAR, source, MPI_COMM_WORLD);
    MPI_Bcast(&mol->atoms[0], atomsize, MPI_INT, source, MPI_COMM_WORLD);
    MPI_Bcast(&mol->masses[0], atomsize, MPI_DOUBLE, source, MPI_COMM_WORLD);
    for(int i=0;i<bondsize;i++)
        MPI_Bcast(&mol->bonds[i][0], 2, MPI_INT, source, MPI_COMM_WORLD);
    MPI_Bcast(&mol->atomtypes[0], atypesize, MPI_INT, source, MPI_COMM_WORLD);
    MPI_Bcast(&mol->bondtypes[0], btypesize, MPI_INT, source, MPI_COMM_WORLD);
    return;
}


void MPIClasses::bcastParticleCoord(Particle* ptcl, int source){
    real c[3]={ptcl->coord[0], ptcl->coord[1], ptcl->coord[2]};
    MPI_Bcast(&c[0], 3, MPI_DOUBLE, source, MPI_COMM_WORLD);
    ptcl->coord=Real3D{c[0], c[1], c[2]};
    return;
}

void MPIClasses::bcastParticleVeloc(Particle* ptcl, int source){
    real c[3]={ptcl->veloc[0], ptcl->veloc[1], ptcl->veloc[2]};
    MPI_Bcast(&c[0], 3, MPI_DOUBLE, source, MPI_COMM_WORLD);
    ptcl->veloc=Real3D{c[0], c[1], c[2]};
    return;
}

void MPIClasses::bcastParticleForce(Particle* ptcl, int source){
    real c[3]={ptcl->force[0], ptcl->force[1], ptcl->force[2]};
    MPI_Bcast(&c[0], 3, MPI_DOUBLE, source, MPI_COMM_WORLD);
    ptcl->force=Real3D{c[0], c[1], c[2]};
    return;
}

void MPIClasses::sendParticleCoord(Particle *ptcl, int dest, int tag){
    real c[3]={ptcl->coord[0], ptcl->coord[1], ptcl->coord[2]};
    MPI_Send(&ptcl->coord[0], 3, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
    return;
}

void MPIClasses::recvParticleCoord( Particle *ptcl, int source, int tag){
    real c[3];
    MPI_Recv(&c[0], 3, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    ptcl->coord=Real3D{c[0], c[1], c[2]};
    return;
}

void MPIClasses::sendParticleVeloc(Particle *ptcl, int dest, int tag){
    real c[3]={ptcl->veloc[0], ptcl->veloc[1], ptcl->veloc[2]};
    MPI_Send(&ptcl->veloc[0], 3, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
    return;
}

void MPIClasses::recvParticleVeloc( Particle *ptcl, int source, int tag){
    real c[3];
    MPI_Recv(&c[0], 3, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    ptcl->veloc=Real3D{c[0], c[1], c[2]};
    return;
}

void MPIClasses::sendParticleForce(Particle *ptcl, int dest, int tag){
    real c[3]={ptcl->force[0], ptcl->force[1], ptcl->force[2]};
    MPI_Send(&ptcl->force[0], 3, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
    return;
}

void MPIClasses::recvParticleForce( Particle *ptcl, int source, int tag){
    real c[3];
    MPI_Recv(&c[0], 3, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    ptcl->force=Real3D{c[0], c[1], c[2]};
    return;
}
Real3D MPIClasses::recvParticleGhostForce( Particle *ptcl, int source, int tag){
    real c[3];
    MPI_Recv(&c[0], 3, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    return Real3D{c[0], c[1], c[2]};
}
