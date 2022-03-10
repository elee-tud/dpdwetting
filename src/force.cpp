#include "force.hpp"

using namespace dpd;
Force::Force(Initialization* init):init(init){
    decomp=init->getDecomposition();
    mpi=init->getMPI();
    inter=init->getInteractions();
    nonbondedlist=inter->getNonbondedInteractions();
    externallist=inter->getExternalInteractions();
    bondedlist=inter->getBondedInteractions();
    mpiptcl=MPIClasses(mpi);
    particles=decomp->getParticles();
    cells=decomp->getCells();


}

Force::Force(Decomposition* decomp, Interactions* inter):decomp(decomp), inter(inter){
    mpi=decomp->getMPI();
    nonbondedlist=inter->getNonbondedInteractions();
    externallist=inter->getExternalInteractions();
    bondedlist=inter->getBondedInteractions();
    mpiptcl=MPIClasses(mpi);
    particles=decomp->getParticles();
    cells=decomp->getCells();


}


void Force::setZeroForce(){
    cells=decomp->getCells();
    for(int i=0;i<cells.size();i++){
        Ivec beads=cells[i]->getBeads();
        for(int j=0;j<beads.size();j++){
            particles[beads[j]]->force=Real3D(0.0, 0.0, 0.0);
            particles[beads[j]]->stress=Rvec(9,0.0);
        }
    }
    return;
}


void Force::calculateForce(){
    setZeroForce();
    for(int i=0;i<nonbondedlist.size();i++){
        nonbondedlist[i]->calculateForce();
    }
    for(int i=0;i<bondedlist.size();i++){
        bondedlist[i]->calculateForce();
    }
    for(int i=0;i<externallist.size();i++){
        externallist[i]->calculateForce();
    }

    reduceForce();
    return;
}

void Force::calculateForce(Real3D com, real** press, real dr){
    setZeroForce();
    for(int i=0;i<nonbondedlist.size();i++){
        nonbondedlist[i]->calculateForce(com, press, dr);
    }
    for(int i=0;i<bondedlist.size();i++){
        bondedlist[i]->calculateForce();
    }
    for(int i=0;i<externallist.size();i++){
        externallist[i]->calculateForce();
    }

    return;
}

void Force::reduceForce(){
    Ivec2D ptcls_to_send=decomp->getProcsOfGhosts();
    Ivec2D ptcls_to_recv=decomp->getProcsOfRealBeads();
    for(int i=0;i<mpi->size();i++){
        if(mpi->rank()!=i){
            int sendsize=ptcls_to_send[i].size();
            int recvsize=ptcls_to_recv[i].size();
            Rvec sendforce=serializeNForces(ptcls_to_send[i], particles);
            Rvec sendstress=serializeNStresses(ptcls_to_send[i], particles);
            Rvec ghostforce(3*recvsize);
            MPI_Sendrecv(&sendforce[0], 3*sendsize, MPI_DOUBLE, i, 1, &ghostforce[0], 3*recvsize, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for(int j=0;j<recvsize;j++){
                particles[ptcls_to_recv[i][j]]->force[0]+=ghostforce[j*3];
                particles[ptcls_to_recv[i][j]]->force[1]+=ghostforce[j*3+1];
                particles[ptcls_to_recv[i][j]]->force[2]+=ghostforce[j*3+2];
            }
            Rvec ghoststress(9*recvsize);
            MPI_Sendrecv(&sendstress[0], 9*sendsize, MPI_DOUBLE, i, 2, &ghoststress[0], 9*recvsize, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for(int j=0;j<recvsize;j++){
                particles[ptcls_to_recv[i][j]]->stress[0]+=ghoststress[j*9];
                particles[ptcls_to_recv[i][j]]->stress[1]+=ghoststress[j*9+1];
                particles[ptcls_to_recv[i][j]]->stress[2]+=ghoststress[j*9+2];
                particles[ptcls_to_recv[i][j]]->stress[3]+=ghoststress[j*9+3];
                particles[ptcls_to_recv[i][j]]->stress[4]+=ghoststress[j*9+4];
                particles[ptcls_to_recv[i][j]]->stress[5]+=ghoststress[j*9+5];
                particles[ptcls_to_recv[i][j]]->stress[6]+=ghoststress[j*9+6];
                particles[ptcls_to_recv[i][j]]->stress[7]+=ghoststress[j*9+7];
                particles[ptcls_to_recv[i][j]]->stress[8]+=ghoststress[j*9+8];
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    return;
}

