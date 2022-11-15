#include "bonded.hpp"

using namespace dpd;
Bonded::Bonded(Topology* topol, Configuration* config, Decomposition* decomp):topol(topol), config(config), decomp(decomp){
    cells=decomp->getCells();
    truecells=decomp->getTrueCells();
    ghostcells=decomp->getGhostCells();
    particles=config->getParticles();
    box=config->getBox();
    pbc=PeriodicBoundary(box);
    control=decomp->getControl();
    mpi=decomp->getMPI();
    err=Error(mpi);
    need_whole.reserve(particles.size());
    need_from.reserve(particles.size());
}

/*Communicatation of the position of the bonded particles which are not in the same domain*/
void Bonded::communicateBondedParticles(){
    need_whole.clear();
    need_from.clear();
    Ivec need_here;
    need_here.reserve(particles.size());        
    Ivec beads=decomp->getBeadsIndexInDomain();
    /*Finding bonded particles not in the same domain*/
    for(int j=0;j<beads.size();j++){
        Ivec index=getParticleIndexForCommunication(particles[beads[j]]);       //Which particle information is needed in this domain
        for(int k=0;k<index.size();k++){
            need_here.push_back(index[k]);
        }
    }
    need_num=need_here.size();
    /*Building a list of particles to be comunicated in the master*/
    if(mpi->isMaster()){
        Ivec need_recv;
        for(int i=0;i<need_num;i++){
            need_whole.push_back(need_here[i]);
            need_from.push_back(MASTER);
        }
        for(int i=1;i<mpi->size();i++){
            MPI_Recv(&need_num, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            need_recv.reserve(need_num);
            if(need_num>0)
                MPI_Recv(&need_recv[0], need_num, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for(int j=0;j<need_num;j++){
                need_whole.push_back(need_recv[j]);
                need_from.push_back(i);
            }
        }
        need_num=need_whole.size();
    }
    else{
        MPI_Send(&need_num, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);     //Sending the number of particles to be needed
        if(need_num>0)
            MPI_Send(&need_here[0], need_num, MPI_INT, MASTER, 1, MPI_COMM_WORLD);  //Sending the particle index to be needed
    }
    /*Broadcast the list of communicating particles*/
    MPI_Bcast(&need_num, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    if(!mpi->isMaster()){
        need_whole=Ivec(need_num);
        need_from=Ivec(need_num);
    }
    
    MPI_Bcast(&need_whole[0], need_num, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&need_from[0], need_num, MPI_INT, MASTER, MPI_COMM_WORLD);

    /*Check which domain has the particle*/
    where=Ivec(need_num);
    for(int i=0;i<need_num;i++){

        if(particles[need_whole[i]]->existsHere()==TRUEPTCL)
            exist=true;
        else
            exist=false;

        if(mpi->isMaster()){
            if(exist){
                where[i]=MASTER;
            }
            for(int j=1;j<mpi->size();j++){
                MPI_Recv(&exist, 1, MPI_CXX_BOOL, j, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if(exist){
                    where[i]=j;
                }
            }
        }
        else{
            MPI_Send(&exist, 1, MPI_CXX_BOOL, MASTER, 1, MPI_COMM_WORLD);
        }
    }
    /*Broadcasting the position of the particles*/
    MPI_Bcast(&where[0], need_num, MPI_INT, MASTER, MPI_COMM_WORLD);
    

    /*Broadcasting the particle positions*/
    for(int i=0;i<need_num;i++){
        if(mpi->rank()==where[i])
            MPI_Send(&(particles[need_whole[i]]->coord[0]), 3, MPI_DOUBLE, need_from[i], need_whole[i], MPI_COMM_WORLD);
        if(mpi->rank()==need_from[i])
            MPI_Recv(&(particles[need_whole[i]]->coord[0]), 3, MPI_DOUBLE, where[i], need_whole[i], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    return;

}





void Bonded::reduceForce(){
    need_num=need_whole.size();
    for(int i=0;i<need_num;i++){
        if(mpi->rank()==need_from[i]){
            /*Sending forces*/
            MPI_Send(&(particles[need_whole[i]]->force[0]), 3, MPI_DOUBLE, where[i], need_whole[i], MPI_COMM_WORLD);
            MPI_Send(&(particles[need_whole[i]]->stress[0]), 9, MPI_DOUBLE, where[i], need_whole[i], MPI_COMM_WORLD);
        }
        if(mpi->rank()==where[i]){
            Real3D force_recv(0.0);
            Rvec stress_recv(9,0.);
            /*Receiving forces*/
            MPI_Recv(&force_recv[0], 3, MPI_DOUBLE, need_from[i], need_whole[i], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&stress_recv[0], 9, MPI_DOUBLE, need_from[i], need_whole[i], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            particles[need_whole[i]]->force+=force_recv;
            /*Adding stress tensor contributed by bonded potential*/
            for(int j=0;j<9;j++){
                particles[need_whole[i]]->stress[j]+=stress_recv[j];        
            }
        }
    }
    return;
}




void Bonded::calculateForce(){
    cells=decomp->getCells();
    truecells=decomp->getTrueCells();
    ghostcells=decomp->getGhostCells();
    pbc=PeriodicBoundary(decomp->getBox());
    int calc_nbonds=0;
    Ivec beads=decomp->getBeadsIndexInDomain();
    communicateBondedParticles();       //Communicating positions of bonded particles
    for(int j=0;j<beads.size();j++){
        calc_nbonds+=calculateParticleForce(particles[beads[j]]);   //Calculating bonded force in the current domain
    }
    reduceForce();
    return;
}

void Bonded::calculateForce(Real3D com, real** press, real dr){
    int calc_nbonds=0;
    cells=decomp->getCells();
    truecells=decomp->getTrueCells();
    ghostcells=decomp->getGhostCells();
    pbc=PeriodicBoundary(decomp->getBox());
    Ivec beads=decomp->getBeadsIndexInDomain();
    for(int j=0;j<beads.size();j++){
        calc_nbonds+=calculateParticleForce(particles[beads[j]], com, press, dr);
    }
    return;
}
