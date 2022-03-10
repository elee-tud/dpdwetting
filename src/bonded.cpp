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
        Ivec index=getParticleIndexForCommunication(particles[beads[j]]);
        for(int k=0;k<index.size();k++){
            need_here.push_back(index[k]);
        }
    }
    need_num=need_here.size();
//    std::cout << "Particle needed for rank " << mpi->rank() << ":" << need_here << std::endl;
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
        MPI_Send(&need_num, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
        if(need_num>0)
            MPI_Send(&need_here[0], need_num, MPI_INT, MASTER, 1, MPI_COMM_WORLD);
    }
//    if(mpi->isMaster()) std::cout << "need="<<need_num << std::endl;
    /*Broadcast the list of communicating particles*/
    MPI_Bcast(&need_num, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    if(!mpi->isMaster()){
        need_whole=Ivec(need_num);
        need_from=Ivec(need_num);
    }
    
    MPI_Bcast(&need_whole[0], need_num, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&need_from[0], need_num, MPI_INT, MASTER, MPI_COMM_WORLD);
//    std::cout << "Rank " << mpi->rank() <<":" << need_num << "," << need_whole << std::endl;
//    if(mpi->isMaster()) std::cout << "what" << need_whole << std::endl;
//    if(mpi->isMaster()) std::cout << "from" << need_from << std::endl;

    /*Check which domain has the particle*/
    where=Ivec(need_num);
//    std::cout << need_num << " in rank " << mpi->rank() << std::endl;
//    if(need_num>0){
        for(int i=0;i<need_num;i++){

            if(particles[need_whole[i]]->existsHere()==TRUEPTCL)
                exist=true;
            else
                exist=false;

//                std::cout << "particle " << need_whole[i] << " exists in rank " << mpi->rank() << std::endl;

            if(mpi->isMaster()){
//                int existnum=0;
                if(exist){
                    where[i]=MASTER;
//                    existnum++;
                }
                for(int j=1;j<mpi->size();j++){
//                    std::cout << &exist << std::endl;
//                    std::cout << "0" << std::endl;
                    MPI_Recv(&exist, 1, MPI_CXX_BOOL, j, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//                    std::cout << "1" << std::endl;
                    if(exist){
                        where[i]=j;
//                        existnum++;
                    }
                }
//                std::cout << existnum << std::endl;
            }
            else{
                MPI_Send(&exist, 1, MPI_CXX_BOOL, MASTER, 1, MPI_COMM_WORLD);
            }
        }
//    if(mpi->isMaster()) std::cout << "where=" << where << std::endl;
        /*Broadcasting the position of the particles*/
        MPI_Bcast(&where[0], need_num, MPI_INT, MASTER, MPI_COMM_WORLD);
        

        /*Broadcasting the particle positions*/
        /*
        for(int i=0;i<need_num;i++){
            MPI_Bcast(&(particles[need_whole[i]]->coord[0]), 3, MPI_DOUBLE, where[i], MPI_COMM_WORLD);
        }
        */
        for(int i=0;i<need_num;i++){
            if(mpi->rank()==where[i])
                MPI_Send(&(particles[need_whole[i]]->coord[0]), 3, MPI_DOUBLE, need_from[i], need_whole[i], MPI_COMM_WORLD);
            if(mpi->rank()==need_from[i])
                MPI_Recv(&(particles[need_whole[i]]->coord[0]), 3, MPI_DOUBLE, where[i], need_whole[i], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

//    }
    return;

}





void Bonded::reduceForce(){
    need_num=need_whole.size();
//    if(mpi->isMaster()){
//        std::cout << mpi->rank() <<  need_whole << need_from << where << std::endl;
//    }
    for(int i=0;i<need_num;i++){
        if(mpi->rank()==need_from[i]){
            MPI_Send(&(particles[need_whole[i]]->force[0]), 3, MPI_DOUBLE, where[i], need_whole[i], MPI_COMM_WORLD);
//            if(need_whole[i]==1274)
//                std::cout << "sent_force=" << particles[1274]->force << "in rank " << mpi->rank() << std::endl;
            MPI_Send(&(particles[need_whole[i]]->stress[0]), 9, MPI_DOUBLE, where[i], need_whole[i]+10000000, MPI_COMM_WORLD);
        }
        if(mpi->rank()==where[i]){
            Real3D force_recv(0.0);
            Rvec stress_recv(9,0.);
            MPI_Recv(&force_recv[0], 3, MPI_DOUBLE, need_from[i], need_whole[i], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&stress_recv[0], 9, MPI_DOUBLE, need_from[i], need_whole[i]+10000000, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            particles[need_whole[i]]->force+=force_recv;
//            if(need_whole[i]==1274)
//                std::cout << "added_force=" << force_recv << "in rank " << mpi->rank() << std::endl;;
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
//    if(mpi->rank()==0)
//        std::cout << "1273 force before: " << particles[1273]->force << std::endl;
//    if(mpi->rank()==4)
//        std::cout << "1274 force before: " << particles[1274]->force << std::endl;
    communicateBondedParticles();
    for(int j=0;j<beads.size();j++){
        calc_nbonds+=calculateParticleForce(particles[beads[j]]);
    }
    reduceForce();
//    if(mpi->rank()==0)
//        std::cout << "1273 force after: " << particles[1273]->force << std::endl;
//    if(mpi->rank()==4)
//        std::cout << "1274 force after: " << particles[1274]->force << std::endl;
    int calc_nbonds_tot=0;
//    MPI_Reduce(&calc_nbonds, &calc_nbonds_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
//    if(mpi->isMaster()) std::cout << calc_nbonds_tot << "," ;
//    MPI_Barrier(MPI_COMM_WORLD);
//    err.missingBonds(calc_nbonds_tot, topol->getNbonds());
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
    /*
    int calc_nbonds_tot=0;
    MPI_Reduce(&calc_nbonds, &calc_nbonds_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    err.missingBonds(calc_nbonds_tot, topol->getNbonds());
    */
    return;
}
