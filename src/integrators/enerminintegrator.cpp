#include "enerminintegrator.hpp"
#include <algorithm>

using namespace dpd;
EnerMinIntegrator::EnerMinIntegrator(Initialization* init):Integrator(init){
    maxdisp=control->getCellCutoff()*0.05;

}

void EnerMinIntegrator::updatePosition(bool save_prev){
    updatePosition();
}

void EnerMinIntegrator::updatePosition(){
    findMaxForce();
    real disp=maxdisp/maxforce;
    mybeads=decomp->getBeadsIndexInDomain();
    real distmax=0.;
    for(int i=0;i<mybeads.size();i++){
        if(!particles[mybeads[i]]->isFrozen()){
            real dist=disp*particles[mybeads[i]]->force.abs();
            if(dist>distmax)
                distmax=dist;
            particles[mybeads[i]]->coord+=disp*particles[mybeads[i]]->force;
        }
    }
    /*
    Rvec recvdist(mpi->size());
    MPI_Gather(&distmax, 1, MPI_DOUBLE, &recvdist[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(mpi->isMaster()){
        distmax=*std::max_element(recvdist.begin(), recvdist.end());
        std::cout << distmax << std::endl;
    }
    */

    return;
}

void EnerMinIntegrator::findMaxForce(){
    mybeads=decomp->getBeadsIndexInDomain();
    maxforce=0.0;
    for(int i=0;i<mybeads.size();i++){
        real force=particles[mybeads[i]]->force.abs();
        if(force>maxforce)
            maxforce=force;
    }
//    std::cout << mpi->rank() << ": " << maxforce << std::endl;
    Rvec recvforce(mpi->size());
    MPI_Gather(&maxforce, 1, MPI_DOUBLE, &recvforce[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /*
    if(mpi->isMaster()){
        Rvec recvforce(mpi->size());
        recvforce[0]=maxforce;
        for(int i=1;i<mpi->rank();i++){
            MPI_Recv(&recvforce[i], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    else{
        MPI_Send(&maxforce, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    */
    if(mpi->isMaster()){
        maxforce=*std::max_element(recvforce.begin(), recvforce.end());
    }
    MPI_Bcast(&maxforce, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    return;
}






