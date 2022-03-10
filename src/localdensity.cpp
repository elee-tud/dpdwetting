#include "localdensity.hpp"

using namespace dpd;
LocalDensity::LocalDensity(Initialization* init):init(init){
    decomp=init->getDecomposition();
    mpi=decomp->getMPI();
    particles=decomp->getParticles();
    cells=decomp->getCells();
    truecells=decomp->getTrueCells();
    ghostcells=decomp->getGhostCells();
    config=init->getConfiguration();
    box=config->getBox();
    pbc=PeriodicBoundary(box);

    topol=init->getTopology();
    Brcut=topol->getBrcut();
    normconst=7.5/PI;
}

LocalDensity::LocalDensity(Decomposition* decomp):decomp(decomp){
    mpi=decomp->getMPI();
    particles=decomp->getParticles();
    cells=decomp->getCells();
    truecells=decomp->getTrueCells();
    ghostcells=decomp->getGhostCells();
    config=decomp->getConfiguration();
    box=config->getBox();
    pbc=PeriodicBoundary(box);

    topol=decomp->getTopology();
    Brcut=topol->getBrcut();
    normconst=7.5/PI;
}

void LocalDensity::setZeroLocalDensity(){
    cells=decomp->getCells();
    truecells=decomp->getTrueCells();
    ghostcells=decomp->getGhostCells();
    box=config->getBox();
    pbc=PeriodicBoundary(box);
    for(int i=0;i<cells.size();i++){
        Ivec beads=cells[i]->getBeads();
        for(int j=0;j<beads.size();j++){
            particles[beads[j]]->density=0.0;
        }
    }
    return;
}

void LocalDensity::calculateLocalDensity(){
    setZeroLocalDensity();
    for(int i=0;i<truecells.size();i++){
        calculateDensityInTheCells(cells[truecells[i]]);
        Ivec nbidx=cells[truecells[i]]->getNeighborCells();
        for(int j=0;j<nbidx.size();j++){
            calculateDensityBetweenCells(cells[truecells[i]], cells[nbidx[j]]);
        }
    }
    reduceLocalDensity();
    redistLocalDensity();
    return;
}
void LocalDensity::calculateDensityInTheCells(Cell* cell){
    Ivec mybeads=cell->getBeads();
    int num_mybeads=mybeads.size();
    for(int i=0;i<num_mybeads;i++){
        for(int j=i+1;j<num_mybeads;j++){
            calculatePairDensity(mybeads[i], mybeads[j]);
        }
    }
    return;
}

void LocalDensity::calculateDensityBetweenCells(Cell* cell1, Cell* cell2){
    Ivec mybeads1=cell1->getBeads();
    Ivec mybeads2=cell2->getBeads();
    int num_mybeads1=mybeads1.size();
    int num_mybeads2=mybeads2.size();
    for(int i=0;i<num_mybeads1;i++){
        for(int j=0;j<num_mybeads2;j++){
            calculatePairDensity(mybeads1[i], mybeads2[j]);
        }
    }
    return;
}

void LocalDensity::calculatePairDensity(int index1, int index2){
    Real3D rij=pbc.getMinimumImageVector(particles[index2]->coord, particles[index1]->coord);
    real dist=rij.abs();
    int bt1=particles[index1]->getParticleType();
    int bt2=particles[index2]->getParticleType();
    if(dist<Brcut[bt1][bt2]){
        real rho=1.-dist/Brcut[bt1][bt2];
        rho=rho*rho/pow(Brcut[bt1][bt2],3.0)*normconst;
        particles[index2]->density+=rho;
        particles[index1]->density+=rho;
    }



    return;
}
void LocalDensity::reduceLocalDensity(){
    Ivec2D ptcls_to_send=decomp->getProcsOfGhosts();
    Ivec2D ptcls_to_recv=decomp->getProcsOfRealBeads();


    for(int i=0;i<mpi->size();i++){
        if(i!=mpi->rank()){
            int sendsize=ptcls_to_send[i].size();
            int recvsize=ptcls_to_recv[i].size();
            Rvec senddensity(sendsize);
            Rvec recvdensity(recvsize);
            for(int j=0;j<sendsize;j++){
                senddensity[j]=particles[ptcls_to_send[i][j]]->density;
            }
            MPI_Sendrecv(&senddensity[0], sendsize, MPI_DOUBLE, i, 1, &recvdensity[0], recvsize, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for(int j=0;j<recvsize;j++){
                particles[ptcls_to_recv[i][j]]->density+=recvdensity[j];
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    return;
}

//Redistribution of local density of ghost particles//
void LocalDensity::redistLocalDensity(){
    Ivec2D ptcls_to_send=decomp->getProcsOfRealBeads();
    Ivec2D ptcls_to_recv=decomp->getProcsOfGhosts();
    for(int i=0;i<mpi->size();i++){
        if(mpi->rank()!=i){
            int sendsize=ptcls_to_send[i].size();
            int recvsize=ptcls_to_recv[i].size();;
            Rvec senddensity(sendsize);
            Rvec recvdensity(recvsize);
            for(int j=0;j<sendsize;j++)
                senddensity[j]=particles[ptcls_to_send[i][j]]->density;
            MPI_Sendrecv(&senddensity[0], sendsize, MPI_DOUBLE, i, 1, &recvdensity[0], recvsize, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for(int j=0;j<recvsize;j++){
                particles[ptcls_to_recv[i][j]]->density=recvdensity[j];
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    return;
}



