#include "decomposition.hpp" 
#include <algorithm>
#include <unistd.h>

using namespace dpd;
Decomposition::Decomposition(Control* control, Configuration* config, SetMPI* mpi):control(control), config(config), mpi(mpi){

    topology=config->getTopology();
    particles=config->getParticles();
    box=config->getBox();
    numprocs=mpi->size();
    rcut=control->getCellCutoff();
    mpiptcl=MPIClasses(mpi);
    saferatio=control->getSafeRatio();
    maxnbeads=static_cast<int>(topology->getNbeads()/numprocs*saferatio)+1;
    if(maxnbeads>topology->getNbeads())
        maxnbeads=topology->getNbeads();
    dumpfrozen=control->getDumpFrozen();



    pbc=PeriodicBoundary(box);
    err=Error(mpi);

    calculateDomainDivisor();
    indexing_domain=Indexing(num_domains);


    myid=mpi->rank();
    my3did=indexing_domain.get3DIndexFromIndex(myid);
    totnum_cells=0;

    nbsearch=indexing_domain.buildNeighborSearchVector();
    calculateDomainLength();
    printDomainLengthInformation();
    calculateCellLength();
    printCellNumberInformation();
    indexing_cell=Indexing(num_cells);
    makeCellLists();

    mybeads.reserve(maxnbeads);
    ghosts=Ivec2D(numprocs, Ivec{});
    realbeads=Ivec2D(numprocs, Ivec{});
    for(int i=0;i<numprocs;i++){
        ghosts[i].reserve(maxnbeads);
    }
    if(mpi->isMaster()) std::cout << std::endl;
    



}

Decomposition::Decomposition(Control* control, Configuration* config, SetMPI* mpi, int numprocs, real rcut):control(control), config(config), mpi(mpi), numprocs(numprocs), rcut(rcut){

    topology=config->getTopology();
    particles=config->getParticles();
    box=config->getBox();
    mpiptcl=MPIClasses(mpi);
    saferatio=control->getSafeRatio();
    maxnbeads=topology->getNbeads();
    if(maxnbeads>topology->getNbeads())
        maxnbeads=topology->getNbeads();
    dumpfrozen=control->getDumpFrozen();



    pbc=PeriodicBoundary(box);
    err=Error(mpi);

    myid=0;
    calculateDomainDivisor();
    indexing_domain=Indexing(num_domains);

    my3did=indexing_domain.get3DIndexFromIndex(myid);

    totnum_cells=0;

    nbsearch=indexing_domain.buildNeighborSearchVector();
    calculateDomainLength();
    calculateCellLength();
    indexing_cell=Indexing(num_cells);
    makeCellLists();

    mybeads.reserve(maxnbeads);
    ghosts=Ivec2D(numprocs, Ivec{});
    realbeads=Ivec2D(numprocs, Ivec{});
    for(int i=0;i<numprocs;i++){
        ghosts[i].reserve(maxnbeads);
    }
    



}

Decomposition::Decomposition(Control* control, Configuration* config, SetMPI* mpi, bool decomp_anal):control(control), config(config), mpi(mpi){

    topology=config->getTopology();
    particles=config->getParticles();
    box=config->getBox();
    if(decomp_anal==true)
        numprocs=1;
    rcut=control->getCellCutoff();
    mpiptcl=MPIClasses(mpi);
    saferatio=control->getSafeRatio();
    maxnbeads=static_cast<int>(topology->getNbeads()/numprocs*saferatio)+1;
    dumpfrozen=control->getDumpFrozen();



    pbc=PeriodicBoundary(box);
    err=Error(mpi);

    calculateDomainDivisor();
    indexing_domain=Indexing(num_domains);


    myid=mpi->rank();
    my3did=indexing_domain.get3DIndexFromIndex(myid);

    nbsearch=indexing_domain.buildNeighborSearchVector();
    calculateDomainLength();
    calculateCellLength();
    indexing_cell=Indexing(num_cells);
    makeCellLists();

    mybeads.reserve(maxnbeads);
    ghosts=Ivec2D(numprocs, Ivec{});
    realbeads=Ivec2D(numprocs, Ivec{});
    for(int i=0;i<numprocs;i++){
        ghosts[i].reserve(maxnbeads);
    }
    if(mpi->isMaster()) std::cout << std::endl;
    



}

void Decomposition::calculateDomainDivisor(){
    /*Assigning divisor of domains in three directions*/
    if(numprocs==1)  num_domains={1,1,1};
    else if(numprocs==2)  num_domains={2,1,1};
    else if(numprocs==4)  num_domains={2,2,1};
    else if(numprocs==6)  num_domains={3,2,1};
    else if(numprocs==8)  num_domains={2,2,2};
    else if(numprocs==10)  num_domains={5,2,1};
    else if(numprocs==12)  num_domains={3,2,2};
    else if(numprocs==16)  num_domains={4,2,2};
    else if(numprocs==18)  num_domains={3,3,2};
    else if(numprocs==24)  num_domains={4,3,2};
    else if(numprocs==32)  num_domains={4,4,2};
    else if(numprocs==48)  num_domains={4,4,3};
    else if(numprocs==60)  num_domains={5,4,3};
    else if(numprocs==64)  num_domains={4,4,4};
    else if(numprocs==72)  num_domains={6,4,3};
    else if(numprocs==84)  num_domains={7,4,3};
    else if(numprocs==96)  num_domains={6,4,4};
    else err.wrongNumberOfCores();
    return;
}

void Decomposition::calculateDomainLength(){

    domain_length=Real3D(0.0);
    domain_min=Real3D(0.0);
    domain_max=Real3D(0.0);
    
    /*Calculating the length of each direction of the domain*/
    domain_length[0]=box[0]/num_domains[0];
    domain_length[1]=box[1]/num_domains[1];
    domain_length[2]=box[2]/num_domains[2];

    /*Calculating the lower boundary of the domain*/
    domain_min[0]=my3did[0]*domain_length[0];
    domain_min[1]=my3did[1]*domain_length[1];
    domain_min[2]=my3did[2]*domain_length[2];
    
    /*Calculating the upper boundary of the domain*/
    domain_max[0]=(my3did[0]+1)*domain_length[0];
    domain_max[1]=(my3did[1]+1)*domain_length[1];
    domain_max[2]=(my3did[2]+1)*domain_length[2];
    return;
}




void Decomposition::printDomainLengthInformation(){
    if(mpi->isMaster()){
        std::cout << " Domain decomposition into " << num_domains[0] << "*" << num_domains[1] << "*" << num_domains[2] << " domains" << std::endl;
        std::cout << " Domain size = " << domain_length[0] << "*" << domain_length[1] << "*" << domain_length[2] << std::endl;
    }
    return;
}
void Decomposition::clearAllBeadInformation(){
    mybeads.clear();        /*Removing all beads from the domain*/
    for(int i=0;i<cells.size();i++){
        cells[i]->clear();          /*Removing all beads from all cells in the domain*/
    }
    return;
}

void Decomposition::allocateBeadsToDomain(){

    int num_beads=topology->getNbeads();
    mybeads.clear();
    for(int i=0;i<num_beads;i++){
        int idx=getNewDomainIndex(particles[i]);        /*Calculating the domain index of a particle*/
        if(idx==myid){          /*If its in my domain*/
            mybeads.push_back(i);           /*It's my particle*/
            particles[i]->setTrue();
        }
        else
            particles[i]->setNone();        /*It's not my particle*/
    } 
    num_mybeads=mybeads.size();
    
    return;
}


void Decomposition::allocateBeadsToCells(){
    for(int i=0;i<num_mybeads;i++){
        int idx=getNewDomainCellIndex(particles[mybeads[i]])[1];        /*Calculating the cell index if my particles*/
        cells[idx]->addBeads(mybeads[i]);           /*Adding the particle to a proper cell*/
    }

    return;
}



int Decomposition::getNeighborProcessorID(Ivec dist){           
    Ivec neighbor=indexing_domain.addIndexToIndex(my3did, dist);        /*Obtaining the index of neighbor domains*/
    return indexing_domain.getIndexFrom3DIndex(neighbor);               /*3D index*/
}


void Decomposition::printTrueBeadsInDomain(){
    std::cout << "Processor No. " << mpi->rank() << ":" << mybeads << std::endl;
    return;
}

bool Decomposition::calculateCellLength(){

    Ivec prev_num_cells=num_cells;
    num_true_cells={0,0,0};
    num_cells={0,0,0};
    cell_length={0.0, 0.0, 0.0};

    num_true_cells[0]=static_cast<int>(domain_length[0]/rcut);      /*Number of true cells along x*/
    cell_length[0]=domain_length[0]/num_true_cells[0];              /*Length of each cell along x*/
    num_true_cells[1]=static_cast<int>(domain_length[1]/rcut);      /*Number of true cells along y*/
    cell_length[1]=domain_length[1]/num_true_cells[1];              /*Length of each cell along y*/
    num_true_cells[2]=static_cast<int>(domain_length[2]/rcut);      /*Number of true cells along y*/
    cell_length[2]=domain_length[2]/num_true_cells[2];              /*Length of each cell along y*/

    /*Adding ghost cells at both boundaries*/
    num_cells[0]=num_true_cells[0]+2;
    num_cells[1]=num_true_cells[1]+2;
    num_cells[2]=num_true_cells[2]+2;
    prev_totnum_cells=totnum_cells;
    totnum_cells=num_cells[0]*num_cells[1]*num_cells[2];        /*Total number of cells*/
    totnum_true_cells=num_true_cells[0]*num_true_cells[0]*num_true_cells[0];        /*Total number of true cells*/
    totnum_ghost_cells=totnum_cells-totnum_true_cells;          /*The number of ghost cells*/
    
    for(int i=prev_totnum_cells; i<totnum_cells;i++){
        cells.push_back(new Cell(myid, num_domains, i, num_cells));     /*Making a whole cell list*/
    }

    for(int i=prev_totnum_cells-1;i>=totnum_cells;i--){
        delete cells[i];
    }
    cells.resize(totnum_cells);
    return prev_num_cells!=num_cells;
}

void Decomposition::printCellNumberInformation(){
    if(mpi->isMaster()){
        std::cout << " (True) Cells in each domain divided into " << num_true_cells[0] << "*" << num_true_cells[1] << "*" << num_true_cells[2] << " cells" << std::endl;
        std::cout << " Cell size = " << cell_length[0] << "*" << cell_length[1] << "*" << cell_length[2] << std::endl;
    }

    return;
}


void Decomposition::makeCellLists(){
    truecells=Ivec{};
    ghostcells=Ivec{};
    realcells=Ivec{};
    for(int i=0;i<cells.size();i++){
        if(cells[i]->isRealCell())      /*If it's a real cell*/
            realcells.push_back(i);
        if(cells[i]->isGhostCell())     /*If it's a ghost cell*/
            ghostcells.push_back(i);
        else
            truecells.push_back(i);     /*If it's a true (not a ghost) cell*/
    }

    return;
}





void Decomposition::communicateGhostBeads(){
    Ivec2D ptcls_to_send(numprocs, Ivec{});
    Ivec2D cellidx_to_send(numprocs, Ivec{});
    for(int i=0;i<numprocs;i++){
        ptcls_to_send[i].reserve(maxnbeads);
        cellidx_to_send[i].reserve(maxnbeads);
    }

    for(int i=0;i<cells.size();i++){
        if(cells[i]->isRealCell()){     /*If its a real cell*/
            Ivec2D ghostindex=cells[i]->getGhostIndex();        /*Finding cell/domain indices of the ghost of this real cell*/
            Ivec beads=cells[i]->getBeads();                    /*Beads in this cell*/

            for(int j=0;j<beads.size();j++){
                for(int k=0;k<ghostindex.size();k++){
                    int dest=ghostindex[k][0];
                    ptcls_to_send[dest].push_back(beads[j]);        /*Building particle list to send*/
                    cellidx_to_send[dest].push_back(ghostindex[k][1]);  /*Making a list of cells that particles to be sent*/
                }
            }
        }
        else if(cells[i]->isGhostCell()){   /*If it's a ghost cell*/
            Ivec beads=cells[i]->getBeads();            /*Finding beads in a cell*/
            for(int j=0;j<beads.size();j++){
                if(particles[beads[j]]->existsHere()!=TRUEPTCL)     /*if the particle is not here*/
                    particles[beads[j]]->setNone();                 /*Removing information*/
            }
            cells[i]->clear();
        }
    }
    Ivec ghostshere=ptcls_to_send[myid];            
    Ivec ghostcells=cellidx_to_send[myid];
    for(int i=0;i<ghostshere.size();i++){
        cells[ghostcells[i]]->addBeads(ghostshere[i]);      /*Adding beads to ghost cells*/
    }


    for(int i=0;i<numprocs;i++){
        if(myid!=i){
            int sendsize=ptcls_to_send[i].size();
            int recvsize=0;
            MPI_Sendrecv(&sendsize, 1, MPI_INT, i, 0, &recvsize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            /*Serialization of the information for communication*/
            Rvec sendcoord=serializeNCoords(ptcls_to_send[i], particles);
            Rvec sendveloc=serializeNVelocs(ptcls_to_send[i], particles);
            Ivec recvbeads(recvsize);
            Ivec recvcells(recvsize);
            Rvec recvcoord(3*recvsize);
            Rvec recvveloc(3*recvsize);
            /*Sending/receiving particle informations to/from other domains*/
            MPI_Sendrecv(&ptcls_to_send[i][0], sendsize, MPI_INT, i, 1, &recvbeads[0], recvsize, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&cellidx_to_send[i][0], sendsize, MPI_INT, i, 2, &recvcells[0], recvsize, MPI_INT, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&sendcoord[0], 3*sendsize, MPI_DOUBLE, i, 3, &recvcoord[0], 3*recvsize, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&sendveloc[0], 3*sendsize, MPI_DOUBLE, i, 4, &recvveloc[0], 3*recvsize, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            /*Deserialization of the information for communication*/
            deserializeNCoords(recvbeads, particles, recvcoord);
            deserializeNVelocs(recvbeads, particles, recvveloc);
            /*Adding information about new particle received*/
            for(int j=0;j<recvsize;j++){
                cells[recvcells[j]]->addBeads(recvbeads[j]);
                particles[recvbeads[j]]->setGhost();
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    matchGhostsRealBeads();
    assignParticleType();
    return;
}

/*This is the subroutine for refreshing domain beads
 * See the subrouting above*/
void Decomposition::refreshDomainBeads(){
    Ivec2D ptcls_to_send(numprocs, Ivec());
    for(int i=0;i<numprocs;i++){
        ptcls_to_send[i].reserve(maxnbeads);
    }

    for(int i=0;i<totnum_cells;i++){
        if(!cells[i]->isGhostCell()){
            Ivec beads=cells[i]->getBeads();

            int sendsize=beads.size();
            for(int j=0;j<sendsize;j++){
                /*It was found that this function causes segmentation fault with wrong destinationn newindex[0]
                 * See function getNewDomainCellIndex(Particle*)*/
                Ivec newindex=getNewDomainCellIndex(particles[beads[j]]);
                if(newindex[0]!=myid){
                    ptcls_to_send[newindex[0]].push_back(beads[j]);
                    cells[i]->removeBeads(beads[j]);
                    removeBeads(beads[j]);
                    particles[beads[j]]->setNone();
                }
                else{
                    if(newindex[1]!=i){
                        cells[i]->removeBeads(beads[j]);
                        cells[newindex[1]]->addBeads(beads[j]);
                    }
                        
                }
            }
        }
    }
    for(int i=0;i<numprocs;i++){
        if(myid!=i){
            int sendsize=ptcls_to_send[i].size();
            int recvsize=0;
            MPI_Sendrecv(&sendsize, 1, MPI_INT, i, 0, &recvsize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            Rvec sendcoord=serializeNCoords(ptcls_to_send[i], particles);
            Rvec sendveloc=serializeNVelocs(ptcls_to_send[i], particles);
            Ivec recvbeads(recvsize);
            Rvec recvcoord(3*recvsize);
            Rvec recvveloc(3*recvsize);
            MPI_Sendrecv(&ptcls_to_send[i][0], sendsize, MPI_INT, i, 1, &recvbeads[0], recvsize, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&sendcoord[0], 3*sendsize, MPI_DOUBLE, i, 2, &recvcoord[0], 3*recvsize, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&sendveloc[0], 3*sendsize, MPI_DOUBLE, i, 3, &recvveloc[0], 3*recvsize, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            deserializeNCoords(recvbeads, particles, recvcoord);
            deserializeNVelocs(recvbeads, particles, recvveloc);
            for(int j=0;j<recvsize;j++){
                int newcellidx=getNewCellIndex(particles[recvbeads[j]]);
                cells[newcellidx]->addBeads(recvbeads[j]);
                addBeads(recvbeads[j]);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return;
}


void Decomposition::assignParticleType(){
    for(int i=0;i<particles.size();i++){    
        particles[i]->setNone();        /*Changing all particles types into none*/
    }
    for(int i=0;i<totnum_cells;i++){
        if(cells[i]->isGhostCell()){        /*If it is a ghost cell*/
            Ivec beads=cells[i]->getBeads();
            for(int j=0;j<beads.size();j++){
                particles[beads[j]]->setGhost();        /*a particle in this cell is a ghost*/
            }
        }
        else if(!cells[i]->isGhostCell()){      /*If it's noa a ghost cell*/
            Ivec beads=cells[i]->getBeads();

            for(int j=0;j<beads.size();j++){
                particles[beads[j]]->setTrue();     /*a particle is a true particle*/
            }
        }
    }

    return;
}


int Decomposition::getNewCellIndex(Particle* ptcl){
    Ivec cellidx(3,0);
    ptcl->setCoord(pbc.getVectorIntoBox(ptcl->coord));
    cellidx[0]=static_cast<int>(ptcl->getCoord()[0]/cell_length[0]);
    cellidx[1]=static_cast<int>(ptcl->getCoord()[1]/cell_length[1]);
    cellidx[2]=static_cast<int>(ptcl->getCoord()[2]/cell_length[2]);
    cellidx[0]=cellidx[0]%num_true_cells[0]+1;
    cellidx[1]=cellidx[1]%num_true_cells[1]+1;
    cellidx[2]=cellidx[2]%num_true_cells[2]+1;
    return indexing_cell.getIndexFrom3DIndex(cellidx);  /*Returning a cell index of the particle*/
}

int Decomposition::getNewDomainIndex(Particle* ptcl){
    Ivec domainidx(3,0);
    ptcl->setCoord(pbc.getVectorIntoBox(ptcl->coord));
    domainidx[0]=static_cast<int>(ptcl->getCoord()[0]/cell_length[0]);
    domainidx[1]=static_cast<int>(ptcl->getCoord()[1]/cell_length[1]);
    domainidx[2]=static_cast<int>(ptcl->getCoord()[2]/cell_length[2]);
    domainidx[0]=domainidx[0]/num_true_cells[0];
    domainidx[1]=domainidx[1]/num_true_cells[1];
    domainidx[2]=domainidx[2]/num_true_cells[2];
    return indexing_domain.getIndexFrom3DIndex(domainidx);  /*Returning a domain index of the particle*/
}

Ivec Decomposition::getNewDomainCellIndex(Particle* ptcl){
    Ivec cellidx(3,0);
    Ivec domainidx(3,0);
    /*Here, sometimes, new coordinate is not properly calculated. The reason has not been found, but if the particle is 
     * once more brought into the box in the PerioidicBoundary::getVectorIntoBox(const Real3D&), the problem is gone.
     * The reason is stil under investigation*/
    ptcl->setCoord(pbc.getVectorIntoBox(ptcl->coord));
    cellidx[0]=static_cast<int>(ptcl->getCoord()[0]/cell_length[0]);
    cellidx[1]=static_cast<int>(ptcl->getCoord()[1]/cell_length[1]);
    cellidx[2]=static_cast<int>(ptcl->getCoord()[2]/cell_length[2]);
    domainidx[0]=cellidx[0]/num_true_cells[0];
    domainidx[1]=cellidx[1]/num_true_cells[1];
    domainidx[2]=cellidx[2]/num_true_cells[2];
    cellidx[0]=cellidx[0]%num_true_cells[0]+1;
    cellidx[1]=cellidx[1]%num_true_cells[1]+1;
    cellidx[2]=cellidx[2]%num_true_cells[2]+1;
    Ivec result={indexing_domain.getIndexFrom3DIndex(domainidx), indexing_cell.getIndexFrom3DIndex(cellidx)};
    /*Returning a set of cell and domain indices of the particle*/
    return result;
}


void Decomposition::addBeads(int index){
    mybeads.push_back(index);       /*Adding a bead to the domain*/
    num_mybeads++;
    return;
}

void Decomposition::addBeads(Ivec indices){
    for(int i=0;i<indices.size();i++)
        mybeads.push_back(indices[i]);      /*Adding multiple beads to the domain*/
    num_mybeads+=indices.size();
    return;
}
void Decomposition::removeBeads(int index){
    mybeads.erase(std::remove(mybeads.begin(), mybeads.end(), index), mybeads.end());       /*Removing a bead from the domain*/
    num_mybeads--;
    return;
}

void Decomposition::removeBeads(Ivec indices){
    for(int i=0;i<indices.size();i++)
        mybeads.erase(std::remove(mybeads.begin(), mybeads.end(), indices[i]), mybeads.end());  /*Removing multiple beads from the domain*/
    num_mybeads-=indices.size();
    return;
}






void Decomposition::printCellInformation(){
    for(int i=0;i<cells.size();i++){
        cells[i]->printBeads();
    }
    return;
}
/*Matching ghost and real beads information*/
void Decomposition::matchGhostsRealBeads(){
    for(int i=0;i<numprocs;i++){
        ghosts[i].clear();
        realbeads[i].clear();
    }

    for(int i=0;i<cells.size();i++){
        if(cells[i]->isGhostCell()){
            int dest=cells[i]->getRealCellIndex()[0];
            Ivec beads=cells[i]->getBeads();
            for(int j=0;j<beads.size();j++){
                ghosts[dest].push_back(beads[j]);
            }
        }
    }
    for(int i=0;i<numprocs;i++){
        if(myid!=i){
            int sendsize=ghosts[i].size();
            int recvsize;
            MPI_Sendrecv(&sendsize, 1, MPI_INT, i, 0, &recvsize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            realbeads[i]=Ivec(recvsize, 0);
            MPI_Sendrecv(&ghosts[i][0], sendsize, MPI_INT, i, 0, &realbeads[i][0], recvsize, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    return;

}
/*
Ivec2D Decomposition::findProcsOfGhosts(){
    for(int i=0;i<numprocs;i++)
        ghosts[i].clear();
    for(int i=0;i<cells.size();i++){
        if(cells[i]->isGhostCell()){
            int dest=cells[i]->getRealCellIndex()[0];
            Ivec beads=cells[i]->getBeads();
            for(int j=0;j<beads.size();j++){
                ghosts[dest].push_back(beads[j]);
            }
        }
    }
    return ghosts;
}


Ivec2D Decomposition::findProcsOfRealBeads(){
    for(int i=0;i<numprocs;i++)
        realbeads[i].clear();
    for(int i=0;i<cells.size();i++){
        if(cells[i]->isRealCell()){
            Ivec2D ghost_cells=cells[i]->getGhostIndex();
            Ivec beads=cells[i]->getBeads();
            for(int j=0;j<beads.size();j++){
                for(int k=0;k<ghost_cells.size();k++){
                    int source=ghost_cells[k][0];
                    realbeads[source].push_back(beads[j]);
                }
            }
        }
    }
    return realbeads;

}
 */   

/*Setting all information is not synced*/
void Decomposition::setUnsynced(){
    is_coord_synced=false;
    is_veloc_synced=false;
    is_force_synced=false;
    is_density_synced=false;
    is_stress_synced=false;
    return;
}

/*Getting information of unfrozen beads*/
Ivec Decomposition::getUnfrozenBeads(){
    int num_mybeads=mybeads.size();
    Ivec myufbeads={};
    for(int i=0;i<num_mybeads;i++){
        if(!particles[mybeads[i]]->isFrozen())
            myufbeads.push_back(mybeads[i]);
    } 
    return myufbeads;
}


/*Gathering positions of true particles all over the domains*/
void Decomposition::gatherCoords(){
    if(!is_coord_synced){
        if(mpi->isMaster()){
            for(int i=1;i<numprocs;i++){
                int recvsize;
                MPI_Recv(&recvsize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                Ivec recvbeads(recvsize);
                MPI_Recv(&recvbeads[0], recvsize, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                Rvec commcoord(3*recvsize);
                MPI_Recv(&commcoord[0], 3*recvsize, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                deserializeNCoords(recvbeads, particles, commcoord);
            }
        }
        else{
            Ivec sendbeads;
            int sendnum;
            if(!dumpfrozen){
                sendbeads=getUnfrozenBeads();
                sendnum=sendbeads.size();
            }
            else{
                sendbeads=mybeads;
                sendnum=mybeads.size();
            }

            MPI_Send(&sendnum, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
            MPI_Send(&sendbeads[0], sendnum, MPI_INT, MASTER, 1, MPI_COMM_WORLD);
            Rvec commcoord=serializeNCoords(sendbeads, particles);
            MPI_Send(&commcoord[0], 3*sendnum, MPI_DOUBLE, MASTER, 2, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    is_coord_synced=true;

    return;
}

/*Gathering velocities of true particles all over the domains*/
void Decomposition::gatherVelocs(){
    if(!is_veloc_synced){

        if(mpi->isMaster()){
            for(int i=1;i<numprocs;i++){
                int recvsize;
                MPI_Recv(&recvsize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                Ivec recvbeads(recvsize);
                MPI_Recv(&recvbeads[0], recvsize, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                Rvec commveloc(3*recvsize);
                MPI_Recv(&commveloc[0], 3*recvsize, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                deserializeNVelocs(recvbeads, particles, commveloc);
            }
        }
        else{
            Ivec sendbeads;
            int sendnum;
            if(!dumpfrozen){
                sendbeads=getUnfrozenBeads();
                sendnum=sendbeads.size();
            }
            else{
                sendbeads=mybeads;
                sendnum=mybeads.size();
            }

            MPI_Send(&sendnum, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
            MPI_Send(&sendbeads[0], sendnum, MPI_INT, MASTER, 1, MPI_COMM_WORLD);
            Rvec commveloc=serializeNVelocs(sendbeads, particles);
            MPI_Send(&commveloc[0], 3*sendnum, MPI_DOUBLE, MASTER, 3, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    is_veloc_synced=true;
    return;
}

/*Gathering forces of true particles all over the domains*/
void Decomposition::gatherForces(){
    if(!is_force_synced){
        if(mpi->isMaster()){
            for(int i=1;i<numprocs;i++){
                int recvsize;
                MPI_Recv(&recvsize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                Ivec recvbeads(recvsize);
                MPI_Recv(&recvbeads[0], recvsize, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                Rvec commforce(3*recvsize);
                MPI_Recv(&commforce[0], 3*recvsize, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                deserializeNForces(recvbeads, particles, commforce);
            }
        }
        else{
            Ivec sendbeads;
            int sendnum;
            if(!dumpfrozen){
                sendbeads=getUnfrozenBeads();
                sendnum=sendbeads.size();
            }
            else{
                sendbeads=mybeads;
                sendnum=mybeads.size();
            }

            MPI_Send(&sendnum, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
            MPI_Send(&sendbeads[0], sendnum, MPI_INT, MASTER, 1, MPI_COMM_WORLD);
            Rvec commforce=serializeNForces(sendbeads, particles);
            MPI_Send(&commforce[0], 3*sendnum, MPI_DOUBLE, MASTER, 3, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    is_force_synced=true;
    return;
}


/*Gathering local densities of true particles all over the domains*/
void Decomposition::gatherDensities(){
    if(!is_density_synced){
        if(mpi->isMaster()){
            for(int i=1;i<numprocs;i++){
                int recvsize;
                MPI_Recv(&recvsize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                Ivec recvbeads(recvsize);
                MPI_Recv(&recvbeads[0], recvsize, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                Rvec commdensity(recvsize);
                MPI_Recv(&commdensity[0], recvsize, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for(int j=0;j<recvsize;j++)
                    particles[recvbeads[j]]->density=commdensity[j];
            }
        }
        else{
            Ivec sendbeads;
            int sendnum;
            if(!dumpfrozen){
                sendbeads=getUnfrozenBeads();
                sendnum=sendbeads.size();
            }
            else{
                sendbeads=mybeads;
                sendnum=mybeads.size();
            }

            MPI_Send(&sendnum, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
            MPI_Send(&sendbeads[0], sendnum, MPI_INT, MASTER, 1, MPI_COMM_WORLD);
            Rvec commdensity(sendnum);
            for(int j=0;j<num_mybeads;j++)
                commdensity[j]=particles[sendbeads[j]]->density;
            MPI_Send(&commdensity[0], sendnum, MPI_DOUBLE, MASTER, 3, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        is_density_synced=true;
    }

    return;
}

/*Gathering particle stresses of true particles all over the domains*/
void Decomposition::gatherStresses(){
    if(!is_stress_synced){
        if(mpi->isMaster()){
            for(int i=1;i<numprocs;i++){
                int recvsize;
                MPI_Recv(&recvsize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                Ivec recvbeads(recvsize);
                MPI_Recv(&recvbeads[0], recvsize, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                Rvec commstress(9*recvsize);
                MPI_Recv(&commstress[0], 9*recvsize, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                deserializeNStresses(recvbeads, particles, commstress);
            }
        }
        else{
            Ivec sendbeads;
            int sendnum;
            if(!dumpfrozen){
                sendbeads=getUnfrozenBeads();
                sendnum=sendbeads.size();
            }
            else{
                sendbeads=mybeads;
                sendnum=mybeads.size();
            }

            MPI_Send(&sendnum, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
            MPI_Send(&sendbeads[0], sendnum, MPI_INT, MASTER, 1, MPI_COMM_WORLD);
            Rvec commstress=serializeNStresses(sendbeads, particles);
            MPI_Send(&commstress[0], 9*sendnum, MPI_DOUBLE, MASTER, 2, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    is_stress_synced=true;
    return;
}

/*Gathering all information*/
void Decomposition::gatherAll(){
    gatherCoords();
    gatherVelocs();
    gatherForces();
    gatherDensities();
    gatherStresses();
    return;
}



/*Resetting box information*/
void Decomposition::resetBox(Real3D newbox){
    box=newbox;
    config->setBox(box);

    pbc.resetBoxSize(box);
    calculateDomainLength();
//    printDomainLengthInformation();
    bool resetcell=calculateCellLength();
    if(resetcell){
        for(int i=0;i<cells.size();i++){
            cells[i]->resetCell(i, num_cells);
        }
    //    printCellNumberInformation();
        indexing_cell=Indexing(num_cells);
        makeCellLists();
    //    if(mpi->isMaster()) printCellInformation();
        refreshBeadsWhenChangingCellNumber();
    }
}

/*Refreshing beads when the cell number in a domain changes*/
void Decomposition::refreshBeadsWhenChangingCellNumber(){
    Ivec2D ptcls_to_send(numprocs, Ivec());
    for(int i=0;i<numprocs;i++){
        ptcls_to_send[i].reserve(maxnbeads);
    }
    
    Ivec ori_beads=mybeads;
    for(int i=0;i<ori_beads.size();i++){
        Particle* ptcl=particles[ori_beads[i]];
        Ivec newindex=getNewDomainCellIndex(ptcl);
        if(newindex[0]!=myid){
            ptcls_to_send[newindex[0]].push_back(ori_beads[i]);
            removeBeads(ori_beads[i]);
            particles[ori_beads[i]]->setNone();
        }
        else{
            cells[newindex[1]]->addBeads(ori_beads[i]);
        }
    }
    for(int i=0;i<numprocs;i++){
        if(myid!=i){
            int sendsize=ptcls_to_send[i].size();
            int recvsize=0;
            MPI_Sendrecv(&sendsize, 1, MPI_INT, i, 0, &recvsize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            Rvec sendcoord=serializeNCoords(ptcls_to_send[i], particles);
            Rvec sendveloc=serializeNVelocs(ptcls_to_send[i], particles);
            Ivec recvbeads(recvsize);
            Rvec recvcoord(3*recvsize);
            Rvec recvveloc(3*recvsize);
            MPI_Sendrecv(&ptcls_to_send[i][0], sendsize, MPI_INT, i, 1, &recvbeads[0], recvsize, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&sendcoord[0], 3*sendsize, MPI_DOUBLE, i, 2, &recvcoord[0], 3*recvsize, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&sendveloc[0], 3*sendsize, MPI_DOUBLE, i, 3, &recvveloc[0], 3*recvsize, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            deserializeNCoords(recvbeads, particles, recvcoord);
            deserializeNVelocs(recvbeads, particles, recvveloc);
            for(int j=0;j<recvsize;j++){
                int newcellidx=getNewCellIndex(particles[recvbeads[j]]);
                cells[newcellidx]->addBeads(recvbeads[j]);
                addBeads(recvbeads[j]);
                particles[recvbeads[j]]->setTrue();
            }
        }
    }
    return;
}


/*Finding where is the true image of the particle*/
void Decomposition::whereIsTheBead(int index){
    /*
    for(int i=0;i<mybeads.size();i++){
        if(particles[mybeads[i]]->getParticleIndex()==index){
            std::cout << "Particle index " << index << " is at " << particles[mybeads[i]]->coord << ", " ;
            std::cout << "Domain ID = " << myid << " (" << my3did << ")" << ", ";
            std::cout << "Cell ID ="  << getNewCellIndex(particles[mybeads[i]]) << " (" <<   indexing_cell.get3DIndexFromIndex(getNewCellIndex(particles[mybeads[i]])) << ")" << std::endl;
        }
   }
   */
    for(int i=0;i<cells.size();i++){
        if(!cells[i]->isGhostCell()){
            Ivec beads=cells[i]->getBeads();
            for(int j=0;j<beads.size();j++){
                if(particles[beads[j]]->getParticleIndex()==index){
                    std::cout << "Particle index " << index << " is at " << particles[beads[j]]->coord << ", " ;
                    std::cout << "Domain ID = " << myid << " (" << my3did << ")" << ", ";
                    std::cout << "Cell ID ="  << cells[i]->getMyCellIndex() << " (" <<   cells[i]->getMyCell3DIndex() << std::endl;
                }
            }
        }
    }
    return;
}


/*Setting the existence of the particle*/
void Decomposition::setExistence(){
    for(int i=0;i<particles.size();i++){
        particles[i]->setNone();
    }
    for(int i=0;i<cells.size();i++){
        if(cells[i]->isGhostCell()){
            Ivec beads=cells[i]->getBeads();
            for(int j=0;j<beads.size();j++){
                particles[beads[j]]->setGhost();
            }
        }
        else if(cells[i]->isRealCell()){
            Ivec beads=cells[i]->getBeads();
            for(int j=0;j<beads.size();j++){
                particles[beads[j]]->setTrue();
            }
        }
    }
    return;
}
    

