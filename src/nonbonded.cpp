#include "nonbonded.hpp"


using namespace dpd;
Nonbonded::Nonbonded(Topology* topol, Configuration* config, Decomposition* decomp):topol(topol), config(config), decomp(decomp){
    cells=decomp->getCells();
    truecells=decomp->getTrueCells();
    ghostcells=decomp->getGhostCells();
    particles=config->getParticles();
    box=config->getBox();
    pbc=PeriodicBoundary(box);
    control=decomp->getControl();
    mpi=decomp->getMPI();
    integrator=control->getIntegrator();
}

void Nonbonded::calculateForce(){
    cells=decomp->getCells();
    truecells=decomp->getTrueCells();
    ghostcells=decomp->getGhostCells();
    pbc=PeriodicBoundary(decomp->getBox());
    for(int i=0;i<truecells.size();i++){
        Ivec nbidx=cells[truecells[i]]->getNeighborCells();
        calculateForceInTheCells(cells[truecells[i]]);
        for(int j=0;j<nbidx.size();j++){
            calculateForceBetweenCells(cells[truecells[i]], cells[nbidx[j]]);
        }
    }
    return;
}


void Nonbonded::calculateForceInTheCells(Cell* cell){
    Ivec mybeads=cell->getBeads();
    int num_mybeads=mybeads.size();
    for(int i=0;i<num_mybeads;i++){
        for(int j=i+1;j<num_mybeads;j++){
            calculatePairForce(mybeads[i], mybeads[j]);
        }
    }
    return;
}


void Nonbonded::calculateForceBetweenCells(Cell* cell1, Cell* cell2){
    Ivec mybeads1=cell1->getBeads();
    Ivec mybeads2=cell2->getBeads();
    int num_mybeads1=mybeads1.size();
    int num_mybeads2=mybeads2.size();
    for(int i=0;i<num_mybeads1;i++){
        for(int j=0;j<num_mybeads2;j++){
//            std::cout << "Cell=("<< topol->getMPI()->rank() << "," << cell1->getMyCell3DIndex() << "," << cell2->getMyCell3DIndex() << ")";
            calculatePairForce(mybeads1[i], mybeads2[j]);
        }
    }
    return;
}



void Nonbonded::calculateForce(Real3D com, real** press, real dr){
    cells=decomp->getCells();
    truecells=decomp->getTrueCells();
    ghostcells=decomp->getGhostCells();
    pbc=PeriodicBoundary(decomp->getBox());
    for(int i=0;i<truecells.size();i++){
        Ivec nbidx=cells[truecells[i]]->getNeighborCells();
        calculateForceInTheCells(cells[truecells[i]], com, press, dr);
        for(int j=0;j<nbidx.size();j++){
            calculateForceBetweenCells(cells[truecells[i]], cells[nbidx[j]], com, press, dr);
        }
    }
    return;
}

void Nonbonded::calculateForceInTheCells(Cell* cell, Real3D com, real** press, real dr){
    Ivec mybeads=cell->getBeads();
    int num_mybeads=mybeads.size();
    for(int i=0;i<num_mybeads;i++){
        for(int j=i+1;j<num_mybeads;j++){
            calculatePairForce(mybeads[i], mybeads[j], com, press, dr);
        }
    }
    return;
}

void Nonbonded::calculateForceBetweenCells(Cell* cell1, Cell* cell2, Real3D com, real** press, real dr){
    Ivec mybeads1=cell1->getBeads();
    Ivec mybeads2=cell2->getBeads();
    int num_mybeads1=mybeads1.size();
    int num_mybeads2=mybeads2.size();
    for(int i=0;i<num_mybeads1;i++){
        for(int j=0;j<num_mybeads2;j++){
//            std::cout << "Cell=("<< topol->getMPI()->rank() << "," << cell1->getMyCell3DIndex() << "," << cell2->getMyCell3DIndex() << ")";
            calculatePairForce(mybeads1[i], mybeads2[j], com, press, dr);
        }
    }
    return;
}

