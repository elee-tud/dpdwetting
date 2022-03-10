#include "externalforce.hpp"

using namespace dpd;
ExternalForce::ExternalForce(Topology* topol, Configuration* config, Decomposition* decomp):topol(topol), config(config), decomp(decomp){
    cells=decomp->getCells();
    particles=config->getParticles();
    box=config->getBox();
    pbc=PeriodicBoundary(box);
    control=decomp->getControl();
    mpi=decomp->getMPI();
}

void ExternalForce::calculateForce(){
    int numcells=cells.size();
    for(int i=0;i<numcells;i++){
        if(!cells[i]->isGhostCell()){
            mybeads=cells[i]->getBeads();
            for(int j=0;j<mybeads.size();j++){
                calculateSingleForce(particles[mybeads[j]]);
            }
        }
    }
    return;
}
    

