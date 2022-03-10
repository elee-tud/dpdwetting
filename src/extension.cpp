#include "extension.hpp"


using namespace dpd;
Extension::Extension(Topology* topol, Configuration* config, Decomposition* decomp):topol(topol), config(config), decomp(decomp){
    control=decomp->getControl();
    mpi=decomp->getMPI();
    cells=decomp->getCells();
    particles=decomp->getParticles();



}
