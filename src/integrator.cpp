#include "integrator.hpp"

using namespace dpd;
Integrator::Integrator(Initialization* init):init(init){
    decomp=init->getDecomposition();
    mpi=init->getMPI();
    control=init->getControl();
    dt=control->getTimeStep();
    hdtsqr=dt*dt/2.0;
    particles=decomp->getParticles();
    mybeads=decomp->getBeadsIndexInDomain();

}
