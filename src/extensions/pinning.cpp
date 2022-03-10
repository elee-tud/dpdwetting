#include "pinning.hpp"
#include <stdlib.h>
#include "../selection.hpp"


using namespace dpd;
Pinning::Pinning(Topology* topol, Configuration* config, Decomposition* decomp):Extension(topol, config, decomp){
    need_prev_position=true;
    need_prev_velocity=false;
    for_position=true;
    for_velocity=true;
    direct=control->getPinningDirection();
    if(control->getPinningRandom()){
        molname=control->getPinningMoleculeName();
        npinpermol=control->getNumPinningPerMolecule();
        Selection select(particles);
        molecules=select.selectMoleculesByMoleculeName(molname);
        pinindex.clear();
        for(int i=0;i<molecules.size();i++){
            int nptclsinmol=molecules[i].size();
            int ridx=rand()%nptclsinmol;
            pinindex.push_back(molecules[i][ridx]->getParticleIndex()-1);
        }
        int msize=pinindex.size();
        MPI_Bcast(&msize, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
        if(!mpi->isMaster())
            pinindex=Ivec(msize, 0);
        MPI_Bcast(&pinindex[0], msize, MPI_INT, MASTER, MPI_COMM_WORLD);

            


    }
    else{
        pinindex=control->getPinningIndex();
    }

    for(int i=0;i<pinindex.size();i++){
        particles[pinindex[i]]->setPinned();
    }
    /*
    if(mpi->rank()==1){
        for(int i=0;i<particles.size();i++){
            if(particles[i]->isPinned())
                std::cout << particles[i]->getParticleIndex() << std::endl;
        }
    }
    */


}

void Pinning::applyExtensionForPosition(int step){
    Ivec beads=decomp->getBeadsIndexInDomain();
    for(int i=0;i<beads.size();i++){
        if(particles[beads[i]]->isPinned()){
            if(direct[0]==1)
                particles[beads[i]]->coord[0]=particles[beads[i]]->prevcoord[0];
            if(direct[1]==1)
                particles[beads[i]]->coord[1]=particles[beads[i]]->prevcoord[1];
            if(direct[2]==1)
                particles[beads[i]]->coord[2]=particles[beads[i]]->prevcoord[2];
        }
    }

    return;
}



void Pinning::applyExtensionForVelocity(int step){
    Ivec beads=decomp->getBeadsIndexInDomain();
    for(int i=0;i<beads.size();i++){
        if(particles[beads[i]]->isPinned()){
            if(direct[0]==1)
                particles[beads[i]]->veloc[0]=0.;
            if(direct[1]==1)
                particles[beads[i]]->veloc[1]=0.;
            if(direct[2]==1)
                particles[beads[i]]->veloc[2]=0.;
        }
    }

    return;
}
