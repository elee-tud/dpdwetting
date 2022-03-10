#include "boxelongation.hpp"
using namespace dpd;
BoxElongation::BoxElongation(Topology *topol, Configuration* config, Decomposition* decomp):Extension(topol, config, decomp){
    need_prev_position=false;
    need_prev_velocity=false;
    for_position=true;
    for_velocity=false;
    freq=control->getDeformationFrequency();
    factor=control->getDeformationFactor();
    invsqrtf=1./sqrt(factor);
}

void BoxElongation::applyExtensionForPosition(int step){
    if(step%freq==0){
        Real3D newbox=decomp->getBox();
        newbox[0]*=factor;
        newbox[1]*=invsqrtf;
        newbox[2]*=invsqrtf;
        Ivec beads=decomp->getBeadsIndexInDomain();
        for(int j=0;j<beads.size();j++){
            particles[beads[j]]->coord[0]*=factor;
            particles[beads[j]]->coord[1]*=invsqrtf;
            particles[beads[j]]->coord[2]*=invsqrtf;
        }
        decomp->resetBox(newbox);
    }


    return;
}

