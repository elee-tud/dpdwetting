#include "uniflowpbc.hpp"
using namespace dpd;
UniflowPBC::UniflowPBC(Topology *topol, Configuration* config, Decomposition* decomp):Extension(topol, config, decomp){
    need_prev_position=false;
    need_prev_velocity=false;
    for_position=true;
    for_velocity=false;
    oldbox=decomp->getBox();
    sheartensor=control->getShearTensor();
    dt=control->getTimeStep();

}

void UniflowPBC::applyExtensionForPosition(int step){
    Real3D newbox;
    newbox[0]=oldbox[0]*exp(sheartensor[0]*dt*step);
    newbox[1]=oldbox[1]*exp(sheartensor[4]*dt*step);
    newbox[2]=oldbox[2]*exp(sheartensor[8]*dt*step);
    decomp->resetBox(newbox);

    return;
}

