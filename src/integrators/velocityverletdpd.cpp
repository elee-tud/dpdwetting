#include "velocityverletdpd.hpp"

using namespace dpd;
VelocityVerletDPD::VelocityVerletDPD(Initialization *init):Integrator(init){


}
void VelocityVerletDPD::updatePosition(bool save_prev){
    if(save_prev){
        mybeads=decomp->getBeadsIndexInDomain();
        for(int i=0;i<mybeads.size();i++){
            if(!particles[mybeads[i]]->isFrozen())
                particles[mybeads[i]]->prevcoord=particles[mybeads[i]]->coord;
        }
    }
    updatePosition();
    return;
}
void VelocityVerletDPD::updatePosition(){
    mybeads=decomp->getBeadsIndexInDomain();
    for(int i=0;i<mybeads.size();i++){
        if(!particles[mybeads[i]]->isFrozen())
            particles[mybeads[i]]->coord+=dt*particles[mybeads[i]]->veloc+hdtsqr*particles[mybeads[i]]->force/particles[mybeads[i]]->getMass();
        
    }
    return;
}

void VelocityVerletDPD::updateVelocity(bool save_prev){
    if(save_prev){
        mybeads=decomp->getBeadsIndexInDomain();
        for(int i=0;i<mybeads.size();i++){
            if(!particles[mybeads[i]]->isFrozen())
                particles[mybeads[i]]->prevveloc=particles[mybeads[i]]->veloc;
        }
    }
    updateVelocity();
    return;
}

void VelocityVerletDPD::updateVelocity(){
    mybeads=decomp->getBeadsIndexInDomain();
    for(int i=0;i<mybeads.size();i++){
        if(!particles[mybeads[i]]->isFrozen())
            particles[mybeads[i]]->veloc+=0.5*dt*particles[mybeads[i]]->force/particles[mybeads[i]]->getMass();
//        std::cout << particles[mybeads[i]]->veloc << std::endl;
    }
    return;
}




