#include "solidwall.hpp"


using namespace dpd;
SolidWall::SolidWall(Topology* topol, Configuration* config, Decomposition* decomp):Extension(topol, config, decomp){
    need_prev_position=true;
    need_prev_velocity=true;
    for_position=true;
    for_velocity=true;
    direct=control->getWallDirection();
    dmin=control->getWallMinPosition();
    dmax=control->getWallMaxPosition();
}

void SolidWall::applyExtensionForPosition(int step){
    for(int i=0;i<cells.size();i++){
        Ivec beads=cells[i]->getBeads();
        for(int j=0;j<beads.size();j++){
            if(!particles[beads[j]]->isFrozen()){
                if(isCrossingWall(particles[beads[j]])){
                    particles[beads[j]]->coord=particles[beads[j]]->prevcoord;
                }
            }
        }
    }
    return;
}

bool SolidWall::isCrossingWall(Particle *ptcl){
    if(direct[0]==1){
        if(ptcl->coord[0]<dmin[0] || ptcl->coord[0]>=dmax[0])
            return true;
    }
    if(direct[1]==1){
        if(ptcl->coord[1]<dmin[1] || ptcl->coord[1]>=dmax[1])
            return true;
    }
    if(direct[2]==1){
        if(ptcl->coord[2]<dmin[2] || ptcl->coord[2]>=dmax[2])
            return true;
    }
    return false;
}



void SolidWall::applyExtensionForVelocity(int step){
    for(int i=0;i<cells.size();i++){
        Ivec beads=cells[i]->getBeads();
        for(int j=0;j<beads.size();j++){
            if(particles[beads[j]]->coord==particles[beads[j]]->prevcoord)
                particles[beads[j]]->veloc=particles[beads[j]]->prevveloc;
        }
    }
    return;
}
