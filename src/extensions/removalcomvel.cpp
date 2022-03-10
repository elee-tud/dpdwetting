#include "removalcomvel.hpp"

using namespace dpd;
RemovalCOMVel::RemovalCOMVel(Topology* topol, Configuration* config, Decomposition* decomp):Extension(topol, config, decomp){
    need_prev_position=false;
    need_prev_position=false;
    for_position=true;
    for_velocity=false;
    rmcomv_freq=control->getRmCOMVelFrequency();
    rmcomv_dir=control->getRmCOMVelDirection();
}

void RemovalCOMVel::applyExtensionForVelocity(int step){
    if(step%rmcomv_freq==0){

        comvel=Real3D{0.0, 0.0, 0.0};
        comvelproc=Real3D{0.0, 0.0, 0.0};
        Ivec mybeads=decomp->getBeadsIndexInDomain();
        for(int i=0;i<mybeads.size();i++){
            comvelproc+=particles[mybeads[i]]->veloc*particles[mybeads[i]]->getMass();
        }
        MPI_Allreduce(&comvelproc[0], &comvel[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        comvel/=topol->getTotalMass();
        for(int i=0;i<mybeads.size();i++){
            if(!particles[mybeads[i]]->isFrozen()){
                if(rmcomv_dir[0])
                    particles[mybeads[i]]->veloc[0]-=comvel[0];
                if(rmcomv_dir[1])
                    particles[mybeads[i]]->veloc[1]-=comvel[1];
                if(rmcomv_dir[2])
                    particles[mybeads[i]]->veloc[2]-=comvel[2];
            }
        }
    }
    return;
}

    
