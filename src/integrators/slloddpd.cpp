#include "slloddpd.hpp"

using namespace dpd;
SllodDPD::SllodDPD(Initialization *init):Integrator(init){
    shearten=control->getShearTensor();
    halfbox=decomp->getBox()/2;
    refpos=Real3D{0.0, 0.0, 0.0};
    refpos[0]=control->getWallMinPosition()[0];
    refpos[1]=control->getWallMinPosition()[1];
    refpos[2]=control->getWallMinPosition()[2];


}
void SllodDPD::updatePosition(bool save_prev){
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
void SllodDPD::updatePosition(){
    mybeads=decomp->getBeadsIndexInDomain();
    halfbox=decomp->getBox()/2;


    for(int i=0;i<mybeads.size();i++){
        if(!particles[mybeads[i]]->isFrozen()){
            Real3D shear_vel(0.0);
            /*
            shear_vel[0]=shearten[0]*(particles[mybeads[i]]->coord[0]-halfbox[0])+shearten[1]*(particles[mybeads[i]]->coord[1]-halfbox[1])+shearten[2]*(particles[mybeads[i]]->coord[2]-halfbox[2]);
            shear_vel[1]=shearten[3]*(particles[mybeads[i]]->coord[0]-halfbox[0])+shearten[4]*(particles[mybeads[i]]->coord[1]-halfbox[1])+shearten[5]*(particles[mybeads[i]]->coord[2]-halfbox[2]);
            shear_vel[2]=shearten[6]*(particles[mybeads[i]]->coord[0]-halfbox[0])+shearten[7]*(particles[mybeads[i]]->coord[1]-halfbox[1])+shearten[8]*(particles[mybeads[i]]->coord[2]-halfbox[2]);
            */
            shear_vel[0]=shearten[0]*(particles[mybeads[i]]->coord[0]-refpos[0])+shearten[1]*(particles[mybeads[i]]->coord[1]-refpos[1])+shearten[2]*(particles[mybeads[i]]->coord[2]-refpos[2]);
            shear_vel[1]=shearten[3]*(particles[mybeads[i]]->coord[0]-refpos[0])+shearten[4]*(particles[mybeads[i]]->coord[1]-refpos[1])+shearten[5]*(particles[mybeads[i]]->coord[2]-refpos[2]);
            shear_vel[2]=shearten[6]*(particles[mybeads[i]]->coord[0]-refpos[0])+shearten[7]*(particles[mybeads[i]]->coord[1]-refpos[1])+shearten[8]*(particles[mybeads[i]]->coord[2]-refpos[2]);
            particles[mybeads[i]]->veloc+=shear_vel;
            Real3D shear_frc(0.0);
            shear_frc[0]=shearten[0]*particles[mybeads[i]]->veloc[0]+shearten[1]*particles[mybeads[i]]->veloc[1]+shearten[2]*particles[mybeads[i]]->veloc[2];
            shear_frc[1]=shearten[3]*particles[mybeads[i]]->veloc[0]+shearten[4]*particles[mybeads[i]]->veloc[1]+shearten[5]*particles[mybeads[i]]->veloc[2];
            shear_frc[2]=shearten[6]*particles[mybeads[i]]->veloc[0]+shearten[7]*particles[mybeads[i]]->veloc[1]+shearten[8]*particles[mybeads[i]]->veloc[2];
            particles[mybeads[i]]->force+=shear_frc*particles[mybeads[i]]->getMass();;

            particles[mybeads[i]]->coord+=dt*particles[mybeads[i]]->veloc+hdtsqr*particles[mybeads[i]]->force/particles[mybeads[i]]->getMass();
        }
        
    }
    return;
}

void SllodDPD::updateVelocity(bool save_prev){
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

void SllodDPD::updateVelocity(){
    mybeads=decomp->getBeadsIndexInDomain();
    for(int i=0;i<mybeads.size();i++){
        if(!particles[mybeads[i]]->isFrozen()){
            Real3D shear_frc(0.0);
            shear_frc[0]=shearten[0]*particles[mybeads[i]]->veloc[0]+shearten[1]*particles[mybeads[i]]->veloc[1]+shearten[2]*particles[mybeads[i]]->veloc[2];
            shear_frc[1]=shearten[3]*particles[mybeads[i]]->veloc[0]+shearten[4]*particles[mybeads[i]]->veloc[1]+shearten[5]*particles[mybeads[i]]->veloc[2];
            shear_frc[2]=shearten[6]*particles[mybeads[i]]->veloc[0]+shearten[7]*particles[mybeads[i]]->veloc[1]+shearten[8]*particles[mybeads[i]]->veloc[2];
            particles[mybeads[i]]->force+=shear_frc*particles[mybeads[i]]->getMass();;
            particles[mybeads[i]]->veloc+=0.5*dt*particles[mybeads[i]]->force/particles[mybeads[i]]->getMass();
        }
//        std::cout << particles[mybeads[i]]->veloc << std::endl;
    }
    return;
}




