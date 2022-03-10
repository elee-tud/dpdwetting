#include "dropvelocity.hpp"
using namespace dpd;
DropVelocity::DropVelocity(InitialSet initset):Property(initset){
    title="  Calculation of radial and axial velocities of a droplet";
    
    nprops=2;
    for_timeevol=true;
    need_position=true;
    need_velocity=true;

    outfile="drop_velocity.out";
    
    outtitle_tevol="#Droplet velocity along axial and radial directions";
    outheader_tevol=Svec{"Time", "v_axi", "v_rad"};

}

void DropVelocity::getSpecificParameters(){
    dz=0.65;

    command->getCommandSingleOption("-dz", dz, &dz);
    surfaceb=control->getWallMinPosition()[2];
    surfacet=control->getWallMaxPosition()[2];
    return;
}

void DropVelocity::initializeVariables(){
    numz=static_cast<int>((surfacet-surfaceb)/dz)+1;
    for(int i=0;i<numz;i++)
        particles_atz.push_back(new ParticleGroup(pbc));
    com=R3vec(numz, Real3D(0.0));
    ref0=R3vec(numz, Real3D(box/2-1.0));
    ref1=R3vec(numz, Real3D(box/2+1.0));
    for(int i=0;i<numz;i++){
        ref0[i][2]=surfaceb+dz*0.5+i*dz;
        ref1[i][2]=ref0[i][2]+dz;
    }

    initializeResultArrays();
    return;
}

void DropVelocity::calculateStep(int step){
    real vaxial=0.;
    real vradial=0.;
    for(int i=0;i<numz;i++){
        particles_atz[i]->clear();
    }
    for(int i=0;i<nbeads;i++){
        int index=static_cast<int>((particles[i]->coord[2]-surfaceb-0.5*dz)/dz);
        particles_atz[index]->addParticle(particles[i]);
        vaxial+=particles[i]->veloc[2];
    }
    vaxial/=nbeads;


    for(int i=0;i<numz;i++){
        int numptcls=particles_atz[i]->size();
        if(numptcls>1){
            particles_atz[i]->setReference(ref0[i], ref1[i]);
            com[i]=particles_atz[i]->calculateCenterOfMass();
            for(int j=0;j<numptcls;j++){
                Real3D dvec=pbc.getMinimumImageVector(particles_atz[i]->at(j)->coord, com[i]).unit();
                real vr=particles_atz[i]->at(j)->veloc[0]*dvec[0]+particles_atz[i]->at(j)->veloc[1]*dvec[1];
                vradial+=vr;

            }
        }
    }
    vradial/=nbeads;
    tevol[0][step]=vaxial;
    tevol[1][step]=vradial;
    return;
}

    
