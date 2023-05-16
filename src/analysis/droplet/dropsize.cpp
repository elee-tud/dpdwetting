#include "dropsize.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
DropSize::DropSize(InitialSet initset):Property(initset){
    title="   Calculation of a droplet size, base diameter, and a contact angle." ;
    if(control->getWallMinPosition()==Rvec(3, -10.0))
        with_wall=false;
    else
        with_wall=true;
    if(with_wall){
        nprops=5;
        outheader_tevol=Svec{"Time", "Radius", "Base_radius", "Second_radius", "contact_angle", "Height"};
    }
    else{
        nprops=2;
        outheader_tevol=Svec{"Time", "Radius", "Height"};
    }

    for_timeevol=true;
    need_position=true;

    outfile="dropsize_linear.out";


    outtitle_tevol="#Droplet size, base radius, contact angle as a function of time";


}



void DropSize::getSpecificParameters(){
    dz=0.65;
    dcrit=3.04;
    dr=0.5;
    hradius=1.5;
    pilheight=0.0;

    command->getCommandSingleOption("-dz", dz, &dz);
    command->getCommandSingleOption("-dd", dcrit, &dcrit);
    command->getCommandSingleOption("-dr", dr, &dr);
    command->getCommandSingleOption("-hr", hradius, &hradius);
    command->getCommandSingleOption("-hr", hradius, &hradius);
    command->getCommandSingleOption("-ph", pilheight, &pilheight);
    if(with_wall){
        surfaceb=control->getWallMinPosition()[2]+pilheight;
        surfacet=control->getWallMaxPosition()[2]-pilheight;
    }
    else{
        surfaceb=0.;
        surfacet=box[2];
    }
    liquidgrps=control->getLiquidGroups();
    return;
}


void DropSize::initializeVariables(){
    if(liquidgrps.size()==0){
        for(int i=0;i<nbeads;i++){
            liquididx.push_back(i);
        }
    }
    else{
        for(int i=0;i<nbeads;i++){
            if(std::find(liquidgrps.begin(), liquidgrps.end(), particles[i]->getMoleculeName())!=liquidgrps.end())
                liquididx.push_back(i);
        }
    }
    nliqptcls=liquididx.size();

    
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
    radius=Rvec(numz);

    numr=static_cast<int>(pow(box[1]*box[1]+box[0]*box[0], 0.5)*1.4142135/dr)+1;
    rdensity=Rvec(numr, 0.);
    initializeResultArrays();
    return;
}

void DropSize::calculateStep(int step){
    for(int i=0;i<numz;i++){
        particles_atz[i]->clear();
    }
    Rvec comxy{0.0, 0.0};
    for(int i=0;i<nliqptcls;i++){
        comxy[0]+=particles[liquididx[i]]->coord[0];
        comxy[1]+=particles[liquididx[i]]->coord[1];
    }
    comxy[0]/=nliqptcls;
    comxy[1]/=nliqptcls;
    Rvec hdensity=Rvec(numz, 0.);
    for(int i=0;i<nliqptcls;i++) {
        int index=static_cast<int>((particles[liquididx[i]]->coord[2]-surfaceb-0.5*dz)/dz);
        if(index>=0 && index<numz){
            particles_atz[index]->addParticle(particles[liquididx[i]]);
        real rdist=pow(pow(particles[liquididx[i]]->coord[0]-comxy[0], 2.)+pow(particles[liquididx[i]]->coord[1]-comxy[1], 2.), 0.5);
        if(rdist<hradius)
            hdensity[index]+=1.;
        }
    }

    for(int i=0;i<numz;i++){
        int numptcls=particles_atz[i]->size();
        radius[i]=0.;
        for(int j=0;j<numr;j++) rdensity[j]=0.;
        
        if(numptcls!=0){
            com[i]=particles_atz[i]->calculateCenterOfMass(ref0[i], ref1[i]);
            for(int j=0;j<numptcls;j++){
                Real3D dvec=pbc.getMinimumImageVector(particles_atz[i]->at(j)->coord, com[i]);
                real dist=pow(dvec[0]*dvec[0]+dvec[1]*dvec[1], 0.5);
                rdensity[int(dist/dr)]+=1.;
            }
            for(int j=0;j<numr;j++)
                rdensity[j]/=(j+0.5)*dr*dr*2*PI*dz;

            for(int j=numr-1;j>=1;j--){
                if(rdensity[j]<dcrit && rdensity[j-1]>=dcrit){
                    radius[i]=dr*(j-(rdensity[j]-dcrit)/(rdensity[j]-rdensity[j-1]));
                    break;
                }
            }
        }

        hdensity[i]/=hradius*hradius*dz*M_PI;
    }
    tevol[0][step]=*std::max_element(radius.begin(), radius.end());


    real height;
    for(int i=numz-1;i>=1;i--){
        height=0.;
        if(hdensity[i]<dcrit && hdensity[i-1]>=dcrit){
            height=dz*(i-(hdensity[i]-dcrit)/(hdensity[i]-hdensity[i-1]));
            break;
        }
    }
    for(int i=0;i<numz-1;i++){
        if(hdensity[i]<dcrit && hdensity[i+1]>=dcrit){
            height-=dz*(i-(hdensity[i]-dcrit)/(hdensity[i+1]-hdensity[i]));
            break;
        }
    }
    if(with_wall){
        tevol[1][step]=radius[0];
        tevol[2][step]=radius[1];
        tevol[4][step]=height;
        if(radius[0]==0)
            tevol[3][step]==180.0;
        else
            tevol[3][step]=atan2(dz, radius[0]-radius[1])*180/PI;

    }
    else{
        tevol[1][step]=height;
    }

    return;

}
