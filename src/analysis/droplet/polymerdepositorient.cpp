#include "polymerdepositorient.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
PolymerDepositOrient::PolymerDepositOrient(InitialSet initset):Property(initset){
    title="   Calculation of polymer orientation w.r.t. radial vector";

    nprops=1;
    for_timeevol=true;
    need_position=true;

    outfile="polymer_depori.out";


    outtitle_tevol="#order parameter of 2D nematics";
    outheader_tevol=Svec{"S"};

}

void PolymerDepositOrient::getSpecificParameters(){
    dz=1.0;
    surfaceb=control->getWallMinPosition()[2];
    command->getCommandSingleOption("-dz", dz, &dz);
    return;
}




void PolymerDepositOrient::initializeVariables(){
    Ivec polidx;

    for(int i=0;i<particles.size();i++){
        if(particles[i]->getMoleculeName().compare("POL")==0){
            int idx=particles[i]->getMoleculeIndex();
            if(std::find(polidx.begin(), polidx.end(), idx)==polidx.end()){
                polidx.push_back(idx);
            }
        }
    }
   
    num_pol=polidx.size();
    for(int i=0;i<num_pol;i++){
        polymers.push_back(new ParticleGroup(pbc));
    }

    for(int i=0;i<particles.size();i++){
        if(particles[i]->getMoleculeName().compare("POL")==0){
            int index=std::distance(polidx.begin(), std::find(polidx.begin(), polidx.end(), particles[i]->getMoleculeIndex()));
            polymers[index]->addParticle(particles[i]);
        }
    }


    dbin=0.1;
    dbegin=0.0;
    ndbin=int((box.sqr()/2)/dbin)+1;
    numvarperstep=num_pol;

    if(liquidgrps.size()==0){
        for(int i=0;i<nbeads;i++){
            liquid.addParticle(particles[i]);
        }
    }
    else{
        for(int i=0;i<nbeads;i++){
            if(std::find(liquidgrps.begin(), liquidgrps.end(), particles[i]->getMoleculeName())!=liquidgrps.end())
                liquid.addParticle(particles[i]);
        }
    }

    nliqptcls=liquididx.size();

    ref0=Real3D(box/2-1.0);
    ref1=Real3D(box/2+1.0);
    ref0[2]=surfaceb+dz*0.5;
    ref1[2]=ref0[2]+dz;

    initializeResultArrays();
    return;
}

void PolymerDepositOrient::calculateStep(int step){
    Real3D com=liquid.calculateCenterOfMass(ref0, ref1);
    for(int i=0;i<num_pol;i++){
        int lpol=polymers[i]->size();
        int hlpol=lpol/2;
        Real3D posdiff=pbc.getMinimumImageVector(polymers[i]->at(hlpol)->coord, com);
        Real3D remol(0.0);
        for(int j=0;j<lpol-1;j++){
            Real3D dr=pbc.getMinimumImageVector(polymers[i]->at(j+1)->coord, polymers[i]->at(j)->coord);
            remol+=dr;
        }
        Real3D v1=posdiff.unit();
        Real3D v2=remol.unit();
        real s=pow((v1[0]*v2[0]+v1[1]*v2[1])/sqrt(v1[0]*v1[0]+v1[1]*v1[1])/sqrt(v2[0]*v2[0]+v2[1]*v2[1]), 2.0)-1;
        tevol[0][step]+=s;
    }
    tevol[0][step]/=num_pol;
    
    return ;
}














