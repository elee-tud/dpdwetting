#include "polymerstretch.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
PolymerStretch::PolymerStretch(InitialSet initset):Property(initset){
    title="   Calculation of polymer stretching factor";

    nprops=3;
    for_timeevol=true;
    need_position=true;

    outfile="polymer_stretch.out";


    outtitle_tevol="#Polymer stretching factor";
    outheader_tevol=Svec{"Time", "<S_x>", "<S_y>", "<S_z>"};

}





void PolymerStretch::initializeVariables(){
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

    for(int i=0;i<num_pol;i++){
        ref.push_back(polymers[i]->at(polymers[i]->size()/2));
    }

    dbin=0.1;
    dbegin=0.0;
    ndbin=int((box.sqr()/2)/dbin)+1;
    numvarperstep=num_pol;

            
    initializeResultArrays();
    return;
}

void PolymerStretch::calculateStep(int step){
    for(int i=0;i<num_pol;i++){
        int nmonpol=polymers[i]->size();
        R3vec newcoord(nmonpol, Real3D(0.));
        newcoord[0]=polymers[i]->at(0)->coord;
        for(int j=0;j<nmonpol-1;j++){
            Real3D pbcvec=pbc.getMinimumImageVector(polymers[i]->at(j+1)->coord, polymers[i]->at(j)->coord);
            newcoord[j+1]=newcoord[j]+pbcvec;
        }

        Real3D max=newcoord[0];
        Real3D min=newcoord[0];
        for(int j=1;j<nmonpol;j++){
            if(newcoord[j][0]>max[0])
                max[0]=newcoord[j][0];
            if(newcoord[j][1]>max[1])
                max[1]=newcoord[j][1];
            if(newcoord[j][2]>max[2])
                max[2]=newcoord[j][2];
            if(newcoord[j][0]<min[0])
                min[0]=newcoord[j][0];
            if(newcoord[j][1]<min[1])
                min[1]=newcoord[j][1];
            if(newcoord[j][2]<min[2])
                min[2]=newcoord[j][2];
        }
        tevol[0][step]+=max[0]-min[0];
        tevol[1][step]+=max[1]-min[1];
        tevol[2][step]+=max[2]-min[2];

    }
    tevol[0][step]/=num_pol;
    tevol[1][step]/=num_pol;
    tevol[2][step]/=num_pol;
    return;
}







