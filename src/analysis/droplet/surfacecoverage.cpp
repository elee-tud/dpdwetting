#include "surfacecoverage.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
SurfaceCoverage::SurfaceCoverage(InitialSet initset):Property(initset){
    title="   Calculation of SurfaceCoverage";

    nprops=4;
    for_timeevol=true;
    need_position=true;

    outfile="surface_coverage.out";


    outtitle_tevol="#Surface Coverage";
    outheader_tevol=Svec{"Time", "coverage", "Diameter", "ads_Solv", "ads_Pol"};

}


void SurfaceCoverage::getSpecificParameters(){
    acrit=1.0;

    dz=1.0;
    dcrit=3.04;
    dr=0.5;

    surfaceb=control->getWallMinPosition()[2];
    surfacet=control->getWallMaxPosition()[2];
    command->getCommandSingleOption("-ac", acrit, &acrit);
    command->getCommandSingleOption("-dz", dz, &dz);
    command->getCommandSingleOption("-dd", dcrit, &dcrit);
    command->getCommandSingleOption("-dr", dr, &dr);
    liquidgrps=control->getLiquidGroups();
    return;
}



void SurfaceCoverage::initializeVariables(){
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

    
//    numz=static_cast<int>((surfacet-surfaceb)/dz)+1;
    particles_bottom=ParticleGroup(pbc);
    com=Real3D(0.0);
    ref0=Real3D(box/2-1.0);
    ref1=Real3D(box/2+1.0);
    ref0[2]=surfaceb+dz*0.5;
    ref1[2]=ref0[2]+dz;
    
    numr=static_cast<int>(pow(box[1]*box[1]+box[0]*box[0], 0.5)*1.4142135/dr)+1;
    rdensity=Rvec(numr, 0.);
    initializeResultArrays();
    return;
}


void SurfaceCoverage::calculateStep(int step){
    particles_bottom.clear();
    for(int i=0;i<nliqptcls;i++){
        int index=static_cast<int>((particles[liquididx[i]]->coord[2]-surfaceb-0.5*dz)/dz);
        if(index==0)
            particles_bottom.addParticle(particles[liquididx[i]]);
    }
    int numptcls=particles_bottom.size();
    real radius=0.;
    for(int i=0;i<numr;i++) rdensity[i]=0.;

    if(numptcls>0){
        Real3D com=particles_bottom.calculateCenterOfMass(ref0, ref1);
        for(int i=0;i<numptcls;i++){
            Real3D dvec=pbc.getMinimumImageVector(particles_bottom.at(i)->coord, com);
            real dist=pow(dvec[0]*dvec[0]+dvec[1]*dvec[1], 0.5);
            rdensity[int(dist/dr)]+=1.;
        }
        for(int i=0;i<numr;i++)
            rdensity[i]/=(i+0.5)*dr*dr*2*PI*dz;

        for(int i=numr-1;i>=1;i--){
            if(rdensity[i]<dcrit && rdensity[i-1]>=dcrit){
                radius=dr*(i-(rdensity[i]-dcrit)/(rdensity[i]-rdensity[i-1]));
                break;
            }
        }

        for(int i=0;i<numptcls;i++){
            Real3D dvec=pbc.getMinimumImageVector(particles_bottom.at(i)->coord, com);
            real dist=pow(dvec[0]*dvec[0]+dvec[1]*dvec[1], 0.5);
            if(dist<=radius){
                if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0)
                    tevol[2][step]+=1.;
                else if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0)
                    tevol[3][step]+=1.;
            }
        }
    }
    tevol[1][step]=radius*2;
    tevol[0][step]=tevol[3][step]/(tevol[2][step]+tevol[3][step]);

    return;
}







