#include "polymerorientation.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
PolymerOrientation::PolymerOrientation(InitialSet initset):Property(initset){
    title="   Calculation of polymer orientation (Second Legendre Polynom.)";

    nprops=3;
    for_timeevol=true;
    need_position=true;

    outfile="polymer_orient.out";


    outtitle_tevol="#Second Legendre Polynomial of polymer vecotrs";
    outheader_tevol=Svec{"L2(costheta)_x", "L2(costheta)_y""L2(costheta)_z"};

}



void PolymerOrientation::getSpecificParameters(){
    molname="POL";
    command->getCommandSingleOption("-mol", molname, &molname);
    std::cout << "Caculation will be done for molecules with name=" << molname << std::endl;

    return;
}



void PolymerOrientation::initializeVariables(){
    Ivec polidx;

    for(int i=0;i<particles.size();i++){
        if(particles[i]->getMoleculeName().compare(molname)==0){
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
        if(particles[i]->getMoleculeName().compare(molname)==0){
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

void PolymerOrientation::calculateStep(int step){
    int nbonds=0;
    for(int i=0;i<num_pol;i++){
        int lpol=polymers[i]->size();
        for(int j=0;j<lpol-1;j++){
            Real3D pvec=pbc.getMinimumImageVector(polymers[i]->at(j+1)->coord, polymers[i]->at(j)->coord);
            real vnorm=pvec.abs();
            tevol[0][step]+=pow(pvec*Real3D{1,0,0}/vnorm, 2.);
            tevol[1][step]+=pow(pvec*Real3D{0,1,0}/vnorm, 2.);
            tevol[2][step]+=pow(pvec*Real3D{0,0,1}/vnorm, 2.);
        }
        nbonds+=lpol-1;
    }
    tevol[0][step]=1.5*tevol[0][step]/nbonds-0.5;
    tevol[1][step]=1.5*tevol[1][step]/nbonds-0.5;
    tevol[2][step]=1.5*tevol[2][step]/nbonds-0.5;





    return;
}







