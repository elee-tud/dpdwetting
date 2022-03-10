#include "polymerbondlength.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
PolymerBondlength::PolymerBondlength(InitialSet initset):Property(initset){
    title="   Calculation of polymer bond length";

    nprops=1;
    for_timeevol=true;
    for_distrib=true;
    need_position=true;
    maxbl=10*control->getCellCutoff();

    outfile="bondlength.out";


    outtitle_tevol="#Time evolution of bond length";
    outheader_tevol=Svec{"Time", "<l>"};
    outtitle_dist="#Distribution of bond lengths";
    outheader_dist=Svec{"R^2", "P(l)"};

}


PolymerBondlength::~PolymerBondlength(){
    for(int i=0;i<num_pol;i++){
        delete polymers[i];
    }
}

void PolymerBondlength::getSpecificParameters(){
    molname="POL";
    command->getCommandSingleOption("-mol", molname, &molname);
    std::cout << "Caculation will be done for molecules with name=" << molname << std::endl;

    return;
}




void PolymerBondlength::initializeVariables(){
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

    dbin=0.01;
    dbegin=0.0;
    ndbin=int(maxbl/dbin)+1;
    numvarperstep=num_pol;

            
    initializeResultArrays();
    return;
}

void PolymerBondlength::calculateStep(int step){
    tevol[0][step]=0.;
    int numbond=0;
    for(int i=0;i<num_pol;i++){
        for(int j=0;j<polymers[i]->size()-1;j++){
            real bl=pbc.getMinimumImageVector(polymers[i]->at(j+1)->coord, polymers[i]->at(j)->coord).abs();
            tevol[0][step]+=bl;
            dist[0][int(bl/dbin)]+=1.;
        }
        numbond+=polymers[i]->size()-1;

    }
    tevol[0][step]/=numbond;
    return;
}







