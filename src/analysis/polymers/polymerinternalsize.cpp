#include "polymerinternalsize.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
PolymerInternalSize::PolymerInternalSize(InitialSet initset):Property(initset){
    title="   Calculation of polymer internal size ";

    nprops=1;
    for_distrib=true;
    need_position=true;

    outfile="polymer_subsize.out";


    outtitle_dist="#Distribution of End-to-end distance for subchain";
    outheader_dist=Svec{"Re^2(s)"};


}




void PolymerInternalSize::getSpecificParameters(){
    molname="POL";
    command->getCommandSingleOption("-mol", molname, &molname);
    std::cout << "Caculation will be done for molecules with name=" << molname << std::endl;

    return;
}


void PolymerInternalSize::initializeVariables(){
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
    dbin=1.;
    dbegin=-0.5;
    ndbin=polymers[0]->size();
    numvarperstep=num_pol;

            
    initializeResultArrays();
    return;
}

void PolymerInternalSize::calculateStep(int step){
    
    for(int i=0;i<num_pol;i++){
        Real3D remol(0.0);
        for(int j=0;j<polymers[i]->size()-1;j++){
            remol+=pbc.getMinimumImageVector(polymers[i]->at(j+1)->coord, polymers[i]->at(j)->coord);
            dist[0][j+1]+=remol.sqr();
        }
    }
    return;
}







