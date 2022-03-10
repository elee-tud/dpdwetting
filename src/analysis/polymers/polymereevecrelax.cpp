#include "polymereevecrelax.hpp" 
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
PolymerEevecRelax::PolymerEevecRelax(InitialSet initset):Property(initset){
    title=  "Calculation of polymer end-to-end vector autocorelation function.";

    nprops=1;
    for_dynamics=true;
    need_position=true;

    outfile="polymer_evacf.out";


    outtitle_dynm="#Auto-correlation function of polymer end-to-end vectors";
    outheader_dynm=Svec{"Time", "ACF(t)"};


}



void PolymerEevecRelax::getSpecificParameters(){
    molname="POL";
    command->getCommandSingleOption("-mol", molname, &molname);
    std::cout << "Caculation will be done for molecules with name=" << molname << std::endl;
    n_mod=0;
    command->getCommandSingleOption("-mod", n_mod, &n_mod);
    return;
}




void PolymerEevecRelax::initializeVariables(){

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
    if(n_mod==0){
        n_mod=polymers[0]->size()-1;
        mon_first=0;
        mon_end=n_mod+mon_first-1;
    }
    else{
        mon_first=static_cast<int>((polymers[0]->size()-n_mod)/2);
        mon_end=n_mod+mon_first-1;
    }

    etoe=R3vec2D(nsteps, R3vec(num_pol, Real3D(0.)));
            
    initializeResultArrays();
    return;
}

void PolymerEevecRelax::calculateStep(int step){
    for(int i=0;i<num_pol;i++){
        int lpol=polymers[i]->size();
        for(int j=mon_first;j<mon_end;j++){
            etoe[step][i]+=pbc.getMinimumImageVector(polymers[i]->at(j+1)->coord, polymers[i]->at(j)->coord);
        }
    }

}

void PolymerEevecRelax::calculateDynamicProperty(){
    for(int i=0;i<ndsteps;i++){
        for(int j=0;j<numavgs[i];j++){
            for(int k=0;k<num_pol;k++){
                tevol_sum[0][i]+=etoe[j+dsteps[i]][k]*etoe[j][k];
            }
        }
    }
    for(int i=0;i<ndsteps;i++){
        tevol_sum[0][i]/=(float)num_pol*numavgs[i];
    }
    return;
}








