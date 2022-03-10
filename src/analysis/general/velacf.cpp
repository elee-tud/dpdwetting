#include "velacf.hpp" 
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
VelocityAutoCorrelation::VelocityAutoCorrelation(InitialSet initset):Property(initset){
    title=  "Velocity auto-correlation function";

    nprops=3;
    for_dynamics=true;
    need_velocity=true;

    outfile="velacf.out";


    outtitle_dynm="#Velocity Auto-correlation function";
    outheader_dynm=Svec{"Time", "C_s(t)", "C_p(t)", "C(t)"};


}







void VelocityAutoCorrelation::initializeVariables(){

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
    pol_vel=R3vec2D(nsteps, R3vec(num_pol, Real3D(0.)));

    for(int i=0;i<particles.size();i++){
        if(particles[i]->getMoleculeName().compare("SOL")==0){
            solidx.push_back(i);
            if(solidx.size()>=MAXSOL)
                break;
        }
    }
    num_sol=solidx.size();
    sol_vel=R3vec2D(nsteps, R3vec(num_sol, Real3D(0.)));

            
    initializeResultArrays();
    return;
}

void VelocityAutoCorrelation::calculateStep(int step){
    for(int i=0;i<num_pol;i++){
        pol_vel[step][i]=ref[i]->veloc;
    }
    for(int i=0;i<num_sol;i++){
        sol_vel[step][i]=particles[solidx[i]]->veloc;
    }
    return;

}

void VelocityAutoCorrelation::calculateDynamicProperty(){
    for(int i=0;i<ndsteps;i++){
        for(int j=0;j<numavgs[i];j++){
            for(int k=0;k<num_pol;k++){
                tevol_sum[0][i]+=pol_vel[j+dsteps[i]][k]*pol_vel[j][k];
                tevol_sum[2][i]+=pol_vel[j+dsteps[i]][k]*pol_vel[j][k];
            }
            for(int k=0;k<num_sol;k++){
                tevol_sum[1][i]+=sol_vel[j+dsteps[i]][k]*sol_vel[j][k];
                tevol_sum[2][i]+=sol_vel[j+dsteps[i]][k]*sol_vel[j][k];
            }
        }
    }
    for(int i=0;i<ndsteps;i++){
        tevol_sum[0][i]/=(float)num_pol*numavgs[i];
        tevol_sum[1][i]/=(float)num_sol*numavgs[i];
        tevol_sum[2][i]/=(float)(num_pol+num_sol)*numavgs[i];
    }
    return;
}


