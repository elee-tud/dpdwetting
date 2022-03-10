#include "msd.hpp" 
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
MeanSquareDisplacement::MeanSquareDisplacement(InitialSet initset):Property(initset){
    title=  "Calculation of Mean Square Displacement";

    nprops=3;
    for_dynamics=true;
    need_position=true;

    outfile="polymer_msd.out";


    outtitle_dynm="#Mean Square Displacement";
    outheader_dynm=Svec{"Time", "g_1(t)", "g_2(t)", "g_3(t)"};


}





void MeanSquareDisplacement::getSpecificParameters(){
    molname="POL";
    command->getCommandSingleOption("-mol", molname, &molname);
    std::cout << "Caculation will be done for molecules with name=" << molname << std::endl;

    return;
}



void MeanSquareDisplacement::initializeVariables(){

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
    mon_pos=R3vec2D(nsteps, R3vec(num_pol, Real3D(0.)));
    com_pos=R3vec2D(nsteps, R3vec(num_pol, Real3D(0.)));
    mon_com=R3vec2D(nsteps, R3vec(num_pol, Real3D(0.)));

            
    initializeResultArrays();
    return;
}

void MeanSquareDisplacement::calculateStep(int step){
    for(int i=0;i<num_pol;i++){
        int lpol=polymers[i]->size();
        mon_pos[step][i]=ref[i]->coord;
        com_pos[step][i]=polymers[i]->calculateCenterOfMass();
        mon_com[step][i]=pbc.getMinimumImageVector(mon_pos[step][i], com_pos[step][i]);
    }
    return;

}

void MeanSquareDisplacement::calculateDynamicProperty(){
    unfoldTrajectory();
    for(int i=0;i<ndsteps;i++){
        for(int j=0;j<numavgs[i];j++){
            for(int k=0;k<num_pol;k++){
                tevol_sum[0][i]+=(mon_pos[j+dsteps[i]][k]-mon_pos[j][k]).sqr();
                tevol_sum[1][i]+=(mon_com[j+dsteps[i]][k]-mon_com[j][k]).sqr();
                tevol_sum[2][i]+=(com_pos[j+dsteps[i]][k]-com_pos[j][k]).sqr();
            }
        }
    }
    for(int i=0;i<ndsteps;i++){
        tevol_sum[0][i]/=(float)num_pol*numavgs[i];
        tevol_sum[1][i]/=(float)num_pol*numavgs[i];
        tevol_sum[2][i]/=(float)num_pol*numavgs[i];
    }
    return;
}


void MeanSquareDisplacement::unfoldTrajectory(){
    std::cout << "Final box=" << box << std::endl;
    Ivec3D bcross(nsteps, Ivec2D(num_pol, Ivec(3)));
    for(int i=1;i<nsteps;i++){
        for(int j=0;j<num_pol;j++){
            Ivec bcr=pbc.isCrossingBox(mon_pos[i][j], mon_pos[i-1][j]);
            bcross[i][j][0]=bcross[i-1][j][0]+bcr[0];
            bcross[i][j][1]=bcross[i-1][j][1]+bcr[1];
            bcross[i][j][2]=bcross[i-1][j][2]+bcr[2];
        }
    }
    for(int i=1;i<nsteps;i++){
        for(int j=0;j<num_pol;j++){
            mon_pos[i][j][0]+=bcross[i][j][0]*box[0];
            mon_pos[i][j][1]+=bcross[i][j][1]*box[1];
            mon_pos[i][j][2]+=bcross[i][j][2]*box[2];
        }
    }
    bcross=Ivec3D(nsteps, Ivec2D(num_pol, Ivec(3)));
    for(int i=1;i<nsteps;i++){
        for(int j=0;j<num_pol;j++){
            Ivec bcr=pbc.isCrossingBox(com_pos[i][j], com_pos[i-1][j]);
            bcross[i][j][0]=bcross[i-1][j][0]+bcr[0];
            bcross[i][j][1]=bcross[i-1][j][1]+bcr[1];
            bcross[i][j][2]=bcross[i-1][j][2]+bcr[2];
        }
    }
    for(int i=1;i<nsteps;i++){
        for(int j=0;j<num_pol;j++){
            com_pos[i][j][0]+=bcross[i][j][0]*box[0];
            com_pos[i][j][1]+=bcross[i][j][1]*box[1];
            com_pos[i][j][2]+=bcross[i][j][2]*box[2];
        }
    }
    return;
}







