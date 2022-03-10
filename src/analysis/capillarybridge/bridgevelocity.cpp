#include "bridgevelocity.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
BridgeVelocity::BridgeVelocity(InitialSet initset):Property(initset){
    title="   Calculation velocity profile along z";
    nprops=3;

    for_distrib=true;
    need_velocity=true;

    outfile="slipvelocity.out";


    outtitle_dist="#Velocity profile along z axis";
    outheader_dist=Svec{"Time", "v(z)", "v_sol(z)", "v_pol(z)"};


}



void BridgeVelocity::getSpecificParameters(){
    ndbin=20;
    command->getCommandSingleOption("-nb", ndbin, &ndbin);
    surfaceb=control->getWallMinPosition()[2];
    surfacet=control->getWallMaxPosition()[2];
   
    return;
}


void BridgeVelocity::initializeVariables(){
    dbin=(surfacet-surfaceb)/ndbin;
    npolatz=Ivec(ndbin, 0); 
    nsolatz=Ivec(ndbin, 0); 

    initializeResultArrays();
    return;
}

void BridgeVelocity::calculateStep(int step){
    for(int i=0;i<particles.size();i++){
        int idx=static_cast<int>((particles[i]->coord[2]-surfaceb)/dbin);
        if(particles[i]->getMoleculeName().compare("POL")==0){
            npolatz[idx]++;
            dist[2][idx]+=particles[i]->veloc[0];
            dist[0][idx]+=particles[i]->veloc[0];
        }
        else if(particles[i]->getMoleculeName().compare("SOL")==0){
            nsolatz[idx]++;
            dist[1][idx]+=particles[i]->veloc[0];
            dist[0][idx]+=particles[i]->veloc[0];
        }

    }
    return;
}
void BridgeVelocity::normalizeResults(){
    for(int i=0;i<ndbin;i++){
        dist_sum[0][i]/=(npolatz[i]+nsolatz[i]);
        dist_sum[1][i]/=nsolatz[i];
        dist_sum[2][i]/=npolatz[i];
    }
    return;
}


