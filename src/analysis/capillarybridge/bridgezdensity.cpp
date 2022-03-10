#include "bridgezdensity.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
BridgeZdensity::BridgeZdensity(InitialSet initset):Property(initset){
    title="   Density Profile along z";
    nprops=3;
    outheader_dist=Svec{"Time", "rho(x)", "rho_S(x)", "rho_P(x)"};

    for_distrib=true;
    need_position=true;

    outfile="zdensity.out";


    outtitle_dist="#Density profile along z";


}



void BridgeZdensity::getSpecificParameters(){
    dbin=0.02;
    dx=10.0;
    command->getCommandSingleOption("-dz", dbin, &dbin);
    command->getCommandSingleOption("-dx", dx, &dx);
    surfaceb=control->getWallMinPosition()[2];
    surfacet=control->getWallMaxPosition()[2];
   
    liquidgrps=control->getLiquidGroups();
    return;
}


void BridgeZdensity::initializeVariables(){
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
    
    xmin=box[0]/2-dx;
    xmax=box[0]/2+dx;
    ndbin=static_cast<int>((surfacet-surfaceb)/dbin)+1;
    initializeResultArrays();
    return;
}

void BridgeZdensity::calculateStep(int step){
    for(int i=0;i<nliqptcls;i++){
        int idxz=static_cast<int>((particles[liquididx[i]]->coord[2]-surfaceb)/dbin);
        if(particles[liquididx[i]]->coord[0]>=xmin && particles[liquididx[i]]->coord[0]<xmax){
            if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0){
                dist[2][idxz]+=1.;
                dist[0][idxz]+=1.;
            }
            else if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0){
                dist[1][idxz]+=1.;
                dist[0][idxz]+=1.;
            }
        }
    }
    return;
}
void BridgeZdensity::normalizeResults(){
    if(mpi->isMaster()){
        for(int i=0;i<nprops;i++){
            for(int j=0;j<ndbin;j++){
                dist_sum[i][j]/=nsteps*dbin*box[1]*dx*2;
            }
        }
    }
    return;
}
