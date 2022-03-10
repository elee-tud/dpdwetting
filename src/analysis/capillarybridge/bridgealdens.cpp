#include "bridgealdens.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
BridgeAdsorblayerDensity::BridgeAdsorblayerDensity(InitialSet initset):Property(initset){
    title="   Density Profile at the first adsorption layer";
    nprops=6;
    outheader_dist=Svec{"Time", "rhob(x)", "rhob_S(x)", "rhob_P(x)", "rhot(x)", "rhot_S(x)", "rhot_P(x)"};

    for_distrib=true;
    need_position=true;

    outfile="adslayer_density.out";


    outtitle_dist="#Density profile along x in the first adsorption layer";


}



void BridgeAdsorblayerDensity::getSpecificParameters(){
    dz=1.45;
    dbin=0.02;
    command->getCommandSingleOption("-dz", dz, &dz);
    command->getCommandSingleOption("-dx", dbin, &dbin);
    surfaceb=control->getWallMinPosition()[2];
    surfacet=control->getWallMaxPosition()[2];
   
    liquidgrps=control->getLiquidGroups();
    return;
}


void BridgeAdsorblayerDensity::initializeVariables(){
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
    zcenter=box[2]/2;
    

    ndbin=static_cast<int>((box[0])/dbin)+1;
    initializeResultArrays();
    return;
}

void BridgeAdsorblayerDensity::calculateStep(int step){
    for(int i=0;i<nliqptcls;i++){
        int idxx=static_cast<int>(particles[liquididx[i]]->coord[0]/dbin);
        if(particles[liquididx[i]]->coord[2]>=surfaceb && particles[liquididx[i]]->coord[2]<surfaceb+dz){
            if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0){
                dist[2][idxx]+=1.;
                dist[0][idxx]+=1.;
            }
            else if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0){
                dist[1][idxx]+=1.;
                dist[0][idxx]+=1.;
            }
        }
        else if(particles[liquididx[i]]->coord[2]<surfacet && particles[liquididx[i]]->coord[2]>=surfacet-dz){
            if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0){
                dist[5][idxx]+=1.;
                dist[3][idxx]+=1.;
            }
            else if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0){
                dist[4][idxx]+=1.;
                dist[3][idxx]+=1.;
            }
        }
    }
    return;
}
void BridgeAdsorblayerDensity::normalizeResults(){
    if(mpi->isMaster()){
        for(int i=0;i<nprops;i++){
            for(int j=0;j<ndbin;j++){
                dist_sum[i][j]/=nsteps*dbin*dz*box[1];
            }
        }
    }
    return;
}
