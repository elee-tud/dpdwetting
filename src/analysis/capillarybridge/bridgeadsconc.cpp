#include "bridgeadsconc.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
BridgeAdsorptionConc::BridgeAdsorptionConc(InitialSet initset):Property(initset){
    title="   Calculation of the polymer concentration on the adsorption layer";
    nprops=2;

    for_distrib=true;
    need_position=true;

    outfile="bridgeadsconc.out";


    outtitle_dist="#Polymer concentration as a function of x";
    outheader_dist=Svec{"Time", "rho_p(x)(bot)", "rho_p(x)(top)"};


}



void BridgeAdsorptionConc::getSpecificParameters(){
    dr=1.0;
    command->getCommandSingleOption("-dr", dr, &dr);
    surfaceb=control->getWallMinPosition()[2];
    surfacet=control->getWallMaxPosition()[2];
   
    liquidgrps=control->getLiquidGroups();
    return;
}


void BridgeAdsorptionConc::initializeVariables(){
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
    
    dbegin=0;
    dbin=dr;
    ndbin=static_cast<int>((box[0])/dbin);
    npolbot=Ivec(ndbin, 0);
    nsolbot=Ivec(ndbin, 0);
    npoltop=Ivec(ndbin, 0);
    nsoltop=Ivec(ndbin, 0);
    initializeResultArrays();
    return;
}

void BridgeAdsorptionConc::calculateStep(int step){
    for(int i=0;i<ndbin;i++){
        npolbot[i]=0;
        nsolbot[i]=0;
        npoltop[i]=0;
        nsoltop[i]=0;
    }

    for(int i=0;i<nliqptcls;i++){
        int idxx=static_cast<int>(particles[liquididx[i]]->coord[0]/dr);
        if(particles[liquididx[i]]->coord[2]>=surfaceb && particles[liquididx[i]]->coord[2]<surfaceb+dr){
            if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0)
                npolbot[idxx]++;
            else if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0)
                nsolbot[idxx]++;
        }
        else if(particles[liquididx[i]]->coord[2]<=surfacet && particles[liquididx[i]]->coord[2]>surfacet-dr){
            if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0)
                npoltop[idxx]++;
            else if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0)
                nsoltop[idxx]++;
        }
    }
    for(int i=0;i<ndbin;i++){
        if(npolbot[i]!=0 || nsolbot[i]!=0)
            dist[0][i]+=static_cast<real>(npolbot[i])/(npolbot[i]+nsolbot[i]);
        if(npoltop[i]!=0 || nsoltop[i]!=0)
            dist[1][i]+=static_cast<real>(npoltop[i])/(npoltop[i]+nsoltop[i]);
    }

}


void BridgeAdsorptionConc::normalizeResults(){
    for(int i=0;i<ndbin;i++){
        dist_sum[0][i]=dist[0][i]/static_cast<real>(nsteps);
        dist_sum[1][i]=dist[1][i]/static_cast<real>(nsteps);
    }



    return;

}
