#include "bridgeslipvel.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
BridgeSlipVelocity::BridgeSlipVelocity(InitialSet initset):Property(initset){
    title="   Calculation of the velocity components as a function of position in z";
    nprops=9;
    outheader_dist=Svec{"z", "v_x(z)", "v_x,S(z)", "v_x,P(z)", "v_y(z)", "v_y,S(z), v_y,P(z)", "v_z(z)", "v_z,S(z)", "v_z,P(z)"};

    for_distrib=true;
    need_position=true;
    need_velocity=true;

    outfile="bridgeslipvel.out";


    outtitle_dist="#Velocity components as a function of z-position";


}



void BridgeSlipVelocity::getSpecificParameters(){
    dbin=1.0;
    dxfc=2.0;
    command->getCommandSingleOption("-dx", dxfc, &dxfc);
    command->getCommandSingleOption("-dz", dbin, &dbin);
    
    surfaceb=control->getWallMinPosition()[2];
    surfacet=control->getWallMaxPosition()[2];
    zsep=surfacet-surfaceb;
   
    liquidgrps=control->getLiquidGroups();
    return;
}


void BridgeSlipVelocity::initializeVariables(){
    dbegin=0.;
    ndbin=static_cast<int>(zsep/dbin)+1;
    
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

    initializeResultArrays();
    nsptclsz=Ivec(ndbin, 0);
    npptclsz=Ivec(ndbin, 0);
    nlptclsz=Ivec(ndbin, 0);
    xcenter=Rvec(ndbin, 0);
    xmax=Rvec(ndbin, 0);
    xmin=Rvec(ndbin, 0);
    


    return;
}

void BridgeSlipVelocity::calculateStep(int step){
    for(int i=0;i<ndbin;i++){
        xcenter[i]=0.;
        nlptclsz[i]=0;
    }
    for(int i=0;i<nliqptcls;i++){
        int zidx=static_cast<int>((particles[liquididx[i]]->coord[2]-surfaceb)/dbin);
        xcenter[zidx]+=particles[liquididx[i]]->coord[0];
        nlptclsz[zidx]++;
    }
    for(int i=0;i<ndbin;i++){
        xcenter[i]/=static_cast<real>(nlptclsz[i]);
        xmin[i]=xcenter[i]-dxfc;
        xmax[i]=xcenter[i]+dxfc;
    }




    for(int i=0;i<nliqptcls;i++){
        int zidx=static_cast<int>((particles[liquididx[i]]->coord[2]-surfaceb)/dbin);
        if(particles[liquididx[i]]->coord[0]>=xmin[zidx] && particles[liquididx[i]]->coord[0]<xmax[zidx] ){
            if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0){
                nsptclsz[zidx]++;
                dist[0][zidx]+=particles[liquididx[i]]->veloc[0];
                dist[1][zidx]+=particles[liquididx[i]]->veloc[0];
                dist[3][zidx]+=particles[liquididx[i]]->veloc[1];
                dist[4][zidx]+=particles[liquididx[i]]->veloc[1];
                dist[6][zidx]+=particles[liquididx[i]]->veloc[2];
                dist[7][zidx]+=particles[liquididx[i]]->veloc[2];
            }
            else if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0){
                npptclsz[zidx]++;
                dist[0][zidx]+=particles[liquididx[i]]->veloc[0];
                dist[2][zidx]+=particles[liquididx[i]]->veloc[0];
                dist[3][zidx]+=particles[liquididx[i]]->veloc[1];
                dist[5][zidx]+=particles[liquididx[i]]->veloc[1];
                dist[6][zidx]+=particles[liquididx[i]]->veloc[2];
                dist[8][zidx]+=particles[liquididx[i]]->veloc[2];
            }
        }

    }


    return;

}

void BridgeSlipVelocity::normalizeResults(){
    for(int i=0;i<ndbin;i++){
        dist_sum[0][i]/=nsptclsz[i]+npptclsz[i];
        dist_sum[1][i]/=nsptclsz[i];
        dist_sum[2][i]/=npptclsz[i];
        dist_sum[3][i]/=nsptclsz[i]+npptclsz[i];
        dist_sum[4][i]/=nsptclsz[i];
        dist_sum[5][i]/=npptclsz[i];
        dist_sum[6][i]/=nsptclsz[i]+npptclsz[i];
        dist_sum[7][i]/=nsptclsz[i];
        dist_sum[8][i]/=npptclsz[i];
    }
    dbegin=(surfaceb-surfacet)/2;
    return;
}
