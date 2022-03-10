#include "bridgevelx.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
BridgeVelocityX::BridgeVelocityX(InitialSet initset):Property(initset){
    title="   Calculation of the velocity profile as a function of x";
    nprops=9;
    outheader_dist=Svec{"x", "v_x(x)", "v_x,S(x)", "v_x,P(x)", "v_y(x)", "v_y,S(x)", "v_y,P(x)", "v_z(x)", "v_z,S(x)", "v_z,P(x)"};

    for_distrib=true;
    need_position=true;
    need_velocity=true;

    outfile="bridgevelx.out";


    outtitle_dist="#Velocity components as a function of x-position at surface";


}



void BridgeVelocityX::getSpecificParameters(){
    dbin=1.0;
    dz=0.91;
    command->getCommandSingleOption("-dx", dbin, &dbin);
    command->getCommandSingleOption("-dz", dz, &dz);
    
    surfaceb=control->getWallMinPosition()[2];
    surfacet=control->getWallMaxPosition()[2];
   
    liquidgrps=control->getLiquidGroups();
    return;
}


void BridgeVelocityX::initializeVariables(){
    dbegin=0.;
    ndbin=static_cast<int>(box[0]/dbin)+1;
    
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
    nsptclsx=Ivec(ndbin, 0);
    npptclsx=Ivec(ndbin, 0);


    return;
}

void BridgeVelocityX::calculateStep(int step){
    for(int i=0;i<nliqptcls;i++){
        if(static_cast<int>((particles[liquididx[i]]->coord[2]-surfaceb)/dz)==0){
            int zidx=static_cast<int>((particles[liquididx[i]]->coord[0])/dbin);
            if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0){
                nsptclsx[zidx]++;
                dist[0][zidx]+=particles[liquididx[i]]->veloc[0];
                dist[1][zidx]+=particles[liquididx[i]]->veloc[0];
                dist[3][zidx]+=particles[liquididx[i]]->veloc[1];
                dist[4][zidx]+=particles[liquididx[i]]->veloc[1];
                dist[6][zidx]+=particles[liquididx[i]]->veloc[2];
                dist[7][zidx]+=particles[liquididx[i]]->veloc[2];
            }
            else if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0){
                npptclsx[zidx]++;
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

void BridgeVelocityX::normalizeResults(){
    for(int i=0;i<ndbin;i++){
        if(nsptclsx[i]>nsteps){
            dist_sum[1][i]/=nsptclsx[i];
            dist_sum[4][i]/=nsptclsx[i];
            dist_sum[7][i]/=nsptclsx[i];
        }
        else if(npptclsx[i]>nsteps){
            dist_sum[2][i]/=npptclsx[i];
            dist_sum[5][i]/=npptclsx[i];
            dist_sum[8][i]/=npptclsx[i];
        }
        else{
            dist_sum[1][i]=0.;
            dist_sum[4][i]=0.;
            dist_sum[7][i]=0.;
            dist_sum[2][i]=0.;
            dist_sum[5][i]=0.;
            dist_sum[8][i]=0.;
        }
        if(nsptclsx[i]+npptclsx[i]>nsteps){
            dist_sum[0][i]/=nsptclsx[i]+npptclsx[i];
            dist_sum[3][i]/=nsptclsx[i]+npptclsx[i];
            dist_sum[6][i]/=nsptclsx[i]+npptclsx[i];
        }
        else{
            dist_sum[0][i]=0.;
            dist_sum[3][i]=0.;
            dist_sum[6][i]=0.;
        }
    }
    return;
}
