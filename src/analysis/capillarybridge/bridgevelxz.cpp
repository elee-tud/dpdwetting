#include "bridgevelxz.hpp"
#include "../../filecontrol.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
BridgeVelocityXZ::BridgeVelocityXZ(InitialSet initset):Property(initset){
    title="   Calculation of the velocity profile as a function of x,z";
    nprops=9;
    outheader_dist=Svec{"x", "z", "v_x", "v_x(S)", "v_x(P)", "v_y", "v_y(S)", "v_y(P)", "v_z", "v_z(S)", "v_z(P)"};

    for_distrib=true;
    need_position=true;
    need_velocity=true;

    outfile="bridgevelxz.out";


    outtitle_dist="#Velocity components as a function of x,z-position";


}



void BridgeVelocityXZ::getSpecificParameters(){
    dx=1.0;
    dz=1.0;
    command->getCommandSingleOption("-dx", dx, &dx);
    command->getCommandSingleOption("-dz", dz, &dz);
    
    surfaceb=control->getWallMinPosition()[2];
    surfacet=control->getWallMaxPosition()[2];
   
    liquidgrps=control->getLiquidGroups();
    return;
}


void BridgeVelocityXZ::initializeVariables(){
    dbegin=0.;
    nx=static_cast<int>(box[0]/dz)+1;
    nz=static_cast<int>(box[2]/dz)+1;
    ndbin=nx*nz;
    
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

void BridgeVelocityXZ::calculateStep(int step){
    for(int i=0;i<nliqptcls;i++){
        int xidx=static_cast<int>((particles[liquididx[i]]->coord[0])/dx);
        int zidx=static_cast<int>((particles[liquididx[i]]->coord[2])/dz);
        int idx=xidx*nz+zidx;
        if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0){
            nsptclsx[idx]++;
            dist[0][idx]+=particles[liquididx[i]]->veloc[0];
            dist[1][idx]+=particles[liquididx[i]]->veloc[0];
            dist[3][idx]+=particles[liquididx[i]]->veloc[1];
            dist[4][idx]+=particles[liquididx[i]]->veloc[1];
            dist[6][idx]+=particles[liquididx[i]]->veloc[2];
            dist[7][idx]+=particles[liquididx[i]]->veloc[2];
        }
        else if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0){
            npptclsx[idx]++;
            dist[0][idx]+=particles[liquididx[i]]->veloc[0];
            dist[2][idx]+=particles[liquididx[i]]->veloc[0];
            dist[3][idx]+=particles[liquididx[i]]->veloc[1];
            dist[5][idx]+=particles[liquididx[i]]->veloc[1];
            dist[6][idx]+=particles[liquididx[i]]->veloc[2];
            dist[8][idx]+=particles[liquididx[i]]->veloc[2];
        }
    }


    return;

}

void BridgeVelocityXZ::normalizeResults(){
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

void BridgeVelocityXZ::writeOutput(){
    if(mpi->isMaster()){
        std::ofstream outstream;
        openFileWithBackup(outdistname, &outstream, false, false);
        outstream << outtitle_dist << std::endl;
        outstream << "#";
        outstream << std::setw(15) << "x";
        for(int i=1;i<11;i++){
            outstream << std::setw(16) << outheader_dist[i];
        }
        outstream << std::endl;

        for(int i=0;i<nx;i++){
            for(int j=0;j<nz;j++){
                int idx=i*nz+j;
                outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << i*dx;
                outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << j*dz;
                for(int k=0;k<9;k++){
                    outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << dist_sum[k][idx];
                }
                outstream << std::endl;
            }
            outstream << std::endl;
        }
        outstream.close();
        std::cout << "Results are written in " << outdistname << "." << std::endl;
    }
    return;
}
    


