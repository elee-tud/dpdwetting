#include "bridgedensxz.hpp"
#include "../../filecontrol.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
BridgeDensityXZ::BridgeDensityXZ(InitialSet initset):Property(initset){
    title="   Calculation of the density profile as a function of x,z";
    nprops=3;
    outheader_dist=Svec{"x", "z", "rho", "rho(S)", "rho(P)"};

    for_distrib=true;
    need_position=true;

    outfile="bridgedensxz.out";


    outtitle_dist="#Particle as a function of x,z-position";


}



void BridgeDensityXZ::getSpecificParameters(){
    dx=0.05;
    dz=0.05;
    command->getCommandSingleOption("-dx", dx, &dx);
    command->getCommandSingleOption("-dz", dz, &dz);
    
    surfaceb=control->getWallMinPosition()[2];
    surfacet=control->getWallMaxPosition()[2];
   
    liquidgrps=control->getLiquidGroups();
    return;
}


void BridgeDensityXZ::initializeVariables(){
    dbegin=0.;
    nx=static_cast<int>(box[0]/dz)+1;
    nz=static_cast<int>((surfacet-surfaceb)/dz)+1;
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

void BridgeDensityXZ::calculateStep(int step){
    for(int i=0;i<nliqptcls;i++){
        int xidx=static_cast<int>((particles[liquididx[i]]->coord[0])/dx);
        int zidx=static_cast<int>((particles[liquididx[i]]->coord[2]-surfaceb)/dz);
        int idx=xidx*nz+zidx;
        if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0){
            nsptclsx[idx]++;
            dist[0][idx]+=1;
            dist[1][idx]+=1;
        }
        else if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0){
            npptclsx[idx]++;
            dist[0][idx]+=1;
            dist[2][idx]+=1;
        }
    }


    return;

}

void BridgeDensityXZ::normalizeResults(){
    for(int j=0;j<nprops;j++){
        for(int i=0;i<ndbin;i++){
            dist_sum[j][i]/=dx*dz*nsteps*box[1];
        }
    }
    return;
}

void BridgeDensityXZ::writeOutput(){
    if(mpi->isMaster()){
        std::ofstream outstream;
        openFileWithBackup(outdistname, &outstream, false, false);
        outstream << outtitle_dist << std::endl;
        outstream << "#";
        outstream << std::setw(15) << "x";
        for(int i=1;i<5;i++){
            outstream << std::setw(16) << outheader_dist[i];
        }
        outstream << std::endl;

        for(int i=0;i<nx;i++){
            for(int j=0;j<nz;j++){
                int idx=i*nz+j;
                outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << i*dx;
                outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << j*dz;
                for(int k=0;k<3;k++){
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
    


