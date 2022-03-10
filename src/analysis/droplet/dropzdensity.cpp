#include "dropzdensity.hpp"
#include <algorithm>


using namespace dpd;

DropZDensity::DropZDensity(InitialSet initset):Property(initset){
    title= "  Calculation of z-density";
    
    nprops=3;
    for_distrib=true;
    need_position=true;

    outfile="zdensity.out";

    outtitle_dist="#Axial solvent, polymer and all particle densities from COM";
    outheader_dist=Svec{"r", "rho(r)", "rho_S(r)", "rho_P(r)"};

}


void DropZDensity::getSpecificParameters(){
    dbin=0.02;
    drcenter=3.0;
    command->getCommandSingleOption("-dr", dbin, &dbin);
    command->getCommandSingleOption("-dc", drcenter, &drcenter);
    surfaceb=control->getWallMinPosition()[2];
    surfacet=control->getWallMaxPosition()[2];
    liquidgrps=control->getLiquidGroups();
    return;
}


void DropZDensity::initializeVariables(){
    center=Rvec{box[0]/2, box[1]/2};
    drsqr=drcenter*drcenter;
    ndbin=static_cast<int>((surfacet-surfaceb)/dbin)+1;
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
    return;

}

void DropZDensity::calculateStep(int step){
    for(int i=0;i<nliqptcls;i++){
        real distance=pow(particles[liquididx[i]]->coord[0]-center[0], 2)+pow(particles[liquididx[i]]->coord[1]-center[1], 2);
        if(distance<=drsqr){
            int index=static_cast<int>((particles[liquididx[i]]->coord[2]-surfaceb)/dbin);
            if(index<ndbin){
                if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0)
                    dist[1][index]+=1;
                else if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0)
                    dist[2][index]+=1;
                dist[0][index]+=1.;
            }
        }
    }


    return;
}

void DropZDensity::normalizeResults(){
    for(int j=0;j<ndbin;j++){
        dist_sum[0][j]/=nsteps*drsqr*PI*dbin;
        dist_sum[1][j]/=nsteps*drsqr*PI*dbin;
        dist_sum[2][j]/=nsteps*drsqr*PI*dbin;
    }
    return;
}




