#include "jumpfreq.hpp" 
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
JumpingFrequency::JumpingFrequency(InitialSet initset):Property(initset){
    title=  "Calculation of Cummulated Jumping Frequency";

    nprops=6;
    for_distrib=true;
    need_position=true;

    outfile="jumpfreq.out";


    outtitle_dist="#Cummulated jumping frequency";
    outheader_dist=Svec{"Time", "P1(t)", "P_1S(t)", "P_1P(t)", "P2(t)", "P_2S(t)", "P_2P(t)"};


}





void JumpingFrequency::getSpecificParameters(){
    refdt=200;
    dsz=0.91;
    dl=0.87;
    command->getCommandSingleOption("-rdt", refdt, &refdt);
    command->getCommandSingleOption("-dz", dsz, &dsz);
    command->getCommandSingleOption("-dl", dl, &dl);
    dlsqr=dl*dl;
    dszsqr=dsz*dsz;
    surfaceb=control->getWallMinPosition()[2];
    surfacet=control->getWallMaxPosition()[2];
    return;
}



void JumpingFrequency::initializeVariables(){

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
    dbegin=begstep;
    dbin=control->getTimeStep()*control->getTrajFrequency();
    ndbin=refdt;
    nliqptcls=liquididx.size();
    ref=R3vec(nliqptcls, Real3D(0.));
    totnsref=0;
    totnpref=0;

    initializeResultArrays();

    return;
}

void JumpingFrequency::calculateStep(int step){
    if(step%refdt==0){
        refstep=step;
        for(int i=0;i<nliqptcls;i++){
            ref[i]=particles[i]->coord;
        }
        srefidx.clear();
        prefidx.clear();
        for(int i=0;i<nliqptcls;i++){
            if(static_cast<int>((particles[liquididx[i]]->coord[2]-surfaceb)/dsz)==0){
                if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0 ){
                    srefidx.push_back(liquididx[i]);
                }
                else if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0 ){
                    prefidx.push_back(liquididx[i]);
                }
            }
        }
        nsref=srefidx.size();
        npref=prefidx.size();
        scummul=Ivec2D(nsref, Ivec(refdt, 0));
        pcummul=Ivec2D(npref, Ivec(refdt, 0));
    }

    else{
        int dt=step-refstep;
        for(int i=0;i<nsref;i++){
            if(scummul[i][dt-1]==0){
                Real3D vec=particles[srefidx[i]]->coord-ref[srefidx[i]];
                real xydist=vec[0]*vec[0]+vec[1]*vec[1];
                real zdist=vec[2]*vec[2];
                if(xydist>=dlsqr && zdist<dszsqr){
                    scummul[i][dt]=1;
                }
                else if(xydist<dlsqr && zdist>=dszsqr){
                    scummul[i][dt]=2;
                }
            }
            else if(scummul[i][dt-1]==1)
                scummul[i][dt]=1;
            else if(scummul[i][dt-1]==2)
                scummul[i][dt]=2;
        }

        for(int i=0;i<npref;i++){
            if(pcummul[i][dt-1]==0){
                Real3D vec=particles[prefidx[i]]->coord-ref[prefidx[i]];
                real xydist=vec[0]*vec[0]+vec[1]*vec[1];
                real zdist=vec[2]*vec[2];
                if(xydist>=dlsqr && zdist<dszsqr){
                    pcummul[i][dt]=1;
                }
                else if(xydist<dlsqr && zdist>=dszsqr){
                    pcummul[i][dt]=2;
                }
            }
            else if(pcummul[i][dt-1]==1)
                pcummul[i][dt]=1;
            else if(pcummul[i][dt-1]==2)
                pcummul[i][dt]=2;
        }

        if(step%refdt==refdt-1){
            for(int i=0;i<nsref;i++){
                for(int j=0;j<refdt;j++){
                    if(scummul[i][j]==1){
                        dist[0][j]++;
                        dist[1][j]++;
                    }
                    else if(scummul[i][j]==2){
                        dist[3][j]++;
                        dist[4][j]++;
                    }
                }
            }
            for(int i=0;i<npref;i++){
                for(int j=0;j<refdt;j++){
                    if(pcummul[i][j]==1){
                        dist[0][j]++;
                        dist[2][j]++;
                    }
                    else if(pcummul[i][j]==2){
                        dist[3][j]++;
                        dist[5][j]++;
                    }
                }
            }
            totnsref+=nsref;
            totnpref+=npref;
        }

    }
    return;

}

void JumpingFrequency::normalizeResults(){
    for(int i=0;i<refdt;i++){
        dist_sum[0][i]/=(real)(totnsref+totnpref);
        dist_sum[1][i]/=(real)totnsref;
        dist_sum[2][i]/=(real)totnpref;
        dist_sum[3][i]/=(real)(totnsref+totnpref);
        dist_sum[4][i]/=(real)totnsref;
        dist_sum[5][i]/=(real)totnpref;
    }
    return;
}




