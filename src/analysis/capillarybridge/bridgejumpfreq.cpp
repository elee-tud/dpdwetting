#include "bridgejumpfreq.hpp" 
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
BridgeJumpingFrequency::BridgeJumpingFrequency(InitialSet initset):Property(initset){
    title=  "Calculation of Cummulated Jumping Frequency";

    nprops=12;
    for_distrib=true;
    need_position=true;

    outfile="jumpfreq.out";


    outtitle_dist="#Cummulated jumping frequency";
    outheader_dist=Svec{"Time", "P1_adv(t)", "P_1S,adv(t)", "P_1,adv(t)", "P2_adv(t)", "P_2S,adv(t)", "P_2P,adv(t)", "P1_rec(t)", "P_1S,rec(t)", "P_1,rec(t)", "P2_rec(t)", "P_2S,rec(t)", "P_2P,rec(t)"};


}





void BridgeJumpingFrequency::getSpecificParameters(){
    refdt=200;
    dsz=1.0;
    dl=0.87;
    dx=0.5;
    dcrit=1.5;
    dcl=3.0;
    command->getCommandSingleOption("-rdt", refdt, &refdt);
    command->getCommandSingleOption("-dsz", dsz, &dsz);
    command->getCommandSingleOption("-dl", dl, &dl);
    command->getCommandSingleOption("-dx", dx, &dx);
    command->getCommandSingleOption("-dc", dcrit, &dcrit);
    command->getCommandSingleOption("-dcl", dcl, &dcl);
    dlsqr=dl*dl;
    surfaceb=control->getWallMinPosition()[2];
    surfacet=control->getWallMaxPosition()[2];
    shearrate=control->getWallShearRate();
    return;
}



void BridgeJumpingFrequency::initializeVariables(){

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
    totnsref_adv=0;
    totnsref_rec=0;
    totnpref_adv=0;
    totnpref_rec=0;
    wdisp=shearrate*(surfacet-surfaceb)/2*dbin;

    
    numx=static_cast<int>((box[0])/dx);
    topdensity=Rvec(numx, 0);
    botdensity=Rvec(numx, 0);
    initializeResultArrays();

    return;
}

void BridgeJumpingFrequency::calculateInterfacePosition(){
    for(int k=0;k<numx;k++){
        topdensity[k]=0.;
        botdensity[k]=0.;
    }
    for(int i=0;i<nliqptcls;i++){
        int idxx=static_cast<int>(particles[liquididx[i]]->coord[0]/dx);
        if(particles[liquididx[i]]->coord[2]>=surfaceb && particles[liquididx[i]]->coord[2]<surfaceb+dsz){
            botdensity[idxx]+=1.;
        }
        if(particles[liquididx[i]]->coord[2]>=surfacet-dsz && particles[liquididx[i]]->coord[2]<surfacet){
            topdensity[idxx]+=1.;
        }
    }
    
    for(int k=0;k<numx;k++){
        topdensity[k]/=dx*box[1]*dsz;
        botdensity[k]/=dx*box[1]*dsz;
    }

    for(int k=0;k<numx-1;k++){
        if(botdensity[k]<dcrit && botdensity[k+1]>=dcrit){
            recinter_bot=k*dx+dx/(botdensity[k+1]-botdensity[k])*(dcrit-botdensity[k]);
            break;
        }
    }
    for(int k=0;k<numx-1;k++){
        if(topdensity[k]<dcrit && topdensity[k+1]>=dcrit){
            advinter_top=k*dx+dx/(topdensity[k+1]-topdensity[k])*(dcrit-topdensity[k]);
            break;
        }
    }
    for(int k=numx-1;k>0;k--){
        if(botdensity[k]<dcrit && botdensity[k-1]>=dcrit){
            advinter_bot=(k-1)*dx+dx/(botdensity[k]-botdensity[k-1])*(dcrit-botdensity[k-1]);
            break;
        }
    }
    for(int k=numx-1;k>0;k--){
        if(topdensity[k]<dcrit && topdensity[k-1]>=dcrit){
            recinter_top=(k-1)*dx+dx/(topdensity[k]-topdensity[k-1])*(dcrit-topdensity[k-1]);
            break;
        }
    }
//    std::cout << advinter_top << "," << advinter_bot << "," << recinter_top << "," << recinter_bot << std::endl;
    advinter_top+=dcl;
    advinter_bot-=dcl;
    recinter_top-=dcl;
    recinter_bot+=dcl;

    return;
}


void BridgeJumpingFrequency::calculateStep(int step){
    calculateInterfacePosition();
    if(step%refdt==0){
        refstep=step;
        for(int i=0;i<nliqptcls;i++){
            ref[i]=particles[i]->coord;
        }
        srefidx_rectop.clear();
        prefidx_rectop.clear();
        srefidx_recbot.clear();
        prefidx_recbot.clear();
        srefidx_advtop.clear();
        prefidx_advtop.clear();
        srefidx_advbot.clear();
        prefidx_advbot.clear();
        for(int i=0;i<nliqptcls;i++){
            if(static_cast<int>((particles[liquididx[i]]->coord[2]-surfaceb)/dsz)==0){
                if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0 ){
                    if(particles[liquididx[i]]->coord[0]<recinter_bot)
                        srefidx_recbot.push_back(liquididx[i]);
                    else if(particles[liquididx[i]]->coord[0]>=advinter_bot)
                        srefidx_advbot.push_back(liquididx[i]);
                }
                else if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0 ){
                    if(particles[liquididx[i]]->coord[0]<recinter_bot)
                        prefidx_recbot.push_back(liquididx[i]);
                    else if(particles[liquididx[i]]->coord[0]>=advinter_bot)
                        prefidx_advbot.push_back(liquididx[i]);
                }
            }
            else if(static_cast<int>((surfacet-particles[liquididx[i]]->coord[2])/dsz)==0){
                if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0 ){
                    if(particles[liquididx[i]]->coord[0]<advinter_top)
                        srefidx_advtop.push_back(liquididx[i]);
                    else if(particles[liquididx[i]]->coord[0]>=recinter_top)
                        srefidx_rectop.push_back(liquididx[i]);
                }
                else if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0 ){
                    if(particles[liquididx[i]]->coord[0]<advinter_top)
                        prefidx_advtop.push_back(liquididx[i]);
                    else if(particles[liquididx[i]]->coord[0]>=recinter_top)
                        prefidx_rectop.push_back(liquididx[i]);
                }
            }
        }
        nsref_rectop=srefidx_rectop.size();
        npref_rectop=prefidx_rectop.size();
        nsref_advtop=srefidx_advtop.size();
        npref_advtop=prefidx_advtop.size();
        nsref_recbot=srefidx_recbot.size();
        npref_recbot=prefidx_recbot.size();
        nsref_advbot=srefidx_advbot.size();
        npref_advbot=prefidx_advbot.size();
        scummul_rectop=Ivec2D(nsref_rectop, Ivec(refdt, 0));
        pcummul_rectop=Ivec2D(npref_rectop, Ivec(refdt, 0));
        scummul_advtop=Ivec2D(nsref_advtop, Ivec(refdt, 0));
        pcummul_advtop=Ivec2D(npref_advtop, Ivec(refdt, 0));
        scummul_recbot=Ivec2D(nsref_recbot, Ivec(refdt, 0));
        pcummul_recbot=Ivec2D(npref_recbot, Ivec(refdt, 0));
        scummul_advbot=Ivec2D(nsref_advbot, Ivec(refdt, 0));
        pcummul_advbot=Ivec2D(npref_advbot, Ivec(refdt, 0));
    }

    else{
        int dt=step-refstep;
        /*Solvent for receding line(bottom)*/
        for(int i=0;i<nsref_recbot;i++){
            if(scummul_recbot[i][dt-1]==STATIC){
                if(static_cast<int>((particles[srefidx_recbot[i]]->coord[2]-surfaceb)/dsz)==0){
                    if(particles[srefidx_recbot[i]]->coord[0]>=recinter_bot)
                        scummul_recbot[i][dt]=OUTINTER;
                    else{
                        Real3D vec=particles[srefidx_recbot[i]]->coord-ref[srefidx_recbot[i]];
                        real distance=(vec[0]+wdisp)*(vec[0]+wdisp)+vec[1]*vec[1];
                        if(distance > dlsqr){
                            scummul_recbot[i][dt]=PARDIFF;
                        }
                    }
                }
                else{
                    scummul_recbot[i][dt]=PERDIFF;
                }
            }
            else if(scummul_recbot[i][dt-1]==PARDIFF)
                scummul_recbot[i][dt]=PARDIFF;
            else if(scummul_recbot[i][dt-1]==PERDIFF)
                scummul_recbot[i][dt]=PERDIFF;
            else if(scummul_recbot[i][dt-1]==OUTINTER)
                scummul_recbot[i][dt]=OUTINTER;
        }
       
        /*Solvent for receding line(top)*/
        for(int i=0;i<nsref_rectop;i++){
            if(scummul_rectop[i][dt-1]==STATIC){
                if(static_cast<int>((surfacet-particles[srefidx_rectop[i]]->coord[2])/dsz)==0){
                    if(particles[srefidx_rectop[i]]->coord[0]<recinter_top)
                        scummul_rectop[i][dt]=OUTINTER;
                    else{
                        Real3D vec=particles[srefidx_rectop[i]]->coord-ref[srefidx_rectop[i]];
                        real distance=(vec[0]-wdisp)*(vec[0]-wdisp)+vec[1]*vec[1];
                        if(distance > dlsqr){
                            scummul_rectop[i][dt]=PARDIFF;
                        }
                    }
                }
                else{
                    scummul_rectop[i][dt]=PERDIFF;
                }
            }
            else if(scummul_rectop[i][dt-1]==PARDIFF)
                scummul_rectop[i][dt]=PARDIFF;
            else if(scummul_rectop[i][dt-1]==PERDIFF)
                scummul_rectop[i][dt]=PERDIFF;
            else if(scummul_rectop[i][dt-1]==OUTINTER)
                scummul_rectop[i][dt]=OUTINTER;
        }

        /*Solvent for advancing line(bottom)*/
        for(int i=0;i<nsref_advbot;i++){
            if(scummul_advbot[i][dt-1]==STATIC){
                if(static_cast<int>((particles[srefidx_advbot[i]]->coord[2]-surfaceb)/dsz)==0){
                    if(particles[srefidx_advbot[i]]->coord[0]<advinter_bot)
                        scummul_advbot[i][dt]=OUTINTER;
                    else{
                        Real3D vec=particles[srefidx_advbot[i]]->coord-ref[srefidx_advbot[i]];
                        real distance=(vec[0]+wdisp)*(vec[0]+wdisp)+vec[1]*vec[1];
                        if(distance > dlsqr){
                            scummul_advbot[i][dt]=PARDIFF;
                        }
                    }
                }
                else{
                    scummul_advbot[i][dt]=PERDIFF;
                }
            }
            else if(scummul_advbot[i][dt-1]==PARDIFF)
                scummul_advbot[i][dt]=PARDIFF;
            else if(scummul_advbot[i][dt-1]==PERDIFF)
                scummul_advbot[i][dt]=PERDIFF;
            else if(scummul_advbot[i][dt-1]==OUTINTER)
                scummul_advbot[i][dt]=OUTINTER;
        }
       
        /*Solvent for advancing line(top)*/
        for(int i=0;i<nsref_advtop;i++){
            if(scummul_advtop[i][dt-1]==STATIC){
                if(static_cast<int>((surfacet-particles[srefidx_advtop[i]]->coord[2])/dsz)==0){
                    if(particles[srefidx_advtop[i]]->coord[0]>=advinter_top)
                        scummul_advtop[i][dt]=OUTINTER;
                    else{
                        Real3D vec=particles[srefidx_advtop[i]]->coord-ref[srefidx_advtop[i]];
                        real distance=(vec[0]-wdisp)*(vec[0]-wdisp)+vec[1]*vec[1];
                        if(distance > dlsqr){
                            scummul_advtop[i][dt]=PARDIFF;
                        }
                    }
                }
                else{
                    scummul_advtop[i][dt]=PERDIFF;
                }
            }
            else if(scummul_advtop[i][dt-1]==PARDIFF)
                scummul_advtop[i][dt]=PARDIFF;
            else if(scummul_advtop[i][dt-1]==PERDIFF)
                scummul_advtop[i][dt]=PERDIFF;
            else if(scummul_advtop[i][dt-1]==OUTINTER)
                scummul_advtop[i][dt]=OUTINTER;
        }


        /*Polymer for receding line(bottom)*/
        for(int i=0;i<npref_recbot;i++){
            if(pcummul_recbot[i][dt-1]==STATIC){
                if(static_cast<int>((particles[prefidx_recbot[i]]->coord[2]-surfaceb)/dsz)==0){
                    if(particles[prefidx_recbot[i]]->coord[0]>=recinter_bot)
                        pcummul_recbot[i][dt]=OUTINTER;
                    else{
                        Real3D vec=particles[prefidx_recbot[i]]->coord-ref[prefidx_recbot[i]];
                        real distance=(vec[0]+wdisp)*(vec[0]+wdisp)+vec[1]*vec[1];
                        if(distance > dlsqr){
                            pcummul_recbot[i][dt]=PARDIFF;
                        }
                    }
                }
                else{
                    pcummul_recbot[i][dt]=PERDIFF;
                }
            }
            else if(pcummul_recbot[i][dt-1]==PARDIFF)
                pcummul_recbot[i][dt]=PARDIFF;
            else if(pcummul_recbot[i][dt-1]==PERDIFF)
                pcummul_recbot[i][dt]=PERDIFF;
            else if(pcummul_recbot[i][dt-1]==OUTINTER)
                pcummul_recbot[i][dt]=OUTINTER;
        }
       
        /*Polymer for receding line(top)*/
        for(int i=0;i<npref_rectop;i++){
            if(pcummul_rectop[i][dt-1]==STATIC){
                if(static_cast<int>((surfacet-particles[prefidx_rectop[i]]->coord[2])/dsz)==0){
                    if(particles[prefidx_rectop[i]]->coord[0]<recinter_top)
                        pcummul_rectop[i][dt]=OUTINTER;
                    else{
                        Real3D vec=particles[prefidx_rectop[i]]->coord-ref[prefidx_rectop[i]];
                        real distance=(vec[0]-wdisp)*(vec[0]-wdisp)+vec[1]*vec[1];
                        if(distance > dlsqr){
                            pcummul_rectop[i][dt]=PARDIFF;
                        }
                    }
                }
                else{
                    pcummul_rectop[i][dt]=PERDIFF;
                }
            }
            else if(pcummul_rectop[i][dt-1]==PARDIFF)
                pcummul_rectop[i][dt]=PARDIFF;
            else if(pcummul_rectop[i][dt-1]==PERDIFF)
                pcummul_rectop[i][dt]=PERDIFF;
            else if(pcummul_rectop[i][dt-1]==OUTINTER)
                pcummul_rectop[i][dt]=OUTINTER;
        }

        /*Polymer for advancing line(bottom)*/
        for(int i=0;i<npref_advbot;i++){
            if(pcummul_advbot[i][dt-1]==STATIC){
                if(static_cast<int>((particles[prefidx_advbot[i]]->coord[2]-surfaceb)/dsz)==0){
                    if(particles[prefidx_advbot[i]]->coord[0]<advinter_bot)
                        pcummul_advbot[i][dt]=OUTINTER;
                    else{
                        Real3D vec=particles[prefidx_advbot[i]]->coord-ref[prefidx_advbot[i]];
                        real distance=(vec[0]+wdisp)*(vec[0]-wdisp)+vec[1]*vec[1];
                        if(distance > dlsqr){
                            pcummul_advbot[i][dt]=PARDIFF;
                        }
                    }
                }
                else{
                    pcummul_advbot[i][dt]=PERDIFF;
                }
            }
            else if(pcummul_advbot[i][dt-1]==PARDIFF)
                pcummul_advbot[i][dt]=PARDIFF;
            else if(pcummul_advbot[i][dt-1]==PERDIFF)
                pcummul_advbot[i][dt]=PERDIFF;
            else if(pcummul_advbot[i][dt-1]==OUTINTER)
                pcummul_advbot[i][dt]=OUTINTER;
        }
       
        /*Polymer for advancing line(top)*/
        for(int i=0;i<npref_advtop;i++){
            if(pcummul_advtop[i][dt-1]==STATIC){
                if(static_cast<int>((surfacet-particles[prefidx_advtop[i]]->coord[2])/dsz)==0){
                    if(particles[prefidx_advtop[i]]->coord[0]>=advinter_top)
                        scummul_advtop[i][dt]=OUTINTER;
                    else{
                        Real3D vec=particles[prefidx_advtop[i]]->coord-ref[prefidx_advtop[i]];
                        real distance=(vec[0]+wdisp)*(vec[0]+wdisp)+vec[1]*vec[1];
                        if(distance > dlsqr){
                            pcummul_advtop[i][dt]=PARDIFF;
                        }
                    }
                }
                else{
                    pcummul_advtop[i][dt]=PERDIFF;
                }
            }
            else if(pcummul_advtop[i][dt-1]==PARDIFF)
                pcummul_advtop[i][dt]=PARDIFF;
            else if(pcummul_advtop[i][dt-1]==PERDIFF)
                pcummul_advtop[i][dt]=PERDIFF;
            else if(pcummul_advtop[i][dt-1]==OUTINTER)
                pcummul_advtop[i][dt]=OUTINTER;
        }



        /*Calculating cummulative fraction*/
        if(step%refdt==refdt-1){
            for(int i=0;i<nsref_advbot;i++){
                for(int j=0;j<refdt;j++){
                    if(scummul_advbot[i][j]==PARDIFF){
                        dist[0][j]++;
                        dist[1][j]++;
                    }
                    else if(scummul_advbot[i][j]==PERDIFF){
                        dist[3][j]++;
                        dist[4][j]++;
                    }
                }
            }

            for(int i=0;i<nsref_advtop;i++){
                for(int j=0;j<refdt;j++){
                    if(scummul_advtop[i][j]==PARDIFF){
                        dist[0][j]++;
                        dist[1][j]++;
                    }
                    else if(scummul_advtop[i][j]==PERDIFF){
                        dist[3][j]++;
                        dist[4][j]++;
                    }
                }
            }

            for(int i=0;i<npref_advbot;i++){
                for(int j=0;j<refdt;j++){
                    if(pcummul_advbot[i][j]==PARDIFF){
                        dist[0][j]++;
                        dist[2][j]++;
                    }
                    else if(pcummul_advbot[i][j]==PERDIFF){
                        dist[3][j]++;
                        dist[5][j]++;
                    }
                }
            }
            for(int i=0;i<npref_advtop;i++){
                for(int j=0;j<refdt;j++){
                    if(pcummul_advtop[i][j]==PARDIFF){
                        dist[0][j]++;
                        dist[2][j]++;
                    }
                    else if(pcummul_advtop[i][j]==PERDIFF){
                        dist[3][j]++;
                        dist[5][j]++;
                    }
                }
            }

            for(int i=0;i<nsref_recbot;i++){
                for(int j=0;j<refdt;j++){
                    if(scummul_recbot[i][j]==PARDIFF){
                        dist[6][j]++;
                        dist[7][j]++;
                    }
                    else if(scummul_recbot[i][j]==PERDIFF){
                        dist[9][j]++;
                        dist[10][j]++;
                    }
                }
            }

            for(int i=0;i<nsref_rectop;i++){
                for(int j=0;j<refdt;j++){
                    if(scummul_rectop[i][j]==PARDIFF){
                        dist[6][j]++;
                        dist[7][j]++;
                    }
                    else if(scummul_rectop[i][j]==PERDIFF){
                        dist[9][j]++;
                        dist[10][j]++;
                    }
                }
            }

            for(int i=0;i<npref_recbot;i++){
                for(int j=0;j<refdt;j++){
                    if(pcummul_recbot[i][j]==PARDIFF){
                        dist[6][j]++;
                        dist[8][j]++;
                    }
                    else if(pcummul_recbot[i][j]==PERDIFF){
                        dist[9][j]++;
                        dist[11][j]++;
                    }
                }
            }
            for(int i=0;i<npref_rectop;i++){
                for(int j=0;j<refdt;j++){
                    if(pcummul_rectop[i][j]==PARDIFF){
                        dist[6][j]++;
                        dist[8][j]++;
                    }
                    else if(pcummul_rectop[i][j]==PERDIFF){
                        dist[9][j]++;
                        dist[11][j]++;
                    }
                }
            }
            totnsref_adv+=nsref_advbot+nsref_advtop;
            totnsref_rec+=nsref_recbot+nsref_rectop;
            totnpref_adv+=npref_advbot+npref_advtop;
            totnpref_rec+=npref_recbot+npref_rectop;
        }

    }
    return;

}



void BridgeJumpingFrequency::normalizeResults(){
    for(int i=0;i<refdt;i++){
        dist_sum[0][i]/=(real)(totnsref_adv+totnpref_adv);
        dist_sum[1][i]/=(real)totnsref_adv;
        dist_sum[2][i]/=(real)totnpref_adv;
        dist_sum[3][i]/=(real)(totnsref_adv+totnpref_adv);
        dist_sum[4][i]/=(real)totnsref_adv;
        dist_sum[5][i]/=(real)totnpref_adv;
        dist_sum[6][i]/=(real)(totnsref_rec+totnpref_rec);
        dist_sum[7][i]/=(real)totnsref_rec;
        dist_sum[8][i]/=(real)totnpref_rec;
        dist_sum[9][i]/=(real)(totnsref_rec+totnpref_rec);
        dist_sum[10][i]/=(real)totnsref_rec;
        dist_sum[11][i]/=(real)totnpref_rec;
    }
    return;
}




