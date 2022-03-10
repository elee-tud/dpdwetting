#include "bridgesize.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
BridgeSize::BridgeSize(InitialSet initset):Property(initset){
    title="   Calculation of a contact angle of a bridge." ;
    nprops=2;
    outheader_tevol=Svec{"Time", "CA_rec", "CA_adv"};

    for_timeevol=true;
    need_position=true;

    outfile="bridgesize.out";


    outtitle_tevol="#Advancing and Receding contact angle of a bridge";


}



void BridgeSize::getSpecificParameters(){
    dz=0.65;
    dcrit=6.10/4;
    dr=0.5;
    pilheight=0.;
    command->getCommandSingleOption("-dz", dz, &dz);
    command->getCommandSingleOption("-dd", dcrit, &dcrit);
    command->getCommandSingleOption("-dr", dr, &dr);
    /*Considering the pillar height*/
    command->getCommandSingleOption("-ph", pilheight, &pilheight);
    surfaceb=control->getWallMinPosition()[2]+pilheight;
    surfacet=control->getWallMaxPosition()[2]-pilheight;
   
    liquidgrps=control->getLiquidGroups();
    return;
}


void BridgeSize::initializeVariables(){
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
    

    numx=static_cast<int>((box[0])/dr);
    numz=static_cast<int>((zcenter-surfaceb-dz*0.5)/dz);
    topdensity=Rvec2D(numz, Rvec(numx, 0));
    botdensity=Rvec2D(numz, Rvec(numx, 0));
    interface=Rvec2D(numz, Rvec(4, 0));
    initializeResultArrays();
    return;
}

void BridgeSize::calculateStep(int step){
    for(int i=0;i<numz;i++){
        interface[i][0]=0.;
        interface[i][1]=0.;
        interface[i][2]=0.;
        interface[i][3]=0.;
        for(int k=0;k<numx;k++){
            topdensity[i][k]=0.;
            botdensity[i][k]=0.;
        }
    }
    for(int i=0;i<nliqptcls;i++){
        int idxx=static_cast<int>(particles[liquididx[i]]->coord[0]/dr);
        real distzfc=particles[liquididx[i]]->coord[2]-zcenter;
        if(distzfc<0){ 
            int idxz=static_cast<int>((particles[liquididx[i]]->coord[2]-surfaceb-0.5*dz)/dz);
            if(idxz>=0 && idxz<numz)
                botdensity[idxz][idxx]+=1.;
        }
        else{
            int idxz=static_cast<int>((surfacet-particles[liquididx[i]]->coord[2]-0.5*dz)/dz);
            if(idxz>=0 && idxz<numz)
                topdensity[idxz][idxx]+=1.;
        }
    }
    for(int i=0;i<numz;i++){
        for(int k=0;k<numx;k++){
            topdensity[i][k]/=dr*box[1]*dz;
            botdensity[i][k]/=dr*box[1]*dz;
        }
    }
    /*
    for(int i=0;i<numx;i++){
        std::cout << i*dr << "   " << botdensity[0][i] << "   " << topdensity[0][i] << "   " << botdensity[1][i] << "   " << topdensity[1][i] << std::endl;
    }
    */
    for(int i=0;i<numz;i++){
        for(int k=0;k<numx-1;k++){
            if(botdensity[i][k]<dcrit && botdensity[i][k+1]>=dcrit){
                interface[i][0]=k*dr+dr/(botdensity[i][k+1]-botdensity[i][k])*(dcrit-botdensity[i][k]);
                break;
            }
        }
        for(int k=0;k<numx-1;k++){
            if(topdensity[i][k]<dcrit && topdensity[i][k+1]>=dcrit){
                interface[i][2]=k*dr+dr/(topdensity[i][k+1]-topdensity[i][k])*(dcrit-topdensity[i][k]);
                break;
            }
        }
        for(int k=numx-1;k>0;k--){
            if(botdensity[i][k]<dcrit && botdensity[i][k-1]>=dcrit){
                interface[i][1]=(k-1)*dr+dr/(botdensity[i][k]-botdensity[i][k-1])*(dcrit-botdensity[i][k-1]);
                break;
            }
        }
        for(int k=numx-1;k>0;k--){
            if(topdensity[i][k]<dcrit && topdensity[i][k-1]>=dcrit){
                interface[i][3]=(k-1)*dr+dr/(topdensity[i][k]-topdensity[i][k-1])*(dcrit-topdensity[i][k-1]);
                break;
            }
        }
    }
    tevol[0][step]=(atan2(dz, interface[1][0]-interface[0][0])+atan2(dz, interface[0][3]-interface[1][3]))*90/PI;
    tevol[1][step]=(atan2(dz, interface[0][1]-interface[1][1])+atan2(dz, interface[1][2]-interface[0][2]))*90/PI;




    return;

}
