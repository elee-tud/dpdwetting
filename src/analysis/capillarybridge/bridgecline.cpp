#include "bridgecline.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
BridgeCline::BridgeCline(InitialSet initset):Property(initset){
    title="   Calculation of a contact line position";
    nprops=4;
    outheader_tevol=Svec{"Time", "x_cl1", "x_cl2", "x_cl3", "x_cl4"};

    for_timeevol=true;
    need_position=true;

    outfile="bridgecontline.out";


    outtitle_tevol="#Contact line position as a function of time (bottom left, top left, bottom right, bottom top)";


}



void BridgeCline::getSpecificParameters(){
    dz=0.8;
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


void BridgeCline::initializeVariables(){
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
    topdensity=Rvec(numx, 0);
    botdensity=Rvec(numx, 0);
    initializeResultArrays();
    return;
}

void BridgeCline::calculateStep(int step){
    for(int k=0;k<numx;k++){
        topdensity[k]=0.;
        botdensity[k]=0.;
    }
    for(int i=0;i<nliqptcls;i++){
        int idxx=static_cast<int>(particles[liquididx[i]]->coord[0]/dr);
        real distzfc=particles[liquididx[i]]->coord[2]-zcenter;
        if(distzfc<0){ 
            int idxz=static_cast<int>((particles[liquididx[i]]->coord[2]-surfaceb)/dz);
            if(idxz==0)
                botdensity[idxx]+=1.;
        }
        else{
            int idxz=static_cast<int>((surfacet-particles[liquididx[i]]->coord[2])/dz);
            if(idxz==0)
                topdensity[idxx]+=1.;
        }
    }
    for(int k=0;k<numx;k++){
        topdensity[k]/=dr*box[1]*dz;
        botdensity[k]/=dr*box[1]*dz;
    }
    /*
    for(int i=0;i<numx;i++){
        std::cout << i*dr << "   " << botdensity[0][i] << "   " << topdensity[0][i] << "   " << botdensity[1][i] << "   " << topdensity[1][i] << std::endl;
    }
    */
    int centerxidx=numx/2;
    for(int i=0;i<4;i++){
        tevol[i][step]=0.;
    }
    //Interface bottom left
    for(int k=centerxidx;k>0;k--){
        if(botdensity[k]>=dcrit && botdensity[k-1]<dcrit){
            tevol[0][step]=dr*(k-(botdensity[k]-dcrit)/(botdensity[k]-botdensity[k-1]));
            break;
        }
    }
    //Interface top left
    for(int k=centerxidx;k>0;k--){
        if(topdensity[k]>=dcrit && topdensity[k-1]<dcrit){
            tevol[1][step]=dr*(k-(topdensity[k]-dcrit)/(topdensity[k]-topdensity[k-1]));
            break;
        }
    }
    //Interface bottom right
    for(int k=centerxidx;k<numx-1;k++){
        if(botdensity[k]>=dcrit && botdensity[k+1]<dcrit){
            tevol[2][step]=dr*(k+(botdensity[k]-dcrit)/(botdensity[k]-botdensity[k+1]));
            break;
        }
    }
    //Interface top right
    for(int k=centerxidx;k<numx-1;k++){
        if(topdensity[k]>=dcrit && topdensity[k+1]<dcrit){
            tevol[3][step]=dr*(k+(topdensity[k]-dcrit)/(topdensity[k]-topdensity[k+1]));
            break;
        }
    }




    return;

}
