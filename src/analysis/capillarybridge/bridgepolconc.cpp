#include "bridgepolconc.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
BridgePolymerConc::BridgePolymerConc(InitialSet initset):Property(initset){
    title="   Calculation of the polymer concentration at the contact line zone";
    nprops=2;
    outheader_tevol=Svec{"Time", "x_rec", "x_adv"};

    for_timeevol=true;
    need_position=true;

    outfile="bridgepolconc.out";


    outtitle_tevol="#Polymer concentration at the contact line zone";


}



void BridgePolymerConc::getSpecificParameters(){
    dz=1.45;
    dcrit=1.3;
    dx=0.5;
    dxclz=2.00;
    command->getCommandSingleOption("-dz", dz, &dz);
    command->getCommandSingleOption("-dd", dcrit, &dcrit);
    command->getCommandSingleOption("-dr", dx, &dx);
    command->getCommandSingleOption("-dcl", dxclz, &dxclz);
    surfaceb=control->getWallMinPosition()[2];
    surfacet=control->getWallMaxPosition()[2];
   
    liquidgrps=control->getLiquidGroups();
    return;
}


void BridgePolymerConc::initializeVariables(){
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
    

    numx=static_cast<int>((box[0])/dx);
    topdensity=Rvec(numx, 0);
    botdensity=Rvec(numx, 0);
    interface=Rvec(4, 0);
    initializeResultArrays();
    return;
}

void BridgePolymerConc::calculateStep(int step){
    interface[0]=0.;
    interface[1]=0.;
    interface[2]=0.;
    interface[3]=0.;
    for(int k=0;k<numx;k++){
        topdensity[k]=0.;
        botdensity[k]=0.;
    }
    for(int i=0;i<nliqptcls;i++){
        int idxx=static_cast<int>(particles[liquididx[i]]->coord[0]/dx);
        if(particles[liquididx[i]]->coord[2]>=surfaceb && particles[liquididx[i]]->coord[2]<surfaceb+dz){
            botdensity[idxx]+=1.;
        }
        else if(particles[liquididx[i]]->coord[2]>=surfacet-dz && particles[liquididx[i]]->coord[2]<surfacet){
            topdensity[idxx]+=1.;
        }
    }
    for(int k=0;k<numx;k++){
        topdensity[k]/=dx*box[1]*dz;
        botdensity[k]/=dx*box[1]*dz;
    }
    /*
    for(int i=0;i<numx;i++){
        std::cout << i*dr << "   " << botdensity[0][i] << "   " << topdensity[0][i] << "   " << botdensity[1][i] << "   " << topdensity[1][i] << std::endl;
    }
    */
    for(int k=0;k<numx-1;k++){
        if(botdensity[k]<dcrit && botdensity[k+1]>=dcrit){
            interface[0]=k*dx+dx/(botdensity[k+1]-botdensity[k])*(dcrit-botdensity[k]);
            break;
        }
    }
    for(int k=0;k<numx-1;k++){
        if(topdensity[k]<dcrit && topdensity[k+1]>=dcrit){
            interface[2]=k*dx+dx/(topdensity[k+1]-topdensity[k])*(dcrit-topdensity[k]);
            break;
        }
    }
    for(int k=numx-1;k>0;k--){
        if(botdensity[k]<dcrit && botdensity[k-1]>=dcrit){
            interface[1]=(k-1)*dx+dx/(botdensity[k]-botdensity[k-1])*(dcrit-botdensity[k-1]);
            break;
        }
    }
    for(int k=numx-1;k>0;k--){
        if(topdensity[k]<dcrit && topdensity[k-1]>=dcrit){
            interface[3]=(k-1)*dx+dx/(topdensity[k]-topdensity[k-1])*(dcrit-topdensity[k-1]);
            break;
        }
    }

    int npolrec=0, npoladv=0;
    int nsolrec=0, nsoladv=0;

    for(int i=0;i<nliqptcls;i++){
        if(particles[liquididx[i]]->coord[2]>=surfaceb && particles[liquididx[i]]->coord[2]<surfaceb+dz){
            if(particles[liquididx[i]]->coord[0]>=interface[0]-dxclz && particles[liquididx[i]]->coord[0]<interface[0]+dxclz){
                if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0)
                    npolrec++;
                else
                    nsolrec++;
            }
            else if(particles[liquididx[i]]->coord[0]>=interface[1]-dxclz && particles[liquididx[i]]->coord[0]<interface[1]+dxclz){
                if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0)
                    npoladv++;
                else
                    nsoladv++;
            }
        }

        if(particles[liquididx[i]]->coord[2]>=surfacet-dz && particles[liquididx[i]]->coord[2]<=surfacet){
            if(particles[liquididx[i]]->coord[3]>=interface[0]-dxclz && particles[liquididx[i]]->coord[0]<interface[3]+dxclz){
                if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0)
                    npolrec++;
                else
                    nsolrec++;
            }
            else if(particles[liquididx[i]]->coord[0]>=interface[2]-dxclz && particles[liquididx[i]]->coord[0]<interface[2]+dxclz){
                if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0)
                    npoladv++;
                else
                    nsoladv++;
            }
        }

    }

    tevol[0][step]=static_cast<real>(npolrec)/(npolrec+nsolrec);
    tevol[1][step]=static_cast<real>(npoladv)/(npoladv+nsoladv);




    return;

}
