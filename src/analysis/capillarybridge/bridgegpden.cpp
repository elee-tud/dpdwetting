#include "bridgegpden.hpp" 
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
BridgeGrooveParticleDensity::BridgeGrooveParticleDensity(InitialSet initset):Property(initset){
    title="   Calculation of the number of particles in a groove";
    nprops=6;
    outheader_tevol=Svec{"Time", "n_S", "n_P", "n_tot", "d_S", "d_P", "d_tot"};

    for_timeevol=true;
    need_position=true;

    outfile="bridgegrvdens.out";


    outtitle_tevol="#The number of particles in a groove";
    dt=control->getTimeStep()*control->getTrajFrequency();


}



void BridgeGrooveParticleDensity::getSpecificParameters(){
    dz=1.3;
    dcrit=6.10/4;
    dr=0.5;
    pilheight=0.;
    pilgap=2.;
    pilgap=1.0;
    command->getCommandSingleOption("-dz", dz, &dz);
    command->getCommandSingleOption("-dd", dcrit, &dcrit);
    command->getCommandSingleOption("-dr", dr, &dr);
    command->getCommandSingleOption("-ph", pilheight, &pilheight);
    command->getCommandSingleOption("-pw", pilwidth, &pilwidth);
    command->getCommandSingleOption("-pg", pilgap, &pilgap);
    grooveb=surfaceb;
    groovet=surfacet;
    surfaceb=control->getWallMinPosition()[2]+pilheight;
    surfacet=control->getWallMaxPosition()[2]-pilheight;
    pilperiod=pilwidth+pilgap; 
    liquidgrps=control->getLiquidGroups();

    srate=control->getWallShearRate();
    platevel=srate*(groovet-grooveb);

    return;
}


void BridgeGrooveParticleDensity::initializeVariables(){
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
    

    topdensity=Rvec(numx, 0);
    botdensity=Rvec(numx, 0);
    interface=Rvec(4, 0);
    initializeResultArrays();
    return;
}

void BridgeGrooveParticleDensity::calculateStep(int step){
    interface[0]=0.;
    interface[1]=0.;
    interface[2]=0.;
    interface[3]=0.;
    for(int k=0;k<numx;k++){
        topdensity[k]=0.;
        botdensity[k]=0.;
    }
    for(int i=0;i<nliqptcls;i++){
        int idxx=static_cast<int>(particles[liquididx[i]]->coord[0]/dr);
        int idxz=static_cast<int>((particles[liquididx[i]]->coord[2]-surfaceb-0.5*dz)/dz);
        if(idxz==0){
                botdensity[idxx]+=1.;
        }
        idxz=static_cast<int>((surfacet-particles[liquididx[i]]->coord[2]-0.5*dz)/dz);
        if(idxz==0){
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
    for(int k=0;k<numx-1;k++){
        if(botdensity[k]<dcrit && botdensity[k+1]>=dcrit){
            interface[0]=k*dr+dr/(botdensity[k+1]-botdensity[k])*(dcrit-botdensity[k]);
                break;
        }
    }
    for(int k=0;k<numx-1;k++){
        if(topdensity[k]<dcrit && topdensity[k+1]>=dcrit){
            interface[2]=k*dr+dr/(topdensity[k+1]-topdensity[k])*(dcrit-topdensity[k]);
            break;
        }
    }
    for(int k=numx-1;k>0;k--){
        if(botdensity[k]<dcrit && botdensity[k-1]>=dcrit){
            interface[1]=(k-1)*dr+dr/(botdensity[k]-botdensity[k-1])*(dcrit-botdensity[k-1]);
            break;
        }
    }
    for(int k=numx-1;k>0;k--){
        if(topdensity[k]<dcrit && topdensity[k-1]>=dcrit){
            interface[3]=(k-1)*dr+dr/(topdensity[k]-topdensity[k-1])*(dcrit-topdensity[k-1]);
            break;
        }
    }
//    std::cout << interface << std::endl;

    for(int i=0;i<nliqptcls;i++){
        if(particles[liquididx[i]]->coord[2]>=grooveb && particles[liquididx[i]]->coord[2] <surfaceb && particles[liquididx[i]]->coord[0]>=interface[0] && particles[liquididx[i]]->coord[0]<interface[1]){
            if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0)
                tevol[0][step]+=1.;
            else if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0)
                tevol[1][step]+=1.;
        }
        else if(particles[liquididx[i]]->coord[2]<groovet && particles[liquididx[i]]->coord[2] >=surfacet && particles[liquididx[i]]->coord[0]>=interface[2] && particles[liquididx[i]]->coord[0]<interface[3]){
            if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0)
                tevol[0][step]+=1.;
            else if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0)
                tevol[1][step]+=1.;
        }
    }
    tevol[2][step]=tevol[0][step]+tevol[1][step];

    pilpos0=fmod(platevel*(begstep+step)*dt, pilperiod);
//    std::cout << begstep << "," << step <<"," << dt << "," << pilpos0 << std::endl;

    int ngrvst0=int((interface[0]+pilpos0)/pilperiod);
    if(fmod(interface[0]+pilpos0, pilperiod)>pilgap)
        ngrvst0++;
    int ngrvst1=int((interface[1]+pilpos0)/pilperiod);
    if(fmod(interface[1]+pilpos0, pilperiod)>pilgap)
        ngrvst1++;
    int ngrvsb=ngrvst1-ngrvst0;
    ngrvst0=int((interface[2]-pilpos0)/pilperiod);
    if(fmod(interface[2]-pilpos0, pilperiod)>pilgap)
        ngrvst0++;
    ngrvst1=int((interface[3]-pilpos0)/pilperiod);
    if(fmod(interface[3]-pilpos0, pilperiod)>pilgap)
        ngrvst1++;
    int ngrvst=ngrvst1-ngrvst0;

    int nygrooves=box[1]/pilperiod;
    real groovevol=((interface[1]-interface[0])*nygrooves+(ngrvsb*(box[1]-pilgap*nygrooves)))*pilgap+((interface[3]-interface[2])*nygrooves+(ngrvst*(box[1]-pilgap*nygrooves)))*pilgap;
    tevol[3][step]=tevol[0][step]/groovevol;
    tevol[4][step]=tevol[1][step]/groovevol;
    tevol[5][step]=tevol[2][step]/groovevol;









    return;

}
