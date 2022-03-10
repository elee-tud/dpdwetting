#include "bridgeclvel.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
BridgeContLineVelocity::BridgeContLineVelocity(InitialSet initset):Property(initset){
    title="   Calculation of the contact line velocity of each component";
    nprops=8;
    outheader_tevol=Svec{"Time", "v_P(rec1)", "x_S(rec1)",  "v_P(rec2)", "x_S(rec2)",  "v_P(adv1)", "v_S(adv1)",  "v_P(adv2)", "v_S(adv2)"};

    for_timeevol=true;
    need_position=true;
    need_velocity=true;

    outfile="bridgeclvel.out";


    outtitle_tevol="#Polymer concentration at the contact line zone";


}



void BridgeContLineVelocity::getSpecificParameters(){
    dz=0.65;
    dcrit=6.10/4;
    dr=0.5;
    drtpl=1.0;
    command->getCommandSingleOption("-dz", dz, &dz);
    command->getCommandSingleOption("-dd", dcrit, &dcrit);
    command->getCommandSingleOption("-dr", dr, &dr);
    command->getCommandSingleOption("-dl", drtpl, &drtpl);
    
    surfaceb=control->getWallMinPosition()[2];
    surfacet=control->getWallMaxPosition()[2];
    shearrate=control->getWallShearRate();
   
    height=surfacet-surfaceb;
    shearvel=shearrate*height/2;
    liquidgrps=control->getLiquidGroups();
    return;
}


void BridgeContLineVelocity::initializeVariables(){
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
    prev_in_cl=Ivec(nliqptcls, 0);
    curr_in_cl=Ivec(nliqptcls, 0);
    zcenter=box[2]/2;
    

    numx=static_cast<int>((box[0])/dr);
    numz=static_cast<int>((zcenter-surfaceb-dz*0.5)/dz);
    topdensity=Rvec(numx, 0);
    botdensity=Rvec(numx, 0);
    interface=Rvec(4, 0);
    initializeResultArrays();
    return;
}

void BridgeContLineVelocity::calculateStep(int step){
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
        real distzfc=particles[liquididx[i]]->coord[2]-zcenter;
        if(distzfc<0){ 
            int idxz=static_cast<int>((particles[liquididx[i]]->coord[2]-surfaceb-0.5*dz)/dz);
            if(idxz==0){
                botdensity[idxx]+=1.;
            }

        }
        else{
            int idxz=static_cast<int>((surfacet-particles[liquididx[i]]->coord[2]-0.5*dz)/dz);
            if(idxz==0){
                topdensity[idxx]+=1.;
            }
        }
    }
    
    for(int k=0;k<numx;k++){
        topdensity[k]/=dr*box[1]*dz;
        botdensity[k]/=dr*box[1]*dz;
    }
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

    int npolrec1=0, npoladv1=0;
    int npolrec2=0, npoladv2=0;
    int nsolrec1=0, nsoladv1=0;
    int nsolrec2=0, nsoladv2=0;
    real drrec1, drrec2, dradv1, dradv2;

    for(int i=0;i<nliqptcls;i++){
        drrec1=sqrt(pow(particles[liquididx[i]]->coord[0]-interface[0], 2.)+pow(particles[liquididx[i]]->coord[2]-surfaceb-0.5*dz, 2.));
        drrec2=sqrt(pow(particles[liquididx[i]]->coord[0]-interface[3], 2.)+pow(surfacet-particles[liquididx[i]]->coord[2]-0.5*dz, 2.));
        dradv1=sqrt(pow(particles[liquididx[i]]->coord[0]-interface[1], 2.)+pow(particles[liquididx[i]]->coord[2]-surfaceb-0.5*dz, 2.));
        dradv2=sqrt(pow(particles[liquididx[i]]->coord[0]-interface[2], 2.)+pow(surfacet-particles[liquididx[i]]->coord[2]-0.5*dz, 2.));

        if(drrec1<drtpl){
            curr_in_cl[i]=1;
            if(step > 0 && prev_in_cl[i]==1){
            
                if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0){
                    npolrec1++;
                    tevol[0][step]+=particles[liquididx[i]]->veloc[0];
                }
                else{
                    nsolrec1++;
                    tevol[1][step]+=particles[liquididx[i]]->veloc[0];
                }
            }
        }
        else if(drrec2<drtpl){
            curr_in_cl[i]=2;
            if(step > 0 && prev_in_cl[i]==2){
            
                if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0){
                    npolrec2++;
                    tevol[2][step]+=particles[liquididx[i]]->veloc[0];
                }
                else{
                    nsolrec2++;
                    tevol[3][step]+=particles[liquididx[i]]->veloc[0];
                }
            }
        }
        else if(dradv1<drtpl){
            curr_in_cl[i]=3;
            if(step > 0 && prev_in_cl[i]==3){
                if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0){
                    npoladv1++;
                    tevol[4][step]+=particles[liquididx[i]]->veloc[0];
                }
                else{
                    nsoladv1++;
                    tevol[5][step]+=particles[liquididx[i]]->veloc[0];
                }
            }
        }
        else if(dradv2<drtpl){
            curr_in_cl[i]=4;
            if(step > 0 && prev_in_cl[i]==4){
                if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0){
                    npoladv2++;
                    tevol[6][step]+=particles[liquididx[i]]->veloc[0];
                }
                else{
                    nsoladv2++;
                    tevol[7][step]+=particles[liquididx[i]]->veloc[0];
                }
            }
        }
        else
            curr_in_cl[i]=0;
    }
    for(int i=0;i<nliqptcls;i++){
        prev_in_cl[i]=curr_in_cl[i];
    }

    if(npolrec1!=0)
        tevol[0][step]/=static_cast<real>(npolrec1);
    else
        tevol[0][step]=0.;
    if(nsolrec1!=0)
        tevol[1][step]/=static_cast<real>(nsolrec1);
    else
        tevol[1][step]=0.;
    if(npolrec2!=0)
        tevol[2][step]/=static_cast<real>(npolrec2);
    else
        tevol[2][step]=0.;
    if(nsolrec2!=0)
        tevol[3][step]/=static_cast<real>(nsolrec2);
    else
        tevol[3][step]=0.;
    if(npoladv1!=0)
        tevol[4][step]/=static_cast<real>(npoladv1);
    else
        tevol[4][step]=0.;
    if(nsoladv1!=0)
        tevol[5][step]/=static_cast<real>(nsoladv1);
    else
        tevol[5][step]=0.;
    if(npoladv2!=0)
        tevol[6][step]/=static_cast<real>(npoladv2);
    else
        tevol[6][step]=0.;
    if(nsoladv2!=0)
        tevol[7][step]/=static_cast<real>(nsoladv2);
    else
        tevol[7][step]=0.;




    return;

}
