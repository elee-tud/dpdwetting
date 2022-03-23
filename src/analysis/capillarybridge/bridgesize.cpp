#include "bridgesize.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
BridgeSize::BridgeSize(InitialSet initset):Property(initset){
    title="   Calculation of a contact angle of a bridge." ;
    nprops=4;
    outheader_tevol=Svec{"Time", "CA_bl", "CA_tl", "CA_br", "CA_tr"};

    for_timeevol=true;
    need_position=true;

    outfile="bridgesize.out";


    outtitle_tevol="#Advancing and Receding contact angle of a bridge";


}



void BridgeSize::getSpecificParameters(){
    dz=0.8;
    dcrit=6.10/4;
    dr=0.5;
    pilheight=0.;
    pmorder=2;
    nfpts=8;

    command->getCommandSingleOption("-dz", dz, &dz);
    command->getCommandSingleOption("-dd", dcrit, &dcrit);
    command->getCommandSingleOption("-dr", dr, &dr);
    /*Considering the pillar height*/
    command->getCommandSingleOption("-ph", pilheight, &pilheight);
    command->getCommandSingleOption("-po", pmorder, &pmorder);
    command->getCommandSingleOption("-np", nfpts, &nfpts);
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
    numz=static_cast<int>((zcenter-surfaceb)/dz);
    if(nfpts>numz){
        nfpts=numz;
        std::cout << "The number of fitting points is larger than available." << std::endl;
        std::cout << "It is set to the maximum available number, " << numz << "." << std::endl;
    }
    topdensity=Rvec2D(nfpts, Rvec(numx, 0));
    botdensity=Rvec2D(nfpts, Rvec(numx, 0));
//    ellipsepts=Rvec2D(4*nfpts, Rvec(2,0));
    initializeResultArrays();
    return;
}

void BridgeSize::calculateStep(int step){
    for(int i=0;i<nfpts;i++){
        for(int k=0;k<numx;k++){
            topdensity[i][k]=0.;
            botdensity[i][k]=0.;
        }
    }
    for(int i=0;i<nliqptcls;i++){
        int idxx=static_cast<int>(particles[liquididx[i]]->coord[0]/dr);
        real distzfc=particles[liquididx[i]]->coord[2]-zcenter;
        if(distzfc<0){ 
//            int idxz=static_cast<int>((particles[liquididx[i]]->coord[2]-surfaceb-0.5*dz)/dz);
            int idxz=static_cast<int>((particles[liquididx[i]]->coord[2]-surfaceb)/dz);
            if(idxz>=0 && idxz<nfpts)
                botdensity[idxz][idxx]+=1.;
        }
        else{
//            int idxz=static_cast<int>((surfacet-particles[liquididx[i]]->coord[2]-0.5*dz)/dz);
            int idxz=static_cast<int>((surfacet-particles[liquididx[i]]->coord[2])/dz);
            if(idxz>=0 && idxz<nfpts)
                topdensity[idxz][idxx]+=1.;
        }
    }
    for(int i=0;i<nfpts;i++){
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
    int centerxidx=numx/2;
    Rvec2D xinterf(4, Rvec(nfpts, 0));
    Rvec2D yinterf(4, Rvec(nfpts, 0));
    for(int i=0;i<nfpts;i++){
        //Interface bottom left
        for(int k=centerxidx;k>0;k--){
            if(botdensity[i][k]>=dcrit && botdensity[i][k-1]<dcrit){
                xinterf[0][i]=dr*(k-(botdensity[i][k]-dcrit)/(botdensity[i][k]-botdensity[i][k-1]));
                yinterf[0][i]=surfaceb+(i+0.5)*dz;
//                ellipsepts[i*4][0]=dr*(k-(botdensity[i][k]-dcrit)/(botdensity[i][k]-botdensity[i][k-1]));
//                ellipsepts[i*4][1]=surfaceb+(i+0.5)*dz;
                break;
            }
        }
        //Interface top left
        for(int k=centerxidx;k>0;k--){
            if(topdensity[i][k]>=dcrit && topdensity[i][k-1]<dcrit){
                xinterf[1][i]=dr*(k-(topdensity[i][k]-dcrit)/(topdensity[i][k]-topdensity[i][k-1]));
                yinterf[1][i]=surfacet-(i+0.5)*dz;
//                ellipsepts[i*4+1][0]=dr*(k-(topdensity[i][k]-dcrit)/(topdensity[i][k]-topdensity[i][k-1]));
//                ellipsepts[i*4+1][1]=surfacet-(i+0.5)*dz;
                break;
            }
        }
        //Interface bottom right
        for(int k=centerxidx;i<numx-1;k++){
            if(botdensity[i][k]>=dcrit && botdensity[i][k+1]<dcrit){
                xinterf[2][i]=dr*(k+(botdensity[i][k]-dcrit)/(botdensity[i][k]-botdensity[i][k+1]));
                yinterf[2][i]=surfaceb+(i+0.5)*dz;
//                ellipsepts[i*4+2][0]=dr*(k+(botdensity[i][k]-dcrit)/(botdensity[i][k]-botdensity[i][k+1]));
//                ellipsepts[i*4+2][1]=surfaceb+(i+0.5)*dz;
                break;
            }
        }
        //Interface top right
        for(int k=centerxidx;i<numx-1;k++){
            if(topdensity[i][k]>=dcrit && topdensity[i][k+1]<dcrit){
                xinterf[3][i]=dr*(k+(topdensity[i][k]-dcrit)/(topdensity[i][k]-topdensity[i][k+1]));
                yinterf[3][i]=surfacet-(i+0.5)*dz;
//                ellipsepts[i*4+3][0]=dr*(k+(topdensity[i][k]-dcrit)/(topdensity[i][k]-topdensity[i][k+1]));
//                ellipsepts[i*4+3][1]=surfacet-(i+0.5)*dz;
                break;
            }
        }

    }

    /*
    for(int i=0;i<numz;i++){
        std::cout << ellipsepts[i*4][0] << " " << ellipsepts[i*4][1] <<std::endl;
        std::cout << ellipsepts[i*4+1][0] << " " << ellipsepts[i*4+1][1] <<std::endl;
        std::cout << ellipsepts[i*4+2][0] << " " << ellipsepts[i*4+2][1] <<std::endl;
        std::cout << ellipsepts[i*4+3][0] << " " << ellipsepts[i*4+3][1] <<std::endl;
    }
    real center_x, center_y, phi, width, height;
    ellipse_fit interfit;
    interfit.set(ellipsepts);
    interfit.fit(center_x, center_y, phi, width, height);
    std::cout << center_x << "," << center_y << "," << phi << "," << width << "," << height << std::endl;
    */

    for(int i=0;i<4;i++){
        PolynomFit polfit(pmorder, xinterf[i], yinterf[i]);
        polfit.fit();
        Rvec coeff=polfit.getResults();
        real xmid=(xinterf[i][0]+xinterf[i][1])/2;
        real slope=0;
        for(int j=1;j<pmorder+1;j++){
            slope+=j*coeff[j]*pow(xmid, j-1);
        }

        int slopesign=sign(slope);
        int delx=1;
        if(i==1||i==2)
            delx=-1;
        tevol[i][step]=atan2(slope*slopesign, delx*slopesign)*180/PI;


    }





    return;

}
