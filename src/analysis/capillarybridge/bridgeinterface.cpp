#include "bridgeinterface.hpp"
#include "../../filecontrol.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
BridgeInterface::BridgeInterface(InitialSet initset):Property(initset){
    title="   Calculation of 3D interfacial points";
    nprops=1;
    outheader_tevol=Svec{"Time", "interface"};

    for_timeevol=true;
    need_position=true;

    outfile="bridgeinterface.out";


    outtitle_tevol="#3D interfacial points";



}



void BridgeInterface::getSpecificParameters(){
    dz=1.3;
    dy=2.0;
    dcrit=6.10/4;
    dr=0.5;
    pilheight=0.;
    command->getCommandSingleOption("-dz", dz, &dz);
    command->getCommandSingleOption("-dy", dy, &dy);
    command->getCommandSingleOption("-dd", dcrit, &dcrit);
    command->getCommandSingleOption("-dr", dr, &dr);
    command->getCommandSingleOption("-ph", pilheight, &pilheight);
    surfaceb=control->getWallMinPosition()[2]+pilheight;
    surfacet=control->getWallMaxPosition()[2]-pilheight;
    dv=dz*dy*dr;
   
    liquidgrps=control->getLiquidGroups();
    return;
}


void BridgeInterface::initializeVariables(){
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
    

    numx=static_cast<int>((box[0])/dr);
    numy=static_cast<int>((box[1])/dz)+1;
    numz=static_cast<int>((surfacet-surfaceb)/dz)+1;
    density=Rvec3D(numy, Rvec2D(numz, Rvec(numx, 0)));
    interface=Rvec3D(numy, Rvec2D(numz, Rvec(2,0)));
    initializeResultArrays();
    return;
}

void BridgeInterface::calculateStep(int step){
    for(int i=0;i<numy;i++){
        for(int j=0;j<numz;j++){
            for(int k=0;k<numx;k++){
                density[i][j][k]=0.;
            }
            interface[i][j][0]=-1;
            interface[i][j][1]=-1;
        }
    }
    
    for(int i=0;i<nliqptcls;i++){
        int idxx=static_cast<int>(particles[liquididx[i]]->coord[0]/dr);
        int idxy=static_cast<int>(particles[liquididx[i]]->coord[1]/dz);
        int idxz=static_cast<int>((particles[liquididx[i]]->coord[2]-surfaceb)/dz);
        density[idxy][idxz][idxx]+=1.;
    }
    for(int i=0;i<numy;i++){
        for(int j=0;j<numz;j++){
            for(int k=0;k<numx;k++){
                density[i][j][k]/=dv;
            }
        }
    }
    int centeridx=numx/2;
    for(int i=0;i<numy;i++){
        for(int j=0;j<numz;j++){
            for(int k=centeridx-1;k>0;k--){
                if(density[i][j][k]>=dcrit && density[i][j][k-1]<dcrit){
                    interface[i][j][0]=dr*(k-(density[i][j][k]-dcrit)/(density[i][j][k]-density[i][j][k-1]));
                }
            }
            for(int k=centeridx;k<numx-1;k++){
                if(density[i][j][k]>=dcrit && density[i][j][k+1]<dcrit){
                    interface[i][j][1]=dr*(k+(density[i][j][k]-dcrit)/(density[i][j][k]-density[i][j][k+1]));
                }
            }
        }
    }
   
    writeResultsForNow(step);


    return;

}

void BridgeInterface::writeResultsForNow(int step){
    int currstep=begstep+step*skipstep;
    std::string prefix=outfile.substr(0, outfile.find("."));
    std::string outfilename=prefix+"_"+std::to_string(currstep)+".out";
    openFileWithBackup(outfilename, &outstream, false, false);
    outstream << outtitle_tevol << std::endl;
    outstream << "# Step=" << currstep << ", Time=" << currstep*dt << std::endl;;
    real y, z;
    for(int i=0;i<numy;i++){
        y=dy*(i+0.5);
        for(int j=0;j<numz;j++){
            z=dz*(j+0.5);
            if(interface[i][j][0]>0){
                outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << interface[i][j][0];
                outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << y;
                outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << z << std::endl;;
            }
        }
        outstream << std::endl;
    }
    for(int i=0;i<numy;i++){
        y=dy*(i+0.5);
        for(int j=0;j<numz;j++){
            z=dz*(j+0.5);
            if(interface[i][j][1]>0){
                outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << interface[i][j][1];
                outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << y;
                outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << z << std::endl;;
            }
        }
        outstream << std::endl;
    }
    outstream.close();
    return;

}

