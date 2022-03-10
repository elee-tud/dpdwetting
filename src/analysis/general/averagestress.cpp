#include "averagestress.hpp"


using namespace dpd;

AverageStress::AverageStress(InitialSet initset):Property(initset){
    title="   Calculation of average stress along x, y, and z direction";

    nprops=6;
    for_timeevol=true;
    need_stress=true;

    outfile="average_stress.out";


    outtitle_tevol="#Average Stress versus time";
    outheader_tevol=Svec{"Time", "s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_zx"};



}


void AverageStress::getSpecificParameters(){
    return;
}

void AverageStress::initializeVariables(){
    initializeResultArrays();
    return;
}

void AverageStress::calculateStep(int step){
    for(int i=0;i<nbeads;i++){
        for(int j=0;j<6;j++){
            tevol[j][step]+=particles[i]->stress[j];
        }
    }
    return;
}


