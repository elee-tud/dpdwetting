#include "rdf.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
RadialDistributionFunction::RadialDistributionFunction(InitialSet initset):Property(initset){
    title="   Calculation of radial distribution function";

    nprops=1;
    for_distrib=true;
    need_position=true;

    outfile="rdf.out";


    outtitle_dist="#Raidal distribution function";
    outheader_dist=Svec{"r", "g(r)"};

}






void RadialDistributionFunction::initializeVariables(){

    maxr=box.abs()/2.;
    dbin=0.01;
    dbegin=0.0;
    ndbin=int(maxr/dbin)+1;

            
    initializeResultArrays();
    if(particles.size()>MAXRDFAVG)
        navg=MAXRDFAVG;
    else
        navg=particles.size();
    return;
}

void RadialDistributionFunction::calculateStep(int step){
    for(int i=0;i<navg;i++){
        for(int j=0;j<particles.size();j++){
            if(i!=j){
                real r=pbc.getMinimumImageVector(particles[i]->coord, particles[j]->coord).abs();
                dist[0][int(r/dbin)]+=1.;
            }
        }
    }
    return;

}
void RadialDistributionFunction::normalizeResults(){
    for(int i=0;i<ndbin;i++){
        dist_sum[0][i]/=nsteps*navg*4*PI*dbin*dbin*dbin*(i+0.5)*(i+0.5)*particles.size()/box[0]/box[1]/box[2];
    }
    return;
}







