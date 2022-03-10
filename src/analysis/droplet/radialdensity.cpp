#include "radialdensity.hpp"


using namespace dpd;

RadialDensity::RadialDensity(InitialSet initset):Property(initset){
    title= "  Calculation radial particle density";
    
    nprops=3;
    for_distrib=true;
    need_position=true;

    outfile="rdensity.out";

    outtitle_dist="#Radial solvent, polymer and all particle densities from COM";
    outheader_dist=Svec{"r", "rho(r)", "rho_S(r)", "rho_P(r)"};

}


void RadialDensity::getSpecificParameters(){
    dbin=0.1;
    command->getCommandSingleOption("-dr", dbin, &dbin);
    return;
}


void RadialDensity::initializeVariables(){
    ptcls=ParticleGroup(particles, pbc);
    ref0=Real3D{box[0]/2-1.0, box[1]/2-1.0, box[2]/2-1.0};
    ref1=Real3D{box[0]/2+1.0, box[1]/2+1.0, box[2]/2+1.0};

    ndbin=static_cast<int>(pow(box[0]*box[0]+box[1]*box[1]+box[2]*box[2], 0.5)*1.4142135/dbin);
    initializeResultArrays();
    return;

}

void RadialDensity::calculateStep(int step){
    com=ptcls.calculateCenterOfMass(ref0, ref1);
    for(int i=0;i<particles.size();i++){
        int idx=static_cast<int>((particles[i]->coord-com).abs()/dbin);
        if(particles[i]->getMoleculeName().compare("POL")==0)
            dist[2][idx]+=1.0;
        else if(particles[i]->getMoleculeName().compare("SOL")==0)
            dist[1][idx]+=1.0;
        dist[0][idx]+=1.0;
    }
    return;
}

void RadialDensity::normalizeResults(){
    for(int j=0;j<ndbin;j++){
        dist_sum[0][j]/=nsteps*4*PI*(0.5+j)*(0.5+j)*dbin*dbin*dbin;
        dist_sum[1][j]/=nsteps*4*PI*(0.5+j)*(0.5+j)*dbin*dbin*dbin;
        dist_sum[2][j]/=nsteps*4*PI*(0.5+j)*(0.5+j)*dbin*dbin*dbin;
    }
    return;
}




