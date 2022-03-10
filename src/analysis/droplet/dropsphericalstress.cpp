#include "dropsphericalstress.hpp"
using namespace dpd;
DropSphericalStress::DropSphericalStress(InitialSet initset):Property(initset){
   
    title="   Calculation of Stress as a function of dist. from com of a dropet";
    
    nprops=3;
    for_distrib=true;
    need_position=true;
    need_multicore=true;

    outfile="spherical_stress.out";
    
    outtitle_dist="#Spherical stress tensor";
    outheader_dist=Svec{"r", "S_rr", "S_tt", "S_pp"};
}

void DropSphericalStress::getSpecificParameters(){
    dbin=0.1;
    command->getCommandSingleOption("-dr", dbin, &dbin);
    MPI_Bcast(&dbin, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    return;
}


void DropSphericalStress::initializeVariables(){
    decomp=initset.decomp;
    config=initset.config;
    inter=initset.interactions;
    force=new Force(decomp, inter);
    localdensity=new LocalDensity(decomp);
    pbc=PeriodicBoundary(config->getBox());
    

    ptcls=ParticleGroup(particles, pbc);

    ndbin=static_cast<int>(pow(box[0]*box[0]+box[1]*box[1]+box[2]*box[2], 0.5)*1.4142135/dbin);
    press=new real* [nprops];
    for(int i=0;i<nprops;i++)
        press[i]=new real[ndbin];
    ref0=Real3D{box[0]/2-1.0, box[1]/2-1.0, box[2]/2-1.0};
    ref1=Real3D{box[0]/2+1.0, box[1]/2+1.0, box[2]/2+1.0};

    initializeResultArrays();
    return;
}

void DropSphericalStress::calculateStep(int step){
    if(mpi->isMaster()){
        com=ptcls.calculateCenterOfMass(ref0, ref1);
    }
    MPI_Bcast(&com[0], 3, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    config->bcastConfiguration();
    decomp->clearAllBeadInformation();
    decomp->allocateBeadsToDomain();
    decomp->allocateBeadsToCells();
    decomp->communicateGhostBeads();
    localdensity->calculateLocalDensity();
    for(int i=0;i<nprops;i++)
        for(int j=0;j<ndbin;j++)
            press[i][j]=0.;
    force->calculateForce(com, press, dbin);
    for(int i=0;i<nprops;i++)
        for(int j=0;j<ndbin;j++)
            dist[i][j]+=press[i][j];


    

    return;
}

void DropSphericalStress::normalizeResults(){
    if(mpi->isMaster()){
        std::cout << "here" << std::endl;
        for(int i=0;i<nprops;i++){
            for(int j=0;j<ndbin;j++){
                dist_sum[i][j]/=nsteps*4*PI*(0.5+j)*(0.5+j)*dbin*dbin;
            }
        }
    }
    return;
}




