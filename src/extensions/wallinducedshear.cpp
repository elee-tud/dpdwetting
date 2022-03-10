#include "wallinducedshear.hpp"

using namespace dpd;

WallInducedShear::WallInducedShear(Topology *topol, Configuration* config, Decomposition* decomp):Extension(topol, config, decomp){
    need_prev_position=false;
    need_prev_velocity=false;
    for_position=true;
    for_velocity=false;
    sheardir=control->getWallShearDir();
    shearrate=control->getWallShearRate();
    group1=control->getShearWallGroup1();
    group2=control->getShearWallGroup2();
    pbc=PeriodicBoundary(decomp->getBox());
    
    findGroupIndex();
    if(sheardir==0){
        vdir=0;
        hdir=1;
    }
    else if(sheardir==1){
        vdir=0;
        hdir=2;
    }
    else if(sheardir==2){
        vdir=1;
        hdir=0;
    }
    else if(sheardir==3){
        vdir=1;
        hdir=2;
    }
    else if(sheardir==4){
        vdir=2;
        hdir=0;
    }
    else if(sheardir==5){
        vdir=2;
        hdir=1;
    }
    dt=control->getTimeStep();
    botwallpos=control->getWallMinPosition()[hdir];
    topwallpos=control->getWallMaxPosition()[hdir];
    height=topwallpos-botwallpos;
    shearvel=shearrate*height/2;
    displ=shearvel*dt;
    assignVelocity();
}

void WallInducedShear::findGroupIndex(){
    if(mpi->isMaster()){
        for(int i=0;i<particles.size();i++){
            if(particles[i]->getMoleculeName().compare(group1)==0)
                group1idx.push_back(i);
            else if(particles[i]->getMoleculeName().compare(group2)==0)
                group2idx.push_back(i);
        }
        group1size=group1idx.size();
        group2size=group2idx.size();
    }
    
    MPI_Bcast(&group1size, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&group2size, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    if(!mpi->isMaster()){
        group1idx=Ivec(group1size, -1);
        group2idx=Ivec(group2size, -1);
    }

    MPI_Bcast(&group1idx[0], group1size, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&group2idx[0], group2size, MPI_INT, MASTER, MPI_COMM_WORLD);


    return;
}

void WallInducedShear::assignVelocity(){
    for(int i=0;i<group1size;i++){
        for(int j=0;j<3;j++){
            if(j==vdir)
                particles[group1idx[i]]->veloc[j]=-shearvel;
            else
                particles[group1idx[i]]->veloc[j]=0;
        }
    }
    for(int i=0;i<group2size;i++){
        for(int j=0;j<3;j++){
            if(j==vdir)
                particles[group2idx[i]]->veloc[j]=shearvel;
            else
                particles[group2idx[i]]->veloc[j]=0;
        }
    }
    return;
}


void WallInducedShear::applyExtensionForPosition(int step){
    for(int i=0;i<group1size;i++){
        particles[group1idx[i]]->coord[vdir]-=displ;
        particles[group1idx[i]]->coord=pbc.getVectorIntoBox(particles[group1idx[i]]->coord);


    }
    for(int i=0;i<group2size;i++){
        particles[group2idx[i]]->coord[vdir]+=displ;
        particles[group2idx[i]]->coord=pbc.getVectorIntoBox(particles[group2idx[i]]->coord);
    }
    return;
}




