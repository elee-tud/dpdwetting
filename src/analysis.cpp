#include "analysis.hpp"

using namespace dpd;
Analysis::Analysis(Initialization* init):init(init){
    decomp=init->getDecomposition();    
    ntot_beads=init->getTopology()->getNbeads();    //The number of total beads in a system
    ntot_unfrozen_beads=init->getTopology()->getNumUnfrozenBeads();     //The number of unfrozen beads in a system
    box=init->getConfiguration()->getBox();     //System box
    boxvolume=box[0]*box[1]*box[2];     //Box volume
    mpi=init->getMPI(); 
    particles=decomp->getParticles();


}
void Analysis::calculateSystemProperties(){
    calculateTemperature();
    calculateStressTensor();
    return;
}

void Analysis::calculateTemperature(){
    mybeads=decomp->getBeadsIndexInDomain();
    real temp_proc=0.;
    for(int i=0;i<mybeads.size();i++){
        if(!particles[mybeads[i]]->isFrozen())      //if particle is not frozen
            temp_proc+=particles[mybeads[i]]->getMass()*(particles[mybeads[i]]->veloc.sqr());   //kinetic energy+=mv^2
    }

    temperature=0.0;
    MPI_Reduce(&temp_proc, &temperature, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);   //Reduce the sum of kinetic energy
    MPI_Barrier(MPI_COMM_WORLD);
    if(mpi->isMaster())
        temperature/=3*ntot_unfrozen_beads;     //T=kinetic energy/(3*particle number)
    return;
}


void Analysis::calculateStressTensor(){
    mybeads=decomp->getBeadsIndexInDomain();
    Rvec virialt_proc(9,0.);
    for(int i=0;i<mybeads.size();i++){
        if(!particles[mybeads[i]]->isFrozen()){
            //Calculate stress tensor
            virialt_proc[0]+=particles[mybeads[i]]->stress[0]/2+particles[mybeads[i]]->getMass()*particles[mybeads[i]]->veloc[0]*particles[mybeads[i]]->veloc[0];   //sigma_xx
            virialt_proc[1]+=particles[mybeads[i]]->stress[1]/2+particles[mybeads[i]]->getMass()*particles[mybeads[i]]->veloc[1]*particles[mybeads[i]]->veloc[1];   //sigma_yy
            virialt_proc[2]+=particles[mybeads[i]]->stress[2]/2+particles[mybeads[i]]->getMass()*particles[mybeads[i]]->veloc[2]*particles[mybeads[i]]->veloc[2];   //sigma_zz
            virialt_proc[3]+=particles[mybeads[i]]->stress[3]/2+particles[mybeads[i]]->getMass()*particles[mybeads[i]]->veloc[0]*particles[mybeads[i]]->veloc[1];   //sigma_xy
            virialt_proc[4]+=particles[mybeads[i]]->stress[4]/2+particles[mybeads[i]]->getMass()*particles[mybeads[i]]->veloc[0]*particles[mybeads[i]]->veloc[2];   //sigma_xz
            virialt_proc[5]+=particles[mybeads[i]]->stress[5]/2+particles[mybeads[i]]->getMass()*particles[mybeads[i]]->veloc[1]*particles[mybeads[i]]->veloc[0];   //sigma_yx
            virialt_proc[6]+=particles[mybeads[i]]->stress[6]/2+particles[mybeads[i]]->getMass()*particles[mybeads[i]]->veloc[1]*particles[mybeads[i]]->veloc[2];   //sigma_yz
            virialt_proc[7]+=particles[mybeads[i]]->stress[7]/2+particles[mybeads[i]]->getMass()*particles[mybeads[i]]->veloc[2]*particles[mybeads[i]]->veloc[0];   //sigma_zx
            virialt_proc[8]+=particles[mybeads[i]]->stress[8]/2+particles[mybeads[i]]->getMass()*particles[mybeads[i]]->veloc[2]*particles[mybeads[i]]->veloc[1];   //sigma_zy


        }
    }
    virialt=Rvec(9, 0.);
    MPI_Reduce(&virialt_proc[0], &virialt[0], 9, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);  //Reduce stresses
    MPI_Barrier(MPI_COMM_WORLD);
    box=decomp->getBox();
    boxvolume=box[0]*box[1]*box[2];
    if(mpi->isMaster()){
        for(int i=0;i<9;i++){
            virialt[i]/=boxvolume;      //calculate stress tensor
        }
        pressure=(virialt[0]+virialt[1]+virialt[2])/3;      //pressure=(sigma_xx+sigma_yy+sigma_zz)/3

    }
    return;
}


