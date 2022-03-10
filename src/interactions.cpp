#include "interactions.hpp"

using namespace dpd;
Interactions::Interactions(Control* control, Topology* topology, Configuration* config, Decomposition* decomp):control(control), topology(topology), config(config), decomp(decomp){
    mpi=decomp->getMPI();
    if(mpi->isMaster()) std::cout << " Interaction types adeed: ";
    addNonbondedInteraction();
    addBondedInteraction();
    addExternalInteraction();
    if(mpi->isMaster()) std::cout << std::endl;
    if(mpi->isMaster()) std::cout << " Extension added: ";
    addExtension();
    if(mpi->isMaster()) std::cout << std::endl;

}


void Interactions::addNonbondedInteraction(){
    if(control->getNonbondedInteraction()==NBDPD){
        nblist.push_back(new DPDForce(topology, config, decomp));
        if(mpi->isMaster()) std::cout << "[DPD force] ";
    }
    else if(control->getNonbondedInteraction()==NBMDPD){
        nblist.push_back(new MDPDForce(topology, config, decomp));
        if(mpi->isMaster()) std::cout << "[MDPD force] ";
    }

    return;

}

void Interactions::addBondedInteraction(){
    if(control->getBondLengthInteraction()==BLHARMONIC){
        bdlist.push_back(new BondLenHarmonic(topology, config, decomp));
        if(mpi->isMaster()) std::cout << "[Harmonic bond length] ";
    }
    if(control->slipSpring()){
        bdlist.push_back(new SlipSpringHarmonic(topology, config, decomp));
        if(mpi->isMaster()) std::cout << "[Slip spring harmonic potential] ";
    }
    return;

}



void Interactions::addExternalInteraction(){
    if(control->getPullingSpringK()!=0.0){
        exlist.push_back(new SpringPulling(topology, config, decomp));
        if(mpi->isMaster()) std::cout << "[Spring pulling toward a point] ";
    }
    if(control->getGravityDirection()!=-1){
        exlist.push_back(new Gravity(topology, config, decomp));
        if(mpi->isMaster()) std::cout << "[Gravitation field] ";
    }
    return;

}

void Interactions::addExtension(){
    if(control->getWallType()==SOLIDWALL){
        extensionlist.push_back(new SolidWall(topology, config, decomp));
        if(mpi->isMaster()) std::cout << "[Solid wall] ";
    }
    if(control->getRmCOMVelFrequency()>0){
        extensionlist.push_back(new RemovalCOMVel(topology, config, decomp));
        if(mpi->isMaster()) std::cout << "[COM motion removal] ";
    }
    if(control->getDeformation()==ELONGATION){
        extensionlist.push_back(new BoxElongation(topology, config, decomp));
        if(mpi->isMaster()) std::cout << "[Isochoric elongational deformation] ";
    }

    if(control->getPinning()){
        extensionlist.push_back(new Pinning(topology, config, decomp));
        if(mpi->isMaster()) std::cout << "[Particle pinning] ";
    }
    if(control->isUniflow()){
        extensionlist.push_back(new UniflowPBC(topology, config, decomp));
        if(mpi->isMaster()) std::cout << "[PBC for uniaxial flow] ";
    }
    if(control->isWallSheared()){
        extensionlist.push_back(new WallInducedShear(topology, config, decomp));
        if(mpi->isMaster()) std::cout << "[Wall-induced shear flow] ";
    }



    if(mpi->isMaster()){
        if(extensionlist.size()==0)
            std::cout << "None." << std::endl;
        else
            std::cout << std::endl;
    }
    return;

}

void Interactions::applyExtensionForPosition(int step){
    for(int i=0;i<extensionlist.size();i++){
        extensionlist[i]->applyExtensionForPosition(step);
    }
    return;
}

void Interactions::applyExtensionForVelocity(int step){
    for(int i=0;i<extensionlist.size();i++){
        extensionlist[i]->applyExtensionForVelocity(step);
    }
    return;
}





