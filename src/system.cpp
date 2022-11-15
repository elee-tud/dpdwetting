#include "system.hpp"
#include <unistd.h>

using namespace dpd;
System::System(Initialization* init):init(init){

    mpi=init->getMPI();
    command=init->getCommand();
    topology=init->getTopology();
    config=init->getConfiguration();
    control=init->getControl();
    decomp=init->getDecomposition();
    interaction=init->getInteractions();
    particles=config->getParticles();
    restart=init->getRestart();
    setIntegrator();
    analysis=new Analysis(init);
    dump=new Dump(init, analysis);
    force=new Force(init);
    if(control->getNonbondedInteraction()==NBMDPD) localdensity=new LocalDensity(init);
    timer=init->getTimer();
    if(control->slipSpring()){
        sspring=init->getSlipSpring();
        if(mpi->isMaster()) std::cout << " Slip-spring algorithm loaded" << std::endl;
    }

    trajfreq=control->getTrajFrequency();
    logfreq=control->getLogFrequency();
    strsfreq=control->getStressFrequency();
    frcfreq=control->getForceFrequency();

    extensionlist=interaction->getExtension();

    resolveDependencyForExtensions();
    startstep=0;
    starttime=0.;
   
}

System::~System(){
    delete dump;
    delete force;
    delete integrator;
    if(control->getNonbondedInteraction()==NBMDPD) delete localdensity;
}

void System::setIntegrator(){
    if(control->getIntegrator()==VELVER){
        integrator=new VelocityVerletDPD(init);
        if(mpi->isMaster()) std::cout << " Velocity-Verlet integrator loaded" << std::endl;
    }
    else if(control->getIntegrator()==EMIN){
        integrator=new EnerMinIntegrator(init);
        if(mpi->isMaster()) std::cout <<" Energy minimization integrator loaded" << std::endl;
    }
    else if(control->getIntegrator()==SLLOD){
        integrator=new SllodDPD(init);
        if(mpi->isMaster()) std::cout <<" DPD Velocity-Verlet with SLLOD eqn. loaded" << std::endl;
    }

}

void System::resolveDependencyForExtensions(){
    need_prev_position=false;
    need_prev_velocity=false;
    for(int i=0;i<extensionlist.size();i++){
        if(extensionlist[i]->needPrevPosition()){
            need_prev_position=true;
        }
        if(extensionlist[i]->needPrevVelocity()){
            need_prev_velocity=true;
        }
    }
}


void System::initializeSimulation(){
    if(control->getRestart()){
        restart->readCheckPoint(&startstep, &starttime);
        config->bcastConfiguration();
        MPI_Bcast(&startstep, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
        MPI_Bcast(&starttime, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    }
    decomp->setUnsynced();
    decomp->allocateBeadsToDomain();
    decomp->allocateBeadsToCells();
    decomp->communicateGhostBeads();
//    decomp->printCellInformation();
    
    if(control->getNonbondedInteraction()==NBMDPD) localdensity->calculateLocalDensity();
   
    if(control->slipSpring())
        sspring->initializeSlipSpring();

    force->calculateForce();
    if(!control->getRestart()) dump->dumpConf(0,0);
    if(control->getRestart() && startstep%frcfreq==0){}
    else dump->dumpForce(0,0);
    if(control->getRestart() && startstep%strsfreq==0){}
    else dump->dumpStress(0,0);
    if(control->getRestart()){}
    else{
        analysis->calculateSystemProperties();
        dump->dumpLog(0,0, timer->getRemainingTime());
    }
    if(mpi->rank()==0){
        std::cout << " Temperature set to T= "<< control->getTemperature()<<std::endl; ;
        if(control->getTempAnnealRate()!=0.)
            std::cout << " Temperature annealed with the rate "<< control->getTempAnnealRate()<<std::endl;
        std::cout << " Initial system state: T= "<< analysis->getTemperature() ;
        std::cout << ",   P= " << analysis->getPressure();
        std::cout << std::endl<<std::endl;
        std::cout << " Starting simulation." << std::endl;
    }
    return;

}

void System::run(){

    int nsteps=control->getTotalSteps();
    for(int i=startstep+1;i<=nsteps;i++){
        /*Current time*/
        real time=i*control->getTimeStep();
        /*Setting the communication flag unsynced*/
        decomp->setUnsynced();
        /*Updating position*/
        integrator->updatePosition(need_prev_position);
        /*Applying extensions for position*/
        interaction->applyExtensionForPosition(i);
        /*Applying extensions for velocity*/
        integrator->updateVelocity(need_prev_velocity);
        /*Refreshing domain bead information*/
        decomp->refreshDomainBeads();
        /*Communicating ghost beads*/
        decomp->communicateGhostBeads();
        /*Calculting local particle density if multi-body DPD*/
        if(control->getNonbondedInteraction()==NBMDPD) localdensity->calculateLocalDensity();
        /*Updating slip-springs if slip-spring is activated*/ 
        if(control->slipSpring() && i%control->getNumDPDSeqSteps()==0){
            sspring->relocation();
            sspring->doMonteCarloSequence(i);
        }
        /*Updating tempearture if the target temperature changes*/
        control->annealingTemperature();

        /*Calculating particle forces*/
        force->calculateForce();
        /*Updating velocities*/
        integrator->updateVelocity();
        /*Applyting extensions for velocities*/
        interaction->applyExtensionForVelocity(i);

        /*Writing logs*/
        if(i%logfreq==0){
            analysis->calculateSystemProperties();
            dump->dumpLog(i, time, timer->getRemainingTime());
        }

        /*Writing configuration trajectory*/
        if(i%trajfreq==0){
            dump->dumpConf(i, time);
            restart->writeCheckPoint(i, time);
        }

        /*Writing virial*/
        if(strsfreq > 0 && i%strsfreq==0)
            dump->dumpStress(i, time);

        /*Writing force*/
        if(frcfreq > 0 && i%frcfreq==0)
            dump->dumpForce(i, time);

        /*Cheking convergence if the integrator is energy minimization*/
        if(control->getIntegrator()==EMIN){
            timer->printProgressRemainingTime(i, integrator->getMaxForce());
            if(integrator->getMaxForce()<control->getTargetForce()){
                if(mpi->isMaster()) std::cout << "Final maximum force is " << integrator->getMaxForce() << std::endl;
                break;
            }
        }
        /*Printing out the timer on a screen*/
        else
            timer->printProgressRemainingTime(i);


    }
    return;
}

