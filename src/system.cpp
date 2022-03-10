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
//    if(mpi->isMaster()) std::cout << "1"<<std::endl;
        real time=i*control->getTimeStep();
//    if(mpi->isMaster()) std::cout << "2"<<std::endl;
        decomp->setUnsynced();
//    if(mpi->isMaster()) std::cout << "3"<<std::endl;
        integrator->updatePosition(need_prev_position);
//    if(mpi->isMaster()) std::cout << "4"<<std::endl;
        interaction->applyExtensionForPosition(i);
//    if(mpi->isMaster()) std::cout << "5"<<std::endl;
        integrator->updateVelocity(need_prev_velocity);
//    if(mpi->isMaster()) std::cout << "6"<<std::endl;
        decomp->refreshDomainBeads();
//    if(mpi->isMaster()) std::cout << "7"<<std::endl;
        decomp->communicateGhostBeads();
//    if(mpi->isMaster()) std::cout << "8"<<std::endl;
        if(control->getNonbondedInteraction()==NBMDPD) localdensity->calculateLocalDensity();
//    if(mpi->isMaster()) std::cout << "9"<<std::endl;
    
//        if(mpi->isMaster()) std::cout << "t=" << i << std::endl;
//        if(particles[1]->existsHere()==TRUEPTCL)
//            std::cout << mpi->rank() << std::endl;

//    if(mpi->isMaster()) std::cout << control->slipSpring() <<"," << control->getNumDPDSeqSteps() <<std::endl;
        if(control->slipSpring() && i%control->getNumDPDSeqSteps()==0){
            sspring->relocation();
//    if(mpi->isMaster()) std::cout << "10"<<std::endl;
            sspring->doMonteCarloSequence(i);
//    if(mpi->isMaster()) std::cout << "11"<<std::endl;
            /*
            if(mpi->isMaster()){
                ParticleList particles=config->getParticles();
                int ntotssbonds=0;
                for(int i=0;i<particles.size();i++){
                    if(particles[i]->getSSBonds().size()>0){
                        for(int j=0;j<particles[i]->getSSBonds().size() ;j++){
                            std::cout << "(" << i << "," << particles[i]->getSSBonds()[j]->getParticleIndex()-1 << ")" << ",";
                                ntotssbonds++;
                        }
                    }
                }
                std::cout << std::endl << "Total " << ntotssbonds << std::endl;;
            }
            */
        }

        force->calculateForce();
//    if(mpi->isMaster()) std::cout << "12"<<std::endl;
        integrator->updateVelocity();
//    if(mpi->isMaster()) std::cout << "13"<<std::endl;
        interaction->applyExtensionForVelocity(i);
//    if(mpi->isMaster()) std::cout << "14"<<std::endl;


        if(i%logfreq==0){
            analysis->calculateSystemProperties();
            dump->dumpLog(i, time, timer->getRemainingTime());
        }

        if(i%trajfreq==0){
            dump->dumpConf(i, time);
            restart->writeCheckPoint(i, time);
        }

        if(strsfreq > 0 && i%strsfreq==0)
            dump->dumpStress(i, time);

        if(frcfreq > 0 && i%frcfreq==0)
            dump->dumpForce(i, time);

        if(control->getIntegrator()==EMIN){
            timer->printProgressRemainingTime(i, integrator->getMaxForce());
            if(integrator->getMaxForce()<control->getTargetForce()){
                if(mpi->isMaster()) std::cout << "Final maximum force is " << integrator->getMaxForce() << std::endl;
                break;
            }
        }
        else
            timer->printProgressRemainingTime(i);


    }
    return;
}

