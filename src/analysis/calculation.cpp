#include "calculation.hpp"

using namespace dpd;

Calculation::Calculation(InitAnalysis* init):init(init){
    prop=init->getProperty();
    traj=init->getTrajectory();
    mpi=init->getMPI();
    control=init->getControl();
    

    nbeads=init->getNbeads();
    nsteps=init->getNsteps();
    begstep=init->getBeginStep();
    endstep=init->getEndStep();
    skipstep=init->getSkipStep();
    timer=init->getTimer();
    traj->openTrajectory();
    calculateProperty();
    traj->closeTrajectory();
    prop->reduceResults();
    prop->normalizeResults();
    prop->generateOutFileName();
    prop->writeOutput();

}

void Calculation::calculateProperty(){
    if(mpi->isMaster()) std::cout << "Starting Analysis..." << std::endl;
    for(int i=0;i<begstep;i++){
        if(!traj->skipTrajectoryStep()&&i<endstep){
            if(mpi->isMaster()) std::cout << "Trajectory length, " << i+1 << ", is shorter than given, " << endstep+1 << "."<< std::endl;
            break;
        }
    }
    for(int i=begstep;i<endstep+1;i++){
        if(i%skipstep!=0){
            if(!traj->skipTrajectoryStep()&&i<endstep){
                if(mpi->isMaster()) std::cout << "Trajectory length, " << i+1 << ", is shorter than given, " << endstep+1 << "."<< std::endl;
                break;
            }
        }
        else{
            if(!traj->readTrajectoryStep()&&i<endstep){
                if(mpi->isMaster()) std::cout << "Trajectory length, " << i+1 << ", is shorter than given, " << endstep+1 << "."<< std::endl;
                break;
            }
            
            prop->resetBox(traj->getBox());
            prop->calculateStep((i-begstep)/skipstep);
            if(prop->getProperty()==TRJTOGRO)
                prop->dumpGro(i, i*skipstep*control->getTimeStep()*control->getTrajFrequency());

            timer->printProgressRemainingTime((i-begstep)/skipstep+1);

        }
    }
    return;
}




