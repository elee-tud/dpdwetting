#include "finalization.hpp"

using namespace dpd;
Finalization::Finalization(Initialization* init, System* system):init(init), system(system){
    config=init->getConfiguration();
    control=init->getControl();
    dump=system->getDump();
    command=init->getCommand();
    dump->dumpConf(control->getTotalSteps(), control->getTotalSteps()*control->getTimeStep(), command->finalconf());
    printFinalInfo();
    timer=init->getTimer();
    timer->printFinalTime();
}

void Finalization::printFinalInfo(){
    if(init->getMPI()->isMaster()){
        std::cout << " [Output Files]" << std::endl;
        std::cout << " Configuration trajectory: "<<  command->trajectory() << std::endl;
        std::cout << " Temperature and pressure: "<<  command->log()<< std::endl;
        std::cout << " Particle force: "<<  command->force()<< std::endl;
        std::cout << " Particle Virial tensor: "<<  command->stress()<< std::endl;
        std::cout << " Final output configuration: "<<  command->finalconf()<< std::endl;
    }
    return;
}

    


