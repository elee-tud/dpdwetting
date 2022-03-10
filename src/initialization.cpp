#include "initialization.hpp"
#include <unistd.h>

using namespace dpd;
Initialization::Initialization(Command* command):command(command){
    /*Creating instances of initializing components*/
    mpi=command->getMPI();
    control=new Control(command, mpi);
    control->setDefaults();
    control->readControl();
    control->setDependency();
    control->bcastControl();
    topol=new Topology(command, control, mpi);
    topol->readTopology();
    topol->bcastTopology();
    topol->bcastMolecules();
    topol->buildTopology();
    topol->checkDependency();
    topol->printSystemInformation();
    config=new Configuration(topol, mpi);
    config->readConfiguration();
    config->bcastConfiguration();
    decomp=new Decomposition(control, config, mpi);
    
    if(control->slipSpring()){
        sspring=new SlipSpring(control, topol, config, decomp);
    }

    interactions=new Interactions(control, topol, config, decomp);
    timer=new Timer(mpi->rank(), control->getTotalSteps());
    MPI_Barrier(MPI_COMM_WORLD);
    restart=new Restart(command, config, topol, decomp, mpi);

}



Initialization::~Initialization(){
    /*Deleting instances*/
    if(control->slipSpring())
        delete sspring;
    delete restart;
    delete timer;
    delete interactions;
    delete decomp;
    delete config;
    delete topol;
    delete control;
}

