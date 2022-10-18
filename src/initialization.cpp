#include "initialization.hpp"
#include <unistd.h>

using namespace dpd;
Initialization::Initialization(Command* command):command(command){
    /*Creating instances of initializing components*/
    mpi=command->getMPI();
    control=new Control(command, mpi);      
    control->setDefaults();             /*Setting-up default values for control variables*/
    control->readControl();             /*Getting control variables*/
    control->setDependency();           /*Resolving dependency of control variables*/
    control->bcastControl();            /*Broadcasting control variables*/
    topol=new Topology(command, control, mpi);
    topol->readTopology();              /*Reading topology*/
    topol->bcastTopology();             /*Broadcasting general topology*/
    topol->bcastMolecules();            /*Broadcasting molecule information*/
    topol->buildTopology();             /*Building molecule topology*/
    topol->checkDependency();           /*Checking topology dependency*/
    topol->printSystemInformation();    /*Printing out system information*/
    config=new Configuration(topol, mpi);
    config->readConfiguration();        /*Reading initial configuration*/
    config->bcastConfiguration();       /*Broadcasting initial configuration*/
    decomp=new Decomposition(control, config, mpi);     /*Domain decomposition*/
    
    if(control->slipSpring()){
        sspring=new SlipSpring(control, topol, config, decomp);     /*Setting-up slip springs*/
    }

    interactions=new Interactions(control, topol, config, decomp);  /*Setting-up interactions*/
    timer=new Timer(mpi->rank(), control->getTotalSteps());         /*Setting-up timer*/
    MPI_Barrier(MPI_COMM_WORLD);
    restart=new Restart(command, config, topol, decomp, mpi);       /*Reading restart*/

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

