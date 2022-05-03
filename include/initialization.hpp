/*************************************************************
 *           Class Definition "Initiailization"
 *************************************************************
 *  This class initializes the program getting command inputs, 
 * reading initial configurations, topology of the system, 
 * setting MPI, interaction types and domain decomposition.
 * **********************************************************/

#ifndef __INITIALIZATION__HPP
#define __INITIALIZATION__HPP


#include "command.hpp"
#include "control.hpp"
#include "topology.hpp"
#include "configuration.hpp"
#include "interactions.hpp"
#include "timer.hpp"
#include "restart.hpp"
#include "slipspring.hpp"

namespace dpd{


class Initialization{
protected:
    SetMPI *mpi;
    Command *command;
    Topology *topol;
    Control *control;
    Configuration *config;
    Decomposition* decomp;
    Interactions* interactions;
    Timer* timer;
    Restart* restart;
    SlipSpring* sspring;

   


public:
    Initialization(){}
    Initialization(Command* command);
    ~Initialization();

    /*Functions to get a pointer of each instance*/
    SetMPI* getMPI(){ return mpi; }
    Command* getCommand(){ return command; }
    Topology* getTopology(){ return topol;}
    Control* getControl(){ return control; }
    Configuration* getConfiguration(){ return config; }
    Decomposition* getDecomposition(){ return decomp; }
    Timer* getTimer(){ return timer; }
    Interactions* getInteractions(){ return interactions; }
    Restart* getRestart(){return restart; }
    SlipSpring* getSlipSpring(){ return sspring; }




};

typedef struct{
    Command* command;
    Control* control;
    Topology* topol;
    Configuration* config;
    SetMPI* mpi;
    Decomposition* decomp;
    Interactions* interactions;
    Restart* restart;
    SlipSpring* sspring;
}InitialSet;
};

#endif
