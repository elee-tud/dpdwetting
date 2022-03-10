#include "main.hpp"

/****************************************************************************
 *        Program for Dissipative Particle Dynamics Simulation
 ****************************************************************************
 * Version: v3.1.0 (updated on 10.Mar 2022)
 * Written by Dr. Eunsang Lee
 * Theoretical Physical Chemistry Dept. TU Darmstadt 
 * Email: e.lee@theo.chemie.tu-darmstadt.de
 ****************************************************************************/
using namespace dpd;
int main(int argc, char* argv[]){
    SetMPI mpi(argc, argv);
    Command command(argc, argv, &mpi);
    if(command.getProgram()==RUN){
        Initialization init(&command);    //Initializing the program
        System system(&init);               //Loading a system class
        system.initializeSimulation();      //Initilizing the simulation setup
        system.run();                       //Running a simulation
        Finalization fin(&init, &system);
    }
    else{
        InitAnalysis init(&command);
        Calculation calc(&init);
    }
    mpi.finalize();
	return 0;
}

