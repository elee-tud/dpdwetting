#include "main.hpp"

/****************************************************************************
 *        Program for Dissipative Particle Dynamics Simulation
 ****************************************************************************
 * Version: v3.1.2 (updated on 07. 07. 2022)
 * Written by Dr. Eunsang Lee
 * Theoretical Physical Chemistry Dept. TU Darmstadt 
 * Email: e.lee@theo.chemie.tu-darmstadt.de
 ****************************************************************************/
using namespace dpd;
int main(int argc, char* argv[]){
    SetMPI mpi(argc, argv);                 //Setting-up MPI
    Command command(argc, argv, &mpi);      //Reading command line arguments
    /*Main module*/
    if(command.getProgram()==RUN){
        Initialization init(&command);      //Initializing the program
        System system(&init);               //Loading a system class
        system.initializeSimulation();      //Initilizing the simulation setup
        system.run();                       //Running a simulation
        Finalization fin(&init, &system);   //Finalization
    }
    /*Analysis module*/
    else{
        InitAnalysis init(&command);        //Initializing the analysis program
        Calculation calc(&init);            //Calculation
    }
    mpi.finalize();                         //Finalizing MPI
	return 0;
}

