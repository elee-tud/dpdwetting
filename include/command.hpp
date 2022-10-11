/*************************************************************
 *           Class "Command"
 *************************************************************
 *  This class manages basic information from the command line
 * including the kind of module to run, and the input and ouput
 * file names from the flag options. It also contains a function
 * of getting special parameters for analysis.
 * **********************************************************/
#ifndef __COMMAND__HPP
#define __COMMAND__HPP

#include "types.hpp"
#include "setmpi.hpp"
#include <string>
#include <fstream>



namespace dpd{

class Command{
protected:
    std::string config_fname;           /*Initial configuration file name*/
    std::string topol_fname;            /*Topology file name*/
    std::string control_fname;          /*Control file name*/
    std::string out_prefix;             /*Prefix for outputs*/
    std::string traj_fname;             /*Trajectory file name*/
    std::string log_fname;              /*Log file name*/
    std::string strs_fname;             /*Stress file name*/
    std::string fnconf_fname;           /*Final configuration file name*/
    std::string frc_fname;              /*force file name*/
    std::string out_fname;              /*Prefix for outputs of analysis modules*/
    std::string ndx_fname;              /*Index file name*/
    std::string ckp_fname;              /*Checkpoint file name*/
    std::string rst_fname;              /*Restart file name*/
    std::string ss_fname;               /*Slip-spring log file name*/
    std::string ssgro_fname;            /*Slip-spring trajectory file name*/

    void printProgramStartingMessage();     /*Printing out a messeage when starting a program*/

    SetMPI *mpi;
    void assignDefaults();                  /*Assigning default values of file names*/
    int argc;                               /*# of arguments*/
    Svec argv;                              /*Arguments*/

    /*varialbes for analysis*/
    int program;                    /*Program Index*/
    int first_step;                 /*First step to analize*/
    int final_step;                 /*Final step to analize*/
    int skip_step;                  /*Frequency to analyze*/

    bool do_restart;                /*if the simulation is restart, true*/


public:
    Command(){}
    Command(int argc, char* argv[], SetMPI *mpi);
    ~Command();
    
    void getProgramOption();            /*Getting the module to run*/
    void getCommandOption();            /*Getting flag options from the command line*/
    void getCommandOptionAnalysis();    /*Getting flag options from the command line for analysis*/
    void getCommandSingleOption(const char opt[], std::string deflt, std::string *value);   /*Getting a flag (string) option from a command line*/
    void getCommandSingleOption(const char opt[], int deflt, int *value);   /*Getting a flag (int) option from a command line*/
    void getCommandSingleOption(const char opt[], real deflt, real *value); /*Getting a flag (float) option from a command line*/
    std::string configuration(){ return config_fname; }
    std::string topology(){ return topol_fname; }
    std::string control(){ return control_fname; }
    std::string trajectory(){ return traj_fname; }
    std::string log(){ return log_fname; }
    std::string stress(){ return strs_fname; }
    std::string force(){ return frc_fname; }
    std::string finalconf(){ return fnconf_fname; }
    std::string output_analysis(){ return out_fname; }
    std::string index(){ return ndx_fname; }
    std::string checkpoint(){ return ckp_fname; }
    std::string restart(){ return rst_fname; }
    std::string slipspring(){ return ss_fname; }
    std::string slipspringtrj(){ return ssgro_fname; }

    bool doRestart(){ return do_restart;}       /*Checking if the simulation is restarted*/

    SetMPI* getMPI(){ return mpi; }

    int firststep(){ return first_step; }
    int finalstep(){ return final_step; }
    int skipstep(){ return skip_step; }
    int getProgram(){ return program; }

};
};
#endif

