#ifndef __COMMAND__HPP
#define __COMMAND__HPP

#include "types.hpp"
#include "setmpi.hpp"
#include <string>
#include <fstream>



namespace dpd{

class Command{
protected:
    std::string config_fname;
    std::string topol_fname;
    std::string control_fname;
    std::string out_prefix;
    std::string traj_fname;
    std::string log_fname;
    std::string strs_fname;
    std::string fnconf_fname;
    std::string frc_fname;
    std::string out_fname;
    std::string ndx_fname;
    std::string ckp_fname;
    std::string rst_fname;
    std::string ss_fname;
    std::string ssgro_fname;

    void printProgramStartingMessage();

    SetMPI *mpi;
    void assignDefaults();
    int argc;
    Svec argv;

    /*varialbes for analysis*/
    int program;
    int first_step;
    int final_step;
    int skip_step;

    bool do_restart;


public:
    Command(){}
    Command(int argc, char* argv[], SetMPI *mpi);
    ~Command();
    
    void getProgramOption();
    void getCommandOption();
    void getCommandOptionAnalysis();
    void getCommandSingleOption(const char opt[], std::string deflt, std::string *value);
    void getCommandSingleOption(const char opt[], int deflt, int *value);
    void getCommandSingleOption(const char opt[], real deflt, real *value);
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

    bool doRestart(){ return do_restart;}

    SetMPI* getMPI(){ return mpi; }

    int firststep(){ return first_step; }
    int finalstep(){ return final_step; }
    int skipstep(){ return skip_step; }
    int getProgram(){ return program; }

};
};
#endif

