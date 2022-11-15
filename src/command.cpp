#include "command.hpp"
#include <cstring>
#include <algorithm>

using namespace dpd;
Command::Command(int argc, char* _argv[], SetMPI* _mpi):argc(argc), mpi(_mpi){
    argv=Svec(argc);
    for(int i=0;i<argc;i++){
        argv[i]=std::string(_argv[i]);
    }

    setAnalysisProgram();       //Assigning the name of program
    getProgramOption();         //Getting program name
    printProgramStartingMessage();      


}

void Command::printProgramStartingMessage(){
    if(mpi->rank()==0){
        std::cout << "***************************************************************************"<< std::endl;

        if(program==RUN)
            std::cout << "*        Program for Dissipative Particle Dynamics Simulation              "<< std::endl;
        else
            std::cout << "*        Program for Analysing DPD trajectory                              "<< std::endl;
        std::cout << "***************************************************************************"<< std::endl;

        std::cout << "* Version: v3.1 (updated on June 2022)"                                  << std::endl; 
        std::cout << "* Written by Dr. Eunsang Lee"                                               << std::endl; 
        std::cout << "* Theoretical Physical Chemistry Dept. TU Darmstadt "                       << std::endl; 
        std::cout << "* Email: e.lee@theo.chemie.tu-darmstadt.de"                                 << std::endl; 
        std::cout << "***************************************************************************"<< std::endl<<std::endl;
    }
    return;
}


void Command::setAnalysisProgram(){
    /*Assining the name of program
     *If you want to add your program, assign the global varialbe for the program in also types.hpp*/
    analysis_program["run"]=RUN;
    analysis_program["dropsize"]=DROPSIZE;
    analysis_program["velocity"]=VELOCITY;
    analysis_program["sphstress"]=SPHERICALSTRESS;
    analysis_program["rdensity"]=RADIALDENSITY;
    analysis_program["avgstress"]=AVGSTRESS;
    analysis_program["polads"]=POLADSORP;
    analysis_program["polsize"]=POLSIZE;
    analysis_program["polevrlx"]=POLEVRLX;
    analysis_program["polorient"]=POLORIENT;
    analysis_program["bondlen"]=BONDLENGTH;
    analysis_program["polstretch"]=POLSTRETCH;
    analysis_program["trjtogro"]=TRJTOGRO;
    analysis_program["polsmsf"]=POLSMSF;
    analysis_program["polsubsize"]=POLSUBSIZE;
    analysis_program["msd"]=POLMSD;
    analysis_program["velacf"]=VELACF;
    analysis_program["surfcov"]=SURFCOV;
    analysis_program["rdf"]=RDF;
    analysis_program["depori"]=DEPORIENT;
    analysis_program["brdgsize"]=BRIDGESIZE;
    analysis_program["brdgvel"]=BRIDGEVEL;
    analysis_program["brdgpc"]=BRIDGEPCONC;
    analysis_program["brdgac"]=BRIDGEADSCONC;
    analysis_program["brdgclvel"]=BRIDGECLVEL;
    analysis_program["dropzd"]=DROPZD;
    analysis_program["surfrdf"]=DROPSRDF;
    analysis_program["jumpfreq"]=JUMPFREQ;
    analysis_program["brdgjpf"]=BRIDGEJF;
    analysis_program["brdgslvel"]=BRIDGESLVEL;
    analysis_program["brdgvelx"]=BRIDGEVELX;
    analysis_program["brdgvelxz"]=BRIDGEVELXZ;
    analysis_program["brdgzd"]=BRIDGEZD;
    analysis_program["brdgxzd"]=BRIDGEDENSXZ;
    analysis_program["brdgclin"]=BRIDGECLINE;
    analysis_program["brdginterf"]=BRIDGEINTERF;
    analysis_program["brdggpd"]=BRIDGEGPDEN;
    analysis_program["nliqcl"]=NUMCLUSTER;
    analysis_program["newprop"]=NEWPROPERTY;
    return;
}






void Command::getProgramOption(){
    bool error=true;
        
    if(mpi->rank()==MASTER){
        if(argc==1){
            std::cout << "[Error] Please give a program option." << std::endl;
            error=true;
            MPI_Abort(MPI_COMM_WORLD, 0);
            exit(0);
        }
        program=analysis_program[argv[1]];
        if(program==0){
            std::cout << "[Error] Program " << argv[1] << " can not be understood." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            exit(0);
        }
        else if(program==RUN)
            getCommandOption();
        else
            getCommandOptionAnalysis();

    }
    MPI_Bcast(&program, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    return;
}



/*Getting options for analysis programs*/
void Command::getCommandOptionAnalysis(){
    bool given_traj=false;
    bool given_strs=false;
    bool given_frc=false;
    if(mpi->rank()==MASTER){

        /*Checking whether unknown option flag is used*/
        /*
        Svec keys={"-c", "-o", "-t", "-p", "-i", "-n", "-x", "-f", "-s"};
        for(int i=2; i<argc; i+=2){
            if(std::find(keys.begin(), keys.end(), argv[i])==keys.end()){
                std::cout << "[Error] The command line option " << argv[i] << " is not understood." << std::endl << std::flush;
                MPI_Abort(MPI_COMM_WORLD, 11);
                exit(0);
            }
        }
        */

       
        /*Reading output prefix*/
        getCommandSingleOption("-i", "default", &out_prefix);
        if(out_prefix.compare("default")){
            traj_fname=out_prefix+".trj";
            log_fname=out_prefix+".log";
            strs_fname=out_prefix+".str";
            frc_fname=out_prefix+".frc";
        }


        /*Reading specific file names*/
        getCommandSingleOption("-c", "conf.gro", &config_fname);
        getCommandSingleOption("-o", "default", &out_fname);
        getCommandSingleOption("-t", "topol.top", &topol_fname);
        getCommandSingleOption("-p", "control.in", &control_fname);
        getCommandSingleOption("-n", "index.ndx", &ndx_fname);
        getCommandSingleOption("-x", "traj.trj", &traj_fname);
        getCommandSingleOption("-f", "force.frc", &frc_fname);
        getCommandSingleOption("-s", "stress.str", &strs_fname);




        std::cout << " Initial configuration from "<< config_fname << std::endl;
        std::cout << " Control parameters from "<< control_fname << std::endl;
        std::cout << " Force field from "<< topol_fname << std::endl;
        std::cout << " Trajectory from "<< traj_fname << std::endl  << std::endl << std::endl;
    }


}


void Command::getCommandOption(){
    do_restart=false;
    if(mpi->rank()==MASTER){
        
        /*Checking whether unknown option flag is used*/
        Svec keys={"-c", "-o", "-t", "-p", "-l", "-x", "-f", "-s", "-r", "-restart", "-ss", "-sg"};
        for(int i=2; i<argc; i+=2){
            if(std::find(keys.begin(), keys.end(), argv[i])==keys.end()){
                std::cout << "[Error] The command line option " << argv[i] << " is not understood." << std::endl << std::flush;
                MPI_Abort(MPI_COMM_WORLD, 11);
                exit(0);
            }
        }

        /*Reading output prefix*/
        getCommandSingleOption("-o", "default", &out_prefix);
        
        if(out_prefix.compare("default")==0){   //As default
            traj_fname="traj.trj";
            strs_fname="stress.str";
            frc_fname="force.frc";
            ckp_fname="check.ckp";
            log_fname="log.out";
            ss_fname="sslog.out";
            ssgro_fname="sstrj.gro";
            fnconf_fname="final.gro";
        }
        else{               //Otherwise
            traj_fname=out_prefix+".trj";
            log_fname=out_prefix+".log";
            strs_fname=out_prefix+".str";
            frc_fname=out_prefix+".frc";
            fnconf_fname=out_prefix+"_final.gro";
            ckp_fname=out_prefix+".ckp";
            ss_fname=out_prefix+".sslog";
            ssgro_fname=out_prefix+"_sstrj.gro";
        }

        /*Reading specific file names*/

        getCommandSingleOption("-c", "conf.gro", &config_fname);
        getCommandSingleOption("-t", "topol.top", &topol_fname);
        getCommandSingleOption("-p", "control.in", &control_fname);
        getCommandSingleOption("-l", "log.out", &log_fname);
        getCommandSingleOption("-x", "traj.trj", &traj_fname);
        getCommandSingleOption("-f", "force.frc", &frc_fname);
        getCommandSingleOption("-s", "stress.str", &strs_fname);
        getCommandSingleOption("-r", "check.ckp", &ckp_fname);
        getCommandSingleOption("-ss", "sslog.out", &ss_fname);
        getCommandSingleOption("-sg", "sstrj.gro", &ssgro_fname);
        do_restart=getCommandSingleOption("-restart", "restart.ckp", &rst_fname);

        if(do_restart){
            if(rst_fname.compare(ckp_fname)==0){
                ckp_fname=ckp_fname+".new";
            }
        }
        std::cout << " Initial configuration from "<< config_fname << std::endl;
        std::cout << " Control parameters from "<< control_fname << std::endl;
        std::cout << " Force field from "<< topol_fname << std::endl << std::endl;
    }


}


bool Command::getCommandSingleOption(const char opt[], std::string deflt, std::string *value){
    if(mpi->rank()==MASTER){
        *value=deflt;
        for(int i=2; i<argc; i+=2){
            if(argv[i].compare(opt)==0){
                *value=argv[i+1];
                return true;
            }
        }
    }
    return false;
}
bool Command::getCommandSingleOption(const char opt[], int deflt, int *value){
    if(mpi->rank()==MASTER){
        *value=deflt;
        for(int i=2; i<argc; i+=2){
            if(argv[i].compare(opt)==0){
                *value=std::stoi(argv[i+1]);
                return true;
            }
            
        }
    }
    return false;
}
bool Command::getCommandSingleOption(const char opt[], real deflt, real *value){
    if(mpi->rank()==MASTER){
        *value=deflt;
        for(int i=2; i<argc; i+=2){
            if(argv[i].compare(opt)==0){
                *value=std::stod(argv[i+1]);
                return true;
            }
        }
    }
    return false;
}

Command::~Command(){
}

void Command::assignDefaults(){
    config_fname="conf.gro";
    topol_fname="topol.top";
    control_fname="control.in";
    ss_fname="sslog.out";
    ndx_fname="index.ndx";
    out_prefix="default";
    first_step=0;
    final_step=-1;
    skip_step=1;
}
