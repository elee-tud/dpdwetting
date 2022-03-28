#include "command.hpp"
#include <cstring>

using namespace dpd;
Command::Command(int argc, char* _argv[], SetMPI* _mpi):argc(argc), mpi(_mpi){
    argv=Svec(argc);
    for(int i=0;i<argc;i++){
        argv[i]=std::string(_argv[i]);
    }
    getProgramOption();
    printProgramStartingMessage();
   /*
   std::string test;
   getCommandSingleOption(argc, argv, std::string("-a"), std::string("default"), &test);
   std::cout << "TEST=" << test << std::endl;
   */


}

void Command::printProgramStartingMessage(){
    if(mpi->rank()==0){
        std::cout << "***************************************************************************"<< std::endl;

        if(program==RUN)
            std::cout << "*        Program for Dissipative Particle Dynamics Simulation              "<< std::endl;
        else
            std::cout << "*        Program for Analysing DPD trajectory                              "<< std::endl;
        std::cout << "***************************************************************************"<< std::endl;

        std::cout << "* Version: v3.0 (updated on July24  2021)"                                  << std::endl; 
        std::cout << "* Written by Dr. Eunsang Lee"                                               << std::endl; 
        std::cout << "* Theoretical Physical Chemistry Dept. TU Darmstadt "                       << std::endl; 
        std::cout << "* Email: e.lee@theo.chemie.tu-darmstadt.de"                                 << std::endl; 
        std::cout << "***************************************************************************"<< std::endl<<std::endl;
    }
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
        if(argv[1].compare("run")==0){
            program=RUN;
            getCommandOption();
            error=false;
        }
        else{
            if(argv[1].compare("dropsize")==0){
                program=DROPSIZE;
                error=false;
            }
            else if(argv[1].compare("velocity")==0){
                program=VELOCITY;
                error=false;
            }
            else if(argv[1].compare("sphstress")==0){
                program=SPHERICALSTRESS;
                error=false;
            }
            else if(argv[1].compare("rdensity")==0){
                program=RADIALDENSITY;
                error=false;
            }

            else if(argv[1].compare("avgstress")==0){
                program=AVGSTRESS;
                error=false;
            }
            else if(argv[1].compare("polads")==0){
                program=POLADSORP;
                error=false;
            }
            else if(argv[1].compare("polsize")==0){
                program=POLSIZE;
                error=false;
            }
            else if(argv[1].compare("polevrlx")==0){
                program=POLEVRLX;
                error=false;
            }
            else if(argv[1].compare("polorient")==0){
                program=POLORIENT;
                error=false;
            }
            else if(argv[1].compare("bondlen")==0){
                program=BONDLENGTH;
                error=false;
            }
            else if(argv[1].compare("polstretch")==0){
                program=POLSTRETCH;
                error=false;
            }
            else if(argv[1].compare("trjtogro")==0){
                program=TRJTOGRO;
                error=false;
            }
            else if(argv[1].compare("polsmsf")==0){
                program=POLSMSF;
                error=false;
            }
            else if(argv[1].compare("polsubsize")==0){
                program=POLSUBSIZE;
                error=false;
            }
            else if(argv[1].compare("msd")==0){
                program=POLMSD;
                error=false;
            }
            else if(argv[1].compare("velacf")==0){
                program=VELACF;
                error=false;
            }
            else if(argv[1].compare("surfcov")==0){
                program=SURFCOV;
                error=false;
            }
            else if(argv[1].compare("rdf")==0){
                program=RDF;
                error=false;
            }
            else if(argv[1].compare("depori")==0){
                program=DEPORIENT;
                error=false;
            }
            else if(argv[1].compare("brdgsize")==0){
                program=BRIDGESIZE;
                error=false;
            }
            else if(argv[1].compare("brdgvel")==0){
                program=BRIDGEVEL;
                error=false;
            }
            else if(argv[1].compare("brdgpc")==0){
                program=BRIDGEPCONC;
                error=false;
            }
            else if(argv[1].compare("brdgac")==0){
                program=BRIDGEADSCONC;
                error=false;
            }
            else if(argv[1].compare("brdgclvel")==0){
                program=BRIDGECLVEL;
                error=false;
            }
            else if(argv[1].compare("dropzd")==0){
                program=DROPZD;
                error=false;
            }
            else if(argv[1].compare("surfrdf")==0){
                program=DROPSRDF;
                error=false;
            }
            else if(argv[1].compare("jumpfreq")==0){
                program=JUMPFREQ;
                error=false;
            }
            else if(argv[1].compare("brdgjpf")==0){
                program=BRIDGEJF;
                error=false;
            }
            else if(argv[1].compare("brdgslvel")==0){
                program=BRIDGESLVEL;
                error=false;
            }
            else if(argv[1].compare("brdgvelx")==0){
                program=BRIDGEVELX;
                error=false;
            }
            else if(argv[1].compare("brdgvelxz")==0){
                program=BRIDGEVELXZ;
                error=false;
            }
            else if(argv[1].compare("brdgald")==0){
                program=BRIDGEALDENS;
                error=false;
            }
            else if(argv[1].compare("brdgzd")==0){
                program=BRIDGEZD;
                error=false;
            }
            else if(argv[1].compare("brdgxzd")==0){
                program=BRIDGEDENSXZ;
                error=false;
            }
            else if(argv[1].compare("brdgcline")==0){
                program=BRIDGECLINE;
                error=false;
            }
            getCommandOptionAnalysis();
        }
        if(error && mpi->isMaster()){
            std::cout << "program " << argv[1] << " can not be understood." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            exit(0);
        }


    }
    MPI_Bcast(&program, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    return;
}

void Command::getCommandOptionAnalysis(){
    bool given_traj=false;
    bool given_strs=false;
    bool given_frc=false;
    out_fname="default";
    if(mpi->rank()==MASTER){
        assignDefaults();
        for(int i=2; i<argc; i+=2){
            if(argv[i].compare("-c")==0)
                config_fname=argv[i+1];
            else if(argv[i].compare("-o")==0)
                out_fname=argv[i+1];
            else if(argv[i].compare("-t")==0)
                topol_fname=argv[i+1];
            else if(argv[i].compare("-p")==0)
                control_fname=argv[i+1];
            else if(argv[i].compare("-i")==0)
                out_prefix=argv[i+1];
            else if(argv[i].compare("-n")==0)
                ndx_fname=argv[i+1];
            else if(argv[i].compare("-x")==0){
                traj_fname=argv[i+1];
                given_traj=true;
            }
            else if(argv[i].compare("-f")==0){
                frc_fname=argv[i+1];
                given_frc=true;
            }
            else if(argv[i].compare("-s")==0){
                strs_fname=argv[i+1];
                given_strs=true;
            }
        }
        if(out_prefix.compare("default")==0){
            if(!given_traj)
                traj_fname="traj.trj";
            if(!given_strs)
                strs_fname="stress.frc";
            if(!given_frc)
                frc_fname="force.str";
        }
        else{
            traj_fname=out_prefix+".trj";
            log_fname=out_prefix+".log";
            strs_fname=out_prefix+".str";
            frc_fname=out_prefix+".frc";
        }
        std::cout << " Initial configuration from "<< config_fname << std::endl;
        std::cout << " Control parameters from "<< control_fname << std::endl;
        std::cout << " Force field from "<< topol_fname << std::endl;
        std::cout << " Trajectory from "<< traj_fname << std::endl  << std::endl << std::endl;
    }


}


void Command::getCommandOption(){
    bool given_traj=false;
    bool given_strs=false;
    bool given_frc=false;
    bool given_ckp=false;
    do_restart=false;
    if(mpi->rank()==MASTER){
        assignDefaults();
        for(int i=2; i<argc; i+=2){
            if(argv[i].compare("-c")==0)
                config_fname=argv[i+1];
            else if(argv[i].compare("-t")==0)
                topol_fname=argv[i+1];
            else if(argv[i].compare("-p")==0)
                control_fname=argv[i+1];
            else if(argv[i].compare("-o")==0)
                out_prefix=argv[i+1];
            else if(argv[i].compare("-l")==0)
                log_fname=argv[i+1];
            else if(argv[i].compare("-x")==0){
                traj_fname=argv[i+1];
                given_traj=true;
            }
            else if(argv[i].compare("-f")==0){
                frc_fname=argv[i+1];
                given_frc=true;
            }
            else if(argv[i].compare("-s")==0){
                strs_fname=argv[i+1];
                given_strs=true;
            }
            else if(argv[i].compare("-r")==0){
                ckp_fname=argv[i+1];
                given_ckp=true;
            }
            else if(argv[i].compare("-restart")==0){
                rst_fname=argv[i+1];
                do_restart=true;
            }
            else if(argv[i].compare("-ss")==0){
                ss_fname=argv[i+1];
            }
            else{
                std::cout << "[Error] The option " << argv[i] << " is not understood." << std::endl << std::flush;
                MPI_Abort(MPI_COMM_WORLD, 11);
                exit(0);
            }
        }
        if(out_prefix.compare("default")==0){
            if(!given_traj)
                traj_fname="traj.trj";
            if(!given_strs)
                strs_fname="stress.str";
            if(!given_frc)
                frc_fname="force.frc";
            if(!given_ckp)
                ckp_fname="check.ckp";
            log_fname="log.out";
            ss_fname="sslog.out";
            ssgro_fname="sstrj.gro";
            fnconf_fname="final.gro";
        }
        else{
            traj_fname=out_prefix+".trj";
            log_fname=out_prefix+".log";
            strs_fname=out_prefix+".str";
            frc_fname=out_prefix+".frc";
            fnconf_fname=out_prefix+"_final.gro";
            ckp_fname=out_prefix+".ckp";
            ss_fname=out_prefix+".sslog";
            ssgro_fname=out_prefix+"_sstrj.gro";
        }
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


void Command::getCommandSingleOption(const char opt[], std::string deflt, std::string *value){
    if(mpi->rank()==MASTER){
        *value=deflt;
        for(int i=2; i<argc; i+=2){
            if(argv[i].compare(opt)==0)
                *value=argv[i+1];
        }
    }
    return;
}
void Command::getCommandSingleOption(const char opt[], int deflt, int *value){
    if(mpi->rank()==MASTER){
        *value=deflt;
        for(int i=2; i<argc; i+=2){
            if(argv[i].compare(opt)==0)
                *value=std::stoi(argv[i+1]);
        }
    }
    return;
}
void Command::getCommandSingleOption(const char opt[], real deflt, real *value){
    if(mpi->rank()==MASTER){
        *value=deflt;
        for(int i=2; i<argc; i+=2){
            if(argv[i].compare(opt)==0)
                *value=std::stod(argv[i+1]);
        }
    }
    return;
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
