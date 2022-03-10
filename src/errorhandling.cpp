#include "errorhandling.hpp"
#include "parsing.hpp"
#include <iostream>

using namespace dpd;
void Error::mismatchedNumberOfParticles(int num_beads_conf, int num_beads){
    bool error=false;
    if(num_beads_conf!=num_beads){
        std::cout << "[Error] The number of beads in the configuration file, " << num_beads_conf << ", is different from defined in the topology file, " << num_beads << "." << std::endl;
        error=true;
    }
    if(error){
        mpi->finalize();
        exit(0);
    }
    return;
}

void Error::wrongNumberOfCores(){
    bool error=false;
    if(mpi->isMaster()){
        std::cout << "[Error] The number of cores should be given among 1, 2, 4, 6, 8, 10, 12, 16, 18, 24, 32, 48, 64." << std::endl;
        error=true;
    }
    if(error){
        mpi->finalize();
        exit(0);
    }
    return;
}

void Error::mismatchMoleculeName(std::string lineingro, Particle* ptcl){
    std::string name=lineingro.substr(5,5);
    if(dpd::trim_copy(lineingro.substr(5,5))!=dpd::trim_copy(ptcl->getMoleculeName())){
        std::cout << "[Error] The molecule names the configuration file and the topology file are not matched." << std::endl;
        mpi->finalize();
        exit(0);
    }

    return;
}

void Error::missingBondedParticle(int index1, int index2){
    bool error=false;
    std::cout << "[Error] Particle No. " << index2 << " is missing for the bonded force with Particle No. " << index1 << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(0);
    return;
}

void Error::missingBonds(int calc_nbonds, int tot_nbonds){
    bool error=false;
    if(mpi->isMaster()){
        if(calc_nbonds!=tot_nbonds){
            std::cout << "[Error] "<< tot_nbonds-calc_nbonds <<" particles are missing to calculate bonded parameter. Use a longer cuttof distance or a shorter time step.\n"<< std::endl ;
            error=true;
        }
    }
    if(error){
        mpi->finalize();
        exit(0);
    }
}
