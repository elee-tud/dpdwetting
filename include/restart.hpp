#ifndef __RESTART__HPP
#define __RESTART__HPP
#include "command.hpp"
#include "configuration.hpp"
#include "topology.hpp"
#include "decomposition.hpp"
#include "setmpi.hpp"
#include <iostream>
#include <fstream>

namespace dpd{

class Restart{
protected:
    Command* command;
    Configuration* config;
    Topology* topol;
    Decomposition *decomp;
    SetMPI* mpi;
    ParticleList particles;
    Real3D box;
    int dumpmult;


    std::string ckp_fname;
    std::ofstream ckpstream;

    std::string rst_fname;
    std::ifstream rststream;

public:
    Restart(){}
    Restart(Command* command, Configuration* config, Topology* topol, Decomposition* decomp, SetMPI* mpi);
    ~Restart(){}


    void writeCheckPoint(int step, real t);
    void readCheckPoint(int *step, real *t);
};
};

#endif
