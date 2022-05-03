#ifndef __DUMP__HPP
#define __DUMP__HPP

#include "setmpi.hpp"
#include "initialization.hpp"
#include "mpiclasses.hpp"
#include "analysis.hpp"
#include <fstream>

namespace dpd{


class Dump{
private:
    Initialization* init;
    Command* command;
    SetMPI* mpi;
    Analysis* analysis;
    Configuration* config;
    Decomposition* decomp;
    Control* control;
    ParticleList particles;
    Real3D box;
    MPIClasses mpiptcl;
   

    bool restart;
    std::string traj_fname;
    std::string log_fname;
    std::string stress_fname;
    std::string force_fname;

    std::ofstream trajstream;
    std::ofstream logstream;
    std::ofstream stressstream;
    std::ofstream forcestream;

    bool dumpbinary;
    int dumpmult;
    bool dumpfrozen;

    int nufbeads;

    void writeLogHeader();

    

public:
    Dump(){}
    Dump(Initialization* init, Analysis* analysis);
    ~Dump();

    void dumpConf(bool on_screen, int step, real t);
    void dumpConf(int step, real t);
    void dumpConf(int step, real t, std::string filename);
    void dumpForce(int step, real t);
    void dumpDensity();
    void dumpLocalDensity();
    void dumpLog(int step, real time, std::string remaint);
    void dumpStress(int step, real time);
};

};

#endif
