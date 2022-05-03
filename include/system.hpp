#ifndef __SYSTEM__HPP
#define __SYSTEM__HPP

#include "initialization.hpp"
#include "dump.hpp"
#include "force.hpp"
#include "localdensity.hpp"
#include "integrators/velocityverletdpd.hpp"
#include "integrators/enerminintegrator.hpp"
#include "integrators/slloddpd.hpp"
#include "analysis.hpp"

namespace dpd{

class System{
private:
    SetMPI* mpi;
    Initialization *init;
    Command* command;
    Topology* topology;
    Control* control;
    Configuration* config;
    Interactions* interaction;
    Decomposition* decomp;
    CellList cells;
    ExtensionList extensionlist;
    ParticleList particles;
    ExtensionList extensions;
    Analysis* analysis;
    Dump* dump;
    Force* force;
    LocalDensity* localdensity;
    Integrator* integrator;
    Timer* timer;
    Restart* restart;
    SlipSpring* sspring;

    int startstep;
    real starttime;
    int trajfreq;
    int logfreq;
    int strsfreq;
    int frcfreq;
    
    bool need_prev_position;
    bool need_prev_velocity;
    void setIntegrator();    
    void resolveDependencyForExtensions();

public:
    System(){}
    System(Initialization *init);
    ~System();

    void initializeSimulation();
    void run();
    Dump* getDump(){ return dump; }

};
};

#endif



