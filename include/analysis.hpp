#ifndef __ANALYSIS__HPP
#define __ANALYSIS__HPP

#include "initialization.hpp"

namespace dpd{

class Analysis{
private:
    Initialization* init;
    Decomposition* decomp;
    SetMPI* mpi;
    ParticleList particles;
    Ivec mybeads;
    Real3D box;
    real boxvolume;
    
    int ntot_beads;
    int ntot_unfrozen_beads;
    real temperature;
    real pressure;
    Rvec virialt;

public:
    Analysis(){}
    Analysis(Initialization* init);
    ~Analysis(){}
   
    void calculateSystemProperties();
    void calculateTemperature();
    void calculateStressTensor();

    real getTemperature(){ return temperature; }
    real getPressure(){ return pressure; }
    Rvec getStressTensor(){ return virialt;}

    
};
};
#endif
