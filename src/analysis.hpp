#ifndef __ANALYSIS__HPP
#define __ANALYSIS__HPP

/****************************************************************************
 * Class Analysis
 *
 * This class is used to anlayze basic properties of system like temperature
 * and pressure on the fly.
 ****************************************************************************/
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
  
    void calculateSystemProperties();   //Function to calculate temperature and pressure
    void calculateTemperature();        //Function to calculate temperature
    void calculateStressTensor();       //Function to calculate stress

    real getTemperature(){ return temperature; }    //Returning temperature value
    real getPressure(){ return pressure; }          //Returning pressure value
    Rvec getStressTensor(){ return virialt;}        //Returning pressure tensor

    
};
};
#endif
