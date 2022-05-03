#ifndef __CALCULATION__HPP
#define __CALCULATION__HPP
#include "initanalysis.hpp"
namespace dpd{
class Calculation{
private:
    InitAnalysis* init;
    Property* prop;
    Trajectory* traj;
    Control *control;
    SetMPI* mpi;
    Timer* timer;

    int nbeads;
    int nsteps;
    int begstep;
    int endstep;
    int skipstep;


public:
    Calculation(){}
    Calculation(InitAnalysis* init);
    ~Calculation(){}
    void calculateProperty();

    


};

};

#endif
