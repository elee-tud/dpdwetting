#ifndef __TRAJECTORY__HPP
#define __TRAJECTORY__HPP

#include "../initialization.hpp"
#include "../errorhandling.hpp"
namespace dpd{

class Trajectory{

private:
    InitialSet initset;
    int type;
    Command *command;
    Configuration *config;
    Control *control;
    Topology *topol;
    Decomposition* decomp;
    SetMPI* mpi;
    std::string traj_fname;
    std::ifstream trajstream;
    ParticleList particles;
    Error err;
    std::string title_traj;
    Real3D box;
    int nbeads;

    bool dumpbinary;
    int dumpmult;

public:
    Trajectory(){}
    Trajectory(InitialSet initset, int type);
    ~Trajectory(){}

    void readFirstStep();
    void openTrajectory();
    bool readTrajectoryStep();
    bool skipTrajectoryStep();
    void closeTrajectory();
    Real3D getBox(){return box;}

};
};

#endif
