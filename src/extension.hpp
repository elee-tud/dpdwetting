#ifndef __EXTENSION__HPP
#define __EXTENSION__HPP

#include "decomposition.hpp"

namespace dpd{

class Extension{
protected:
    Topology* topol;
    Configuration* config;
    Decomposition* decomp;
    Control* control;
    SetMPI* mpi;
    bool need_prev_position;
    bool need_prev_velocity;
    bool for_position;
    bool for_velocity;
    CellList cells;
    ParticleList particles;
public:
    Extension(){}
    Extension(Topology* topol, Configuration* config, Decomposition* decomp);
    ~Extension(){}

    virtual void applyExtensionForPosition(int step){}
    virtual void applyExtensionForVelocity(int step){}

    bool needPrevPosition(){ return need_prev_position; }
    bool needPrevVelocity(){ return need_prev_velocity; }
    bool isForPosition(){ return for_position; }
    bool isForVelocity(){ return for_velocity; }

};

typedef std::vector<Extension*> ExtensionList;
};
#endif



