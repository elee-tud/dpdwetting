#ifndef __FORCE__HPP
#define __FORCE__HPP
#include "initialization.hpp"

namespace dpd{

class Force{
protected:
    Initialization* init;
    Decomposition* decomp;
    Interactions* inter;
    NonbondedList nonbondedlist;
    BondedList bondedlist;
    ExternalList externallist;
    CellList cells;
    ParticleList particles;
    SetMPI* mpi;
    MPIClasses mpiptcl;

public:
    Force(){}
    Force(Initialization* init);
    Force(Decomposition* decomp, Interactions* inter);
    ~Force(){}
    
    
    void setZeroForce();
    void calculateForce();
    void calculateForce(Real3D com, real** press, real dr);

    void reduceForce();
};




};


#endif
