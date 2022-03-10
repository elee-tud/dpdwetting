#ifndef __NONBONDED__HPP
#define __NONBONDED__HPP

#include "decomposition.hpp"
#include "control.hpp"
#include "stress.hpp"

namespace dpd{

class Nonbonded{
protected:
    Topology* topol;
    Configuration* config;
    Decomposition* decomp;
    Real3D box;
    PeriodicBoundary pbc;
    ParticleList particles;
    CellList cells;
    Ivec truecells;
    Ivec ghostcells;
    Control* control;
    SetMPI* mpi;
    int integrator;


public:
    Nonbonded(){}
    Nonbonded(Topology* topol, Configuration* config, Decomposition* decomp);
    ~Nonbonded(){}

    void calculateForce();
    void calculateForceInTheCells(Cell* cell);
    void calculateForceBetweenCells(Cell* cell1, Cell* cell2);
    virtual void calculatePairForce(int index1, int index2)=0;
    
    void calculateForce(Real3D com, real** press, real dr);
    void calculateForceInTheCells(Cell* cell, Real3D com, real** press, real dr);
    void calculateForceBetweenCells(Cell* cell1, Cell* cell2, Real3D com, real** press, real dr);
    virtual void calculatePairForce(int index1, int index2, Real3D com, real** press, real dr)=0;



};

typedef std::vector<Nonbonded*> NonbondedList;
};

#endif
