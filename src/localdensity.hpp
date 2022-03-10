#ifndef __LOCALDENSITY__HPP
#define __LOCALDENSITY__HPP
#include "initialization.hpp"

namespace dpd{

class LocalDensity{
private:
    Initialization* init;
    Decomposition* decomp;
    Topology* topol;
    SetMPI* mpi;
    Configuration* config;
    Rvec2D Brcut;
    ParticleList particles;
    CellList cells;
    PeriodicBoundary pbc;
    Real3D box;
    Ivec truecells;
    Ivec ghostcells;

    real normconst;
    void setZeroLocalDensity();
    void calculateDensityInTheCells(Cell* cell1);
    void calculateDensityBetweenCells(Cell* cell1, Cell* cell2);
    void reduceLocalDensity();
    void redistLocalDensity();
public:
    LocalDensity(){}
    LocalDensity(Initialization* init);
    LocalDensity(Decomposition* decomp);
    ~LocalDensity(){}


    void calculateLocalDensity();
    void calculatePairDensity(int index1, int index2);

};

};
#endif
