/*************************************************************
 *            Class Decomposition
 *
 *  This class is to set up the domain decomposition of the 
 * system. 
 ************************************************************/

#ifndef __DECOMPOSITION__HPP
#define __DECOMPOSITION__HPP

#include "types.hpp"
#include "control.hpp"
#include "configuration.hpp"
#include "errorhandling.hpp"
#include "periodicboundary.hpp"
#include "indexing.hpp"
#include "cell.hpp"
#include "mpiclasses.hpp"

namespace dpd{

class Decomposition{
private:

    SetMPI* mpi;
    Control* control;
    Topology* topology;
    Configuration* config;


    PeriodicBoundary pbc;       //Periodic boundary 
    ParticleList particles;     //List of particle pointers
    Error err;                  //Error handling
    Indexing indexing_domain;   //Indexing class for domain decomposition
    Indexing indexing_cell;     //Indexing class for celllist
    MPIClasses mpiptcl;         //MPI_Particle class

    Real3D box;                 //Box dimension
    int numprocs;               //The number of total processors
    Ivec num_domains;           //The number of domains in each direction

    int myid;                   //Processor index
    Ivec my3did;                //3D process index
    real saferatio;             //Ratio to reserve the memory for particles wrt. average number of particles in a domain
    int maxnbeads;              //Maximum num. of beads in each processor

    Real3D domain_length;       //The length of each domain
    Real3D domain_min;          //The minimum position of domain
    Real3D domain_max;          //The maximum position of domain
    int num_mybeads;            //The number of beads in this processor
    Ivec mybeads;               //Bead indices in this processor
    
    Ivec2D ghosts;              //Destination(processorid, cellid) of ghost beads
    Ivec2D realbeads;           //Destination(processorid, cellid) of real beads

    Ivec truecells;           //List of true cell indices
    Ivec ghostcells;            //List of ghost cell indices
    Ivec realcells;             //List of real cell indices

    Ivec neighbor_proc;         //List of neighbor processor indices
    Ivec2D nbsearch;            //List of neighbor directions

    
    CellList cells;             //List of cell pointers

    real rcut;                  //Cutoff distance for a cell
    Real3D cell_length;         //Cell length in each direction
    Ivec num_true_cells;        //The 3D number of true cells
    Ivec num_cells;             //The 3D number of all cells in this processor
    int totnum_cells;           //The number of all cells
    int totnum_ghost_cells;     //The number of ghost cells
    int totnum_true_cells;      //The number of true cells
    int prev_totnum_cells;      //The number of true cells

    bool is_coord_synced;
    bool is_veloc_synced;
    bool is_force_synced;
    bool is_density_synced;
    bool is_stress_synced;

    bool dumpfrozen;


    void calculateDomainDivisor();
    void calculateDomainLength();
    bool calculateCellLength();
    void makeCellLists();
    Ivec getUnfrozenBeads();

public:
    Decomposition(){}
    Decomposition(Control* control, Configuration* config, SetMPI* mpi);
    Decomposition(Control* control, Configuration* config, SetMPI* mpi, int numprocs, real rcut);
    Decomposition(Control* control, Configuration* config, SetMPI* mpi, bool decomp_anal);
    ~Decomposition(){}

    Control* getControl(){ return control; }
    Configuration* getConfiguration(){ return config; }
    Topology* getTopology(){ return topology; }
    SetMPI* getMPI(){ return mpi;}
    ParticleList& getParticles(){ return particles; }
    void clearAllBeadInformation();
    void allocateBeadsToDomain();
    void allocateBeadsToCells();
    void assignParticleType();

    Ivec& getNumberDomains(){ return num_domains; }
    int getNumberOfBeadsInDomain(){ return num_mybeads; }
    Ivec& getBeadsIndexInDomain(){ return mybeads; }
    Ivec& getNeighborDomains(){ return neighbor_proc;}
    Real3D& getDomainLength(){ return domain_length;}
    Real3D& getDomainMinimumPosition(){ return domain_min; }
    Real3D& getDomainMaximumPosition(){ return domain_max; }
    int& getMyProcessorID(){ return myid; }
    Ivec& getMy3DProcessorID(){ return my3did; }
    int getNeighborProcessorID(Ivec dist);
    void addBeads(int index);
    void addBeads(Ivec indices);
    void removeBeads(int index);
    void removeBeads(Ivec indices);
    void printTrueBeadsInDomain();
    void refreshDomainBeads();
    void communicateGhostBeads();
    int getNewCellIndex(Particle* ptcl);
    int getNewDomainIndex(Particle* ptcl);
    Ivec getNewDomainCellIndex(Particle* ptcl);
    void matchGhostsRealBeads();
//    Ivec2D findProcsOfGhosts(){ return ghosts; }
//    Ivec2D findProcsOfRealBeads(){ return realbeads; }
    Ivec2D getProcsOfGhosts(){ return ghosts; }
    Ivec2D getProcsOfRealBeads(){ return realbeads; }


    Ivec& getNumberOfCells(){ return num_cells; }
    Ivec& getNumberOfTrueCells(){ return num_true_cells; }
    Real3D& getCellLength(){ return cell_length; }
    int& getTotalNumberOfCells(){ return totnum_cells;}
    int& getTotalNumberOfGhostCells(){ return totnum_ghost_cells;}
    int& getTotalNumberOfTrueCells(){ return totnum_true_cells;}
    CellList& getCells(){ return cells;}
    Ivec& getTrueCells(){ return truecells; }
    Ivec& getGhostCells(){ return ghostcells;}
    Ivec& getRealCells(){ return realcells; }
    void printCellInformation();
    void setUnsynced();
    void gatherCoords();
    void gatherVelocs();
    void gatherForces();
    void gatherDensities();
    void gatherStresses();
    void gatherAll();
    inline bool isCoordSynced(){ return is_coord_synced; }
    inline bool isVelocSynced(){ return is_veloc_synced; }
    inline bool isForceSynced(){ return is_force_synced; }
    inline bool isDensitySynced(){ return is_density_synced; }
    inline bool isStressSynced(){ return is_stress_synced; }
   
    Real3D getBox(){ return box; }
    void resetBox(Real3D newbox);
    void whereIsTheBead(int index);
    void refreshBeadsWhenChangingCellNumber();
    
    void printDomainLengthInformation();
    void printCellNumberInformation();

    void setExistence();
};
};
#endif
