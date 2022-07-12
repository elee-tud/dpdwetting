#ifndef __SLIPSPRING__HPP
#define __SLIPSPRING__HPP

#include "control.hpp"
#include "topology.hpp"
#include "configuration.hpp"
#include "decomposition.hpp"

#define INTRASS 0
#define INTERSS 1

#define MAXSS 5

namespace dpd{


class SlipSpring{
private:
    SetMPI* mpi;
    Command* command;
    Control* control;
    Configuration* config;
    Decomposition* decomp;
    Decomposition* ssdecomp;
    Topology* topol;
    ParticleList particles;
    PeriodicBoundary pbc;

    Ivec endindex;
    int numends;

    int ntotptcls;
    Ivec mybeads;
    int totnumss;
    int nmcsteps;
    int ndpdsteps;
    real ssl0;
    real ssk;
    real ssmaxd;
    real ssmind;

    int numintrass;
    int sum_numintrass;


    int sspoltype;
    real intrassbias;
    int pollength;

    Real3D box; 
    Ivec2D allsprings;
    Ivec2D springs;
    int numss;
    CellList cells;


    real ssenergy_calc;
    real sum_ssenergy_calc;

    Ivec startidx;
    Ivec nssproc;
    real ssenergy;
    real deltae;
    real temp;
    real sum_ssenergy;
    unsigned long accepted;
    unsigned long sum_accepted;
    int relocated;
    int sum_relocated;
    std::string ssfname;
    std::string sstrjfname;
    std::ofstream ssstream;
    std::ofstream sstrjstream;
    bool restart;


public:
    SlipSpring(){}
    SlipSpring(Control* control, Topology* topol, Configuration* config, Decomposition* decomp);
    ~SlipSpring();

    void getPolymerLength();
    void initializeSlipSpring();
    int choosePairRandom(int index1);
    void saveCellIndexOfParticles();
    void bcastSprings();
    void calculateTotalEnergy();

    inline real springEnergy(real r){
        return ssk*(r-ssl0)*(r-ssl0)/2.;
    }
    void doMonteCarloSequence(int dpdstep);
    void doMonteCarloStep();
    void reduceSprings();
    void addSpringsToConfiguration();
    void removeAllSpringsFromConfiguration();
    void relocation();
    void relocationForLinearPolymer();
    void relocationForRingPolymer();
    void findPseudoEnds();
    void writeLogHeader();
    void dumpLog(int step);
    void dumpTrajectory(int step);
};
}
    

#endif
