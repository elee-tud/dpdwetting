#ifndef __ERRORHANDLING__HPP
#define __ERRORHANDLING__HPP
#include "setmpi.hpp"
#include "particle.hpp"

namespace dpd{

class Error{
private:
    SetMPI* mpi;

public:
    Error(){}
    Error(SetMPI* mpi):mpi(mpi){}
    ~Error(){}

    void mismatchedNumberOfParticles(int num_beads_conf, int nbeads);
    void wrongNumberOfCores();
    void mismatchMoleculeName(std::string lineingro, Particle* ptcl);
    void missingBondedParticle(int index1, int index2);
    void missingBonds(int calc_nbonds, int tot_nbonds);
    
};

};

#endif
