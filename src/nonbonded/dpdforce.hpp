#ifndef __DPDFORCE__HPP
#define __DPDFORCE__HPP

#include "../nonbonded.hpp"
#include <random>

namespace dpd{

class DPDForce:public Nonbonded{
protected:
    Rvec2D B;
    Rvec2D Brcut;
    real gamma;
    real temp;
    real q;
    real invsqrtdt;
    int randseed;
    std::normal_distribution<real> normdist;
    std::default_random_engine generator;




public:
    DPDForce(){}
    DPDForce(Topology* topol, Configuration* config, Decomposition* decomp);
    ~DPDForce(){}
    virtual void calculatePairForce(int index1, int index2);
    virtual void calculatePairForce(int index1, int index2, Real3D com, real** press, real dr);
};
};

#endif
