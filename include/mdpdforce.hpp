#ifndef __MDPDFORCE__HPP
#define __MDPDFORCE__HPP

#include "../nonbonded.hpp"
#include <random>

namespace dpd{

class MDPDForce:public Nonbonded{
protected:
    Rvec2D A;
    Rvec2D Arcut;
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
    MDPDForce(){}
    MDPDForce(Topology* topol, Configuration* config, Decomposition* decomp);
    ~MDPDForce(){}
    virtual void calculatePairForce(int index1, int index2);
    virtual void calculatePairForce(int index1, int index2, Real3D com, real** press, real dr);
};
};

#endif
