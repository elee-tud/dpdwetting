#ifndef __ENERMININTEGRATOR__HPP
#define __ENERMININTEGRATOR__HPP

#include "../integrator.hpp"

namespace dpd{

class EnerMinIntegrator: public Integrator{
private:
    real maxforce;
    real maxdisp;
    real startdisp;
    real incratio;
    void findMaxForce();

public:
    EnerMinIntegrator(){}
    EnerMinIntegrator(Initialization* init);
    ~EnerMinIntegrator(){}

    void updatePosition();
    void updatePosition(bool save_prev);
    void updateVelocity(){}
    void updateVelocity(bool save_prev){}

    real getMaxForce(){ return maxforce; }
};
};


#endif
