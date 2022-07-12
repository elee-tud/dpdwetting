#ifndef __VELOCITYVERLETDPD__HPP
#define __VELOCITYVERLETDPD__HPP

#include "../integrator.hpp"

namespace dpd{

class VelocityVerletDPD : public Integrator{
private:

public:
    VelocityVerletDPD(){}
    VelocityVerletDPD(Initialization* init);
    ~VelocityVerletDPD(){}

    void updatePosition();
    void updatePosition(bool save_prev);
    void updateVelocity();
    void updateVelocity(bool save_prev);

    real getMaxForce(){ return 0.;}
};
};

#endif
