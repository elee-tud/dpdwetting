#ifndef __SLLODDPD__HPP
#define __SLLODDPD__HPP

#include "../integrator.hpp"

namespace dpd{

class SllodDPD : public Integrator{
private:
    Rvec shearten;
    Real3D halfbox;
    Real3D refpos;

public:
    SllodDPD(){}
    SllodDPD(Initialization* init);
    ~SllodDPD(){}

    void updatePosition();
    void updatePosition(bool save_prev);
    void updateVelocity();
    void updateVelocity(bool save_prev);

    real getMaxForce(){ return 0.;}
};
};

#endif
