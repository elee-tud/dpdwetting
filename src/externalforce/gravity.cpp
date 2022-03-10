#include "gravity.hpp"

using namespace dpd;
Gravity::Gravity(Topology* topol, Configuration* config, Decomposition* decomp):ExternalForce(topol, config, decomp){
    field=control->getGravityField();
    direct=control->getGravityDirection();

}

void Gravity::calculateSingleForce(Particle* ptcl){
    ptcl->force[direct]+=-field/ptcl->getMass();
    return;
}

