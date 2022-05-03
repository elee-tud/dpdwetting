#ifndef __STERSS__HPP
#define __STERSS__HPP
#include "particle.hpp"
#include "types.hpp"

namespace dpd{

void calculateStress(Particle* ptcl1, Particle* ptcl2, Real3D rij, Real3D fij);
void calculateSphericalStress(Particle* ptcl1, Particle* ptcl2, Real3D rij, Real3D fij, Real3D com, real** press, real dr);
};

#endif
