#include "springpulling.hpp"

using namespace dpd;
SpringPulling::SpringPulling(Topology* topol, Configuration* config, Decomposition* decomp):ExternalForce(topol, config, decomp){
    springk=control->getPullingSpringK();
    center=control->getPullingCoord();
    direct=control->getPullingDirect();
    dir.reserve(3);
    if(springk!=0.0 && center==Real3D(-256.0))
            center=box/2;
    if(direct==0)
        dir=Ivec{0,1,2};
    else if(direct==1)
        dir=Ivec{0};
    else if(direct==2)
        dir=Ivec{1};
    else if(direct==3)
        dir=Ivec{2};
    else if(direct==4)
        dir=Ivec{0, 1};
    else if(direct==5)
        dir=Ivec{0, 2};
    else if(direct==6)
        dir=Ivec{1, 2};
    numdir=dir.size();
}

void SpringPulling::calculateSingleForce(Particle* ptcl){
    Real3D rij=pbc.getMinimumImageVector(ptcl->coord, center);
    for(int i=0;i<numdir;i++)
        ptcl->force[dir[i]]+=-springk*rij[dir[i]];
    return;
}

