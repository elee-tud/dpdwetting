#ifndef __SOLIDWALL__HPP
#define __SOLIDWALL__HPP

#include "../extension.hpp"

namespace dpd{

class SolidWall:public Extension{
private:
    Ivec direct;
    Rvec dmin;
    Rvec dmax;

    bool isCrossingWall(Particle* ptcl); 
public:
    SolidWall(){}
    SolidWall(Topology* topol, Configuration* config, Decomposition* decomp);
    ~SolidWall(){}

    void applyExtensionForPosition(int step);
    void applyExtensionForVelocity(int step);
};
};
#endif




