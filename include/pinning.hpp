#ifndef __PINNING__HPP
#define __PINNING__HPP

#include "../extension.hpp"

namespace dpd{

class Pinning:public Extension{
private:
    Ivec direct;
    Ivec pinindex;
    Particle2DList molecules;
    std::string molname;
    int npinpermol;


public:
    Pinning(){}
    Pinning(Topology* topol, Configuration* config, Decomposition* decomp);
    ~Pinning(){}

    void applyExtensionForPosition(int step);
    void applyExtensionForVelocity(int step);
};
};
#endif




