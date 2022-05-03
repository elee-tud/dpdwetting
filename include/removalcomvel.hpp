#ifndef __REMOVALCOMVEL__HPP
#define __REMOVALCOMVEL__HPP

#include "../extension.hpp"
namespace dpd{

class RemovalCOMVel:public Extension{
private:
    int rmcomv_freq;
    Ivec rmcomv_dir;
    Real3D comvel;
    Real3D comvelproc;

public:
    RemovalCOMVel(){}
    RemovalCOMVel(Topology* topol, Configuration* config, Decomposition* decomp);
    ~RemovalCOMVel(){}

    void applyExtensionForPosition(int step){ return; }
    void applyExtensionForVelocity(int step);
};
};
#endif
