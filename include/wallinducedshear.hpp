#ifndef __WALLINDUCEDSHEAR__HPP
#define __WALLINDUCEDSHEAR__HPP


#include "../extension.hpp"
#include "../periodicboundary.hpp"
namespace dpd{

class WallInducedShear:public Extension{
private:
    int sheardir;
    real shearrate;
    real shearvel;
    std::string group1;
    std::string group2;
    Ivec group1idx;
    Ivec group2idx;
    int group1size;
    int group2size;
    real botwallpos;
    real topwallpos;

    real height;
    int vdir;
    int hdir;
    real dt;
    real displ;

    PeriodicBoundary pbc;



public:
    WallInducedShear(){}
    WallInducedShear(Topology* topol, Configuration* config, Decomposition* decomp);
    ~WallInducedShear(){}

    void findGroupIndex();
    void assignVelocity();

    virtual void applyExtensionForPosition(int step);
};
};

#endif
