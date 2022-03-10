#ifndef __UNIFLOWPBC__HPP
#define __UNIFLOWPBC__HPP


#include "../extension.hpp"
namespace dpd{

class UniflowPBC:public Extension{
private:
    int freq;
    real factor;
    real invsqrtf;
    Real3D oldbox;
    Rvec sheartensor;
    real dt;


public:
    UniflowPBC(){}
    UniflowPBC(Topology* topol, Configuration* config, Decomposition* decomp);
    ~UniflowPBC(){}

    virtual void applyExtensionForPosition(int step);


};
};

#endif
