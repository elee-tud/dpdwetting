#ifndef __BOXELONGATION__HPP
#define __BOXELONGATION__HPP


#include "../extension.hpp"
namespace dpd{

class BoxElongation:public Extension{
private:
    int freq;
    real factor;
    real invsqrtf;


public:
    BoxElongation(){}
    BoxElongation(Topology* topol, Configuration* config, Decomposition* decomp);
    ~BoxElongation(){}

    virtual void applyExtensionForPosition(int step);


};
};

#endif
