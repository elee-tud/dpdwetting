#ifndef __INTERACTIONS__HPP
#define __INTERACTIONS__HPP

#include "control.hpp"
#include "topology.hpp"
#include "nonbonded/mdpdforce.hpp"
#include "nonbonded/dpdforce.hpp"
#include "bonded/bondlenharmonic.hpp"
#include "bonded/slipspringharmonic.hpp"
#include "externalforce/springpulling.hpp"
#include "externalforce/gravity.hpp"
#include "extensions/solidwall.hpp"
#include "extensions/removalcomvel.hpp"
#include "extensions/boxelongation.hpp"
#include "extensions/pinning.hpp"
#include "extensions/uniflowpbc.hpp"
#include "extensions/wallinducedshear.hpp"
namespace dpd{


class Interactions{
private:
    Control* control;
    Topology* topology;
    Configuration* config;
    Decomposition* decomp;
    SetMPI* mpi;
    NonbondedList nblist;
    BondedList bdlist;
    ExternalList exlist;
    ExtensionList extensionlist;
    
    
public:
    Interactions(){}
    Interactions(Control* control, Topology* topology, Configuration* config, Decomposition* decomp);
    ~Interactions(){}

    void addNonbondedInteraction();
    void addBondedInteraction();
    void addExternalInteraction();
    void addExtension();
    void applyExtensionForPosition(int step);
    void applyExtensionForVelocity(int step);
    NonbondedList& getNonbondedInteractions(){ return nblist; }
    BondedList& getBondedInteractions(){ return bdlist; }
    ExternalList& getExternalInteractions(){ return exlist; }
    ExtensionList& getExtension(){ return extensionlist;}
};
};


#endif
