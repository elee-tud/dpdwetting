#ifndef __INITANALYSIS__HPP
#define __INITANALYSIS__HPP
#include "../initialization.hpp"
#include "../command.hpp"
#include "trajectory.hpp"

#include "droplet/dropsize.hpp"
#include "droplet/dropvelocity.hpp"
#include "droplet/dropsphericalstress.hpp"
#include "droplet/radialdensity.hpp"
#include "droplet/dropzdensity.hpp"
#include "droplet/surfacerdf.hpp"
#include "droplet/jumpfreq.hpp"
#include "droplet/polymerdepositorient.hpp"
#include "droplet/polymeradsorption.hpp"
#include "droplet/surfacecoverage.hpp"
#include "droplet/polymerstretch.hpp"

#include "capillarybridge/bridgesize.hpp"
#include "capillarybridge/bridgevelocity.hpp"
#include "capillarybridge/bridgepolconc.hpp"
#include "capillarybridge/bridgeadsconc.hpp"
#include "capillarybridge/bridgeclvel.hpp"
#include "capillarybridge/bridgejumpfreq.hpp"
#include "capillarybridge/bridgeslipvel.hpp"
#include "capillarybridge/bridgevelx.hpp"
#include "capillarybridge/bridgevelxz.hpp"
#include "capillarybridge/bridgealdens.hpp"
#include "capillarybridge/bridgezdensity.hpp"
#include "capillarybridge/bridgedensxz.hpp"

#include "polymers/polymersize.hpp"
#include "polymers/polymereevecrelax.hpp"
#include "polymers/polymerorientation.hpp"
#include "polymers/polymerbondlength.hpp"
#include "polymers/msd.hpp"
#include "polymers/polymersmsf.hpp"
#include "polymers/polymerinternalsize.hpp"


#include "general/averagestress.hpp"
#include "general/trjtogro.hpp"
#include "general/velacf.hpp"
#include "general/rdf.hpp"


namespace dpd{

class InitAnalysis:public Initialization{
private:
    int program;
    Property* prop;
    InitialSet initset;
    Trajectory* traj;


    bool multicore;
    
    int nbeads;
    int nsteps;
    int begstep;
    int endstep;
    int skipstep;


    void loadProperty();
    void checkNumberOfCoresForProgram();

public:
    InitAnalysis(){}
    InitAnalysis(Command* command);
    ~InitAnalysis();

    int getNbeads(){ return nbeads; }
    int getNsteps(){ return nsteps; }
    int getBeginStep(){ return begstep; }
    int getEndStep(){ return endstep; }
    int getSkipStep(){ return skipstep; }
    bool needMultiCore(){ return multicore; }

    Property* getProperty(){ return prop; }
    Trajectory* getTrajectory(){ return traj; }




};
};


#endif
