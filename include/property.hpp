#ifndef __PROPERTY__HPP 
#define __PROPERTY__HPP

#define MAXAVG 100000
#include "../initialization.hpp"

namespace dpd{
class Property{
protected:
    InitialSet initset;
    ParticleList particles;
    Command *command;
    SetMPI* mpi;
    Control* control;
    Topology* topol;
    Configuration* config;
    Interactions* inter;
    Real3D box;
    PeriodicBoundary pbc;

    int property;
    bool for_timeevol;
    bool for_distrib;
    bool for_dynamics;
    bool need_multicore;
    bool need_position;
    bool need_velocity;
    bool need_stress;
    bool need_force;

    std::string title;
    int nprops;
    
    real dbin;
    real dbegin;
    int ndbin;
    int numvarperstep;

    std::string outfile;
    std::string outtitle_tevol;
    std::string outtitle_dist;
    std::string outtitle_dynm;
    Svec outheader_tevol;
    Svec outheader_dist;
    Svec outheader_dynm;

    int nbeads;
    int nsteps;
    int begstep;
    int endstep;
    int skipstep;

    float dt;


    Ivec dsteps;
    int ndsteps;
    Ivec numavgs;

    Rvec2D tevol;
    Rvec2D tevol_sum;

    Rvec2D dist;
    Rvec2D dist_sum;
    

    Svec varname;
    
    std::string outname;
    std::string outtimename;
    std::string outdistname;
    std::string outdynmname;
    std::string outtimeheader;
    std::string outdistheader;
    std::string outdynmheader;
    void generateOutputTitleHeader();
    std::string generateOutputHeader(Svec hd);

    void writeTimeEvolution();
    void writeDistribution();
    void writeDynamics();

public:
    Property(){}
    Property(InitialSet initset);
    ~Property(){}
  

    int getProperty(){return property;}
    void printPropertyTitle();
    void getGeneralParameters();
    virtual void getSpecificParameters()=0;
    virtual void initializeVariables()=0;
    void initializeResultArrays();
    virtual void calculateStep(int step)=0;
    void generateDeltaSteps(); 
    virtual void reduceResults();
    virtual void normalizeResults();
    virtual void calculateDynamicProperty(){}
    virtual void generateOutFileName();
    virtual void writeOutput();

    virtual void dumpGro(int step, real t){}


    void reduceTimeEvolution();
    void reduceDistribution();
    void resetBox(Real3D newbox);



    bool forTimeEvolution(){ return for_timeevol; }
    bool forDistribution(){ return for_distrib; }
    bool forDynamics(){ return for_dynamics; }
    bool needPosition(){ return need_position;}
    bool needVelocity(){ return need_velocity;}
    bool needStress(){ return need_stress; }
    bool needForce(){ return need_force; }
    bool needMultiCore(){ return need_multicore; }

    int getTrajectoryType();

    int getNbeads(){ return nbeads; }
    int getNsteps(){ return nsteps; }
    int getBeginStep(){ return begstep; }
    int getEndStep(){ return endstep; }
    int getSkipStep(){ return skipstep; }

    std::string getOutFileName(){ return outfile; }
    std::string getOutTimeTitle(){ return outtitle_tevol;}
    std::string getOutDistTitle(){ return outtitle_dist; }
    std::string getOutDynmTitle(){ return outtitle_dynm; }
    Svec getOutTimeHeader(){ return outheader_tevol; }
    Svec getOutDistHeader(){ return outheader_dist; }
    Svec getOutDynmHeader(){ return outheader_dynm; }
    SetMPI* getMPI(){ return mpi; }
    Command* getCommand(){ return command; }
    Control* getControl(){return control; }
    Rvec2D getTimeEvolution(){ return tevol_sum; }
    Rvec2D getDistribution(){ return dist_sum; }
    real getDistFirst(){ return dbegin; }
    real getDistBin(){ return dbin; }

    Ivec getDeltaSteps(){ return dsteps; }


};
};


#endif
