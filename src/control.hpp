/*************************************************************
 *           Class Control
 *
 *  This class contains functions and variables for control 
 *  variables. 
 *************************************************************/
#ifndef __CONTROL__HPP
#define __CONTROL__HPP

#include <fstream>
#include "real3d.hpp"
#include "setmpi.hpp"
#include "command.hpp"

namespace dpd{


class Control{
protected:
    Command* command;
    SetMPI *mpi;

    std::string control_fname;
    std::ifstream constream;

    /* Control variables
     * If one wants to add varialbes to read from the
     * control file, add variables here 
     * read it from Control::readControl(). They
     * have to be broadcasted in Control:bcastControl()
     * If necessary, set the default values in 
     * Control::setDefaults().
     */
    bool restart;
    real temperature;           /*Simulation temperature*/
    real dt;                    /*Integration time step*/
    real gamma;                 /*Friction constant*/
    real lambda;                /*Coupling constant*/
    real cellcut;               /*Cell cutoff distance*/
    real saferatio;             /*Ratio to multiply for the number of particles in each cell*/
    int randomseed;             /*Seed for random number generator*/
    int integrator;             /*Type of integrator (see types.hpp)*/
    int totalsteps;             /*Total integration steps*/

    bool dumpbinary;            /*True, if trajectory is written binary*/
    int trajfreq;               /*Frequency to write trajectory*/
    int stressfreq;             /*Frequency to write stress*/
    int logfreq;                /*Frequency to write log*/
    int forcefreq;              /*Frequency to write force*/

    int nonbonded;              /*Type of nonbonded potential (see types.hpp)*/
    int bondlength;             /*Type of bonded potential (see types.hpp)*/

    int walltype;               /*Type of solid wall (see types.hpp)*/
    Ivec walldirect;            /*Direction of wall*/
    Rvec wallmin;               /*Position of bottom walls*/
    Rvec wallmax;               /*Position of top walls*/

    int rmcomvfreq;             /*Frequency to remove center-of-mass motion*/
    Ivec rmcomvdir;             /*Direction of removing COM motion*/

    bool dumpfrozen;            /*True, if the trajectory of frozen particles is saved*/

    real targetforce;           /*Target maximum force of energy minimization */

   
    /*Control variables for external fields*/
    real pull_springk;          /*Spring constant of particle pulling*/
    Real3D pull_coord;          /*Position of the center of spring pulling*/
    int pull_direct;            /*Direction of spring pulling*/

    /*Control variables for gravitation field*/
    real gravity_field;         /*Gravity acceleration*/
    int gravity_direct;         /*Direction of gravity*/


    /*Control varialbes for box deformation*/
    int deformation;            /*Type of box deformation*/
    int deformfreq;             /*Frequency to deform the box*/
    real deformfactor;          /*Deformation factor*/

    /*Control variables for particle pinning*/
    Ivec pindirect;             /*Direction of the pinning*/
    std::string pinmolname;     /*Molecule group name to be pinned*/
    Ivec pinindex;              /*Indices of the particles to be pinned*/
    bool pinrandom;             /*True, if the pinning is done to a random particle*/
    int npinpermol;             /*The number of pinned particles per molecule*/


    /*control variables for uniaxial flow*/
    Ivec uniflow;               /*Direction of uniaxial flow*/
    real flowrate;              /*Flow rate*/
    Rvec sheartensor;           /*Shear tensor*/
    void buildShearTensor();    /*Function to construct the shear tensor*/

    /*control variables for shear flow*/
    int sheardir;               /*Direction of shear flow*/
    real shearrate;             /*Shear rate*/
    Svec liquidgrps;            /*Molecule groups to be applied*/
    
    /*control variables for wall-induced shear*/
    std::string shrwallgrp1, shrwallgrp2;       /*Two group names of the walls to be sheared*/
    int wsheardir;                              /*Direction of the shear*/
    real wshearrate;                            /*Wall shear rate*/

    
    Svec tokens;
    std::string line;

    

    /*control variables for slip spring*/
    bool slipspring;                /*True, if slip-springs exist*/
    int numss;                      /*The number of slip-springs*/
    int nsteps_ssmc;                /*The number of MC steps of slip-springs in each sequence*/
    int nsteps_ssdpd;               /*The number of DPD steps in each sequence*/
    real ssl0;                      /*Equilibrium bond length of the slip-spring*/
    real ssk;                       /*Force constance of the slip-spring*/
    real sscutmax;                  /*Maximum distance of the slip-spring to be relocated*/
    real sscutmin;                  /*Minimum distance of the slip-spring to be relocated*/
    int sspoltype;                  /*Type of slip-spring relocation for different polymer topologies*/
    real intrassbias;               /*Bias for intra slip-springs*/


    /*control variables for temperature annealing*/
    real tempannrate;


public:
    Control(){}
    Control(Command* command, SetMPI *mpi);
    Control(std::string control_fname, SetMPI *mpi);
    ~Control(){};

    Command* getCommand(){return command;}
    void setDefaults();             /*Setting-up default values*/
    bool getValues(std::string option, int numvalues, Svec& values);
    void readControl();             /*Reading the control variables*/
    void setDependency();           /*Resolving dependencies*/
    void bcastControl();            /*Broadcasting control parameters*/

    bool getRestart(){return restart;}

    real getTemperature(){ return temperature; }
    real getTempAnnealRate(){ return tempannrate; }
    real getTimeStep(){ return dt;}
    real getGamma(){ return gamma; }
    real getLambda(){ return lambda; }
    real getCellCutoff(){ return cellcut; }
    real getSafeRatio(){ return saferatio;}
    int getRandomSeed(){ return randomseed; }
    int getIntegrator(){ return integrator; }
    int getTotalSteps(){ return totalsteps; }
    int getTrajFrequency(){ return trajfreq; }
    int getLogFrequency(){ return logfreq; }
    int getStressFrequency(){ return stressfreq; }
    int getForceFrequency(){ return forcefreq; }

    int getNonbondedInteraction(){ return nonbonded; }
    int getBondLengthInteraction(){ return bondlength; }
    real getTargetForce(){ return targetforce; }

    bool getDumpFrozen() { return dumpfrozen; }
    bool getDumpBinary(){ return dumpbinary; }

    real getPullingSpringK(){ return pull_springk; }
    Real3D getPullingCoord(){ return pull_coord; }
    int getPullingDirect(){ return pull_direct; }

    real getGravityField(){ return gravity_field; }
    int getGravityDirection(){ return gravity_direct; }

    int getWallType(){ return walltype; }
    Ivec getWallDirection(){ return walldirect; }
    Rvec getWallMinPosition(){ return wallmin; }
    Rvec getWallMaxPosition(){ return wallmax; }

    int getRmCOMVelFrequency(){ return rmcomvfreq; }
    Ivec getRmCOMVelDirection(){ return rmcomvdir; }


    int getDeformation(){ return deformation; }
    int getDeformationFrequency(){ return deformfreq; }
    real getDeformationFactor(){ return deformfactor; }

    bool getPinning();
    Ivec& getPinningDirection(){ return pindirect; }
    std::string& getPinningMoleculeName(){ return pinmolname; }
    Ivec& getPinningIndex(){ return pinindex; }
    bool& getPinningRandom(){ return pinrandom; }
    int& getNumPinningPerMolecule(){ return npinpermol; }

    bool isUniflow();
    Rvec getShearTensor(){ return sheartensor; }

    bool isWallSheared();
    int getWallShearDir(){ return wsheardir; }
    real getWallShearRate(){ return wshearrate; }
    std::string getShearWallGroup1(){ return shrwallgrp1; }
    std::string getShearWallGroup2(){ return shrwallgrp2; }


    bool slipSpring(){ return slipspring; }
    int getNumSlipSprings(){ return numss;}
    int getNumMCSeqSteps(){ return nsteps_ssmc; }
    int getNumDPDSeqSteps(){ return nsteps_ssdpd; }
    real getSlipSpringLength(){ return ssl0; }
    real getSlipSpringForceConst(){ return ssk; }
    real getSlipSpringMaxCutoff(){ return sscutmax; }
    real getSlipSpringMinCutoff(){ return sscutmin; }
    int getSlipSpringPolymerType(){ return sspoltype; }
    int getIntraSlipSpringBias(){ return intrassbias; }


    Svec getLiquidGroups(){ return liquidgrps; }
    void annealingTemperature();


};
}


#endif
