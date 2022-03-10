#ifndef __TOPOLOGY__HPP
#define __TOPOLOGY__HPP

#include <fstream>
#include "types.hpp"
#include "setmpi.hpp"
#include "molecule.hpp"
#include "command.hpp"
#include "control.hpp"
#include "particle.hpp"

namespace dpd{

class Topology{
protected:
    Command *command;
    Control *control;
    SetMPI* mpi;
    std::string topol_fname;
    std::ifstream topstream;
    int nmol_types;
    int nbead_types;
    int nbond_types;
    int nbeads;
    int nbeads_unfrozen;
    int nbonds;
    real totmass;

    int nonbonded;
    int bondlength;

    Svec bead_types;
    Svec bond_types;
    Rvec2D A;
    Rvec2D Arcut;
    Rvec2D B;
    Rvec2D Brcut;
    Rvec bondk;
    Rvec bondl;

    MoleculeList molecules;
    std::string line;
    Svec tokens;

    Ivec2D interbonds;

    ParticleList particles;
    
    
    void openFile();
    void closeFile(){ topstream.close(); }
    std::streampos seekArgument(std::string argument);
    void readMoleculeTypes();
    void readBeadTypes();
    void readBondTypes();
    void readAtomsInMolecule();
    void readBondsInMolecule();
    void readBondsBetweenMolecules();
    void readNonbondedParameters();
    void readBondLengthParameters();
    void exitErrorMissingArgument(std::streampos position, std::string argument);

public:
    Topology(){}
    Topology(Command* command, Control* control, SetMPI* mpi);
    ~Topology(){}

    void readTopology();
    void bcastTopology();
    void bcastMolecules();
    void buildTopology();
    void printSystemInformation();

    void checkDependency();
    
    Rvec2D& getA(){return A; }
    Rvec2D& getArcut(){return Arcut; }
    Rvec2D& getB(){return B; }
    Rvec2D& getBrcut(){return Brcut; }

    Rvec& getBondK(){ return bondk; }
    Rvec& getBondL(){ return bondl; }



    int getNbeads(){ return nbeads;}
    int getNumUnfrozenBeads(){ return nbeads_unfrozen;}
    int getNbonds(){ return nbonds;}
    int getNbondTypes(){ return nbond_types;}
    real getTotalMass(){ return totmass; }
    Command* getCommand(){ return command; }
    Control* getControl(){ return control; }
    SetMPI* getMPI(){ return mpi; }

    ParticleList& getParticles(){ return particles; }
    MoleculeList& getMolecules(){ return molecules; }

};

};
#endif
