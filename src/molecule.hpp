/****************************************************************************
 *                             Class Molecule
 *
 * This class contains information of molecules and their topology
 ****************************************************************************/


#ifndef __MOLECULE__HPP
#define __MOLECULE__HPP
#include "types.hpp"

namespace dpd{

class Molecule{
private:
    std::string name;   /*name of the moleucle*/
    int num_mols;       /*the number of molecules*/
    int num_atoms;      /*the number of atoms in a molecule*/
    int num_bonds;      /*the number of bonds in a molecule*/
    Ivec atoms;         /*the list of particle indices in a molecule*/
    Rvec masses;        /*the list of particle masses in a molecule*/
    Ivec2D bonds;       /*the list of bonds in a molecule*/
    Ivec atomtypes;     /*the list of particle types in a molecule*/
    Svec atomtypestrings;       /*the list of particle types in string*/
    Ivec bondtypes;             /*the list of types of bonds in a molecule*/
    bool is_frozen;             /*is this molecule frozen?*/

public:
    Molecule(){}            /*Constructor*/
    Molecule(std::string name, int num_mols);

    ~Molecule(){}           /*Destructor*/

    /*Assigning values to private variables*/
    void setName(std::string _name){ name=_name;}
    void setNMolecules(int _num_mols){ num_mols=_num_mols;}

    /*Returning values of private variables*/
    std::string getName(){ return name; }
    int getNatoms(){ return atoms.size(); }
    int getNbonds(){ return bonds.size(); }
    int getNmolecules() {return num_mols; }
    Ivec& getAtoms(){ return atoms; }
    Ivec& getAtomTypes(){ return atomtypes; }
    Svec& getAtomTypeStrings(){ return atomtypestrings; }
    Rvec& getMasses(){ return masses; }
    void addAtoms(int index, std::string atypestring, int atype, real mass);
    void addBonds(int index1, int index2, int btype);
    Ivec2D& getBonds() {return  bonds;}
    Ivec& getBondTypes(){ return bondtypes; }
    bool isFrozen(){ return is_frozen; }

    void printMolecule();       /*printing out the molecule information on a screen*/
    void setFrozen();           /*setting as frozen*/
    
    friend class MPIClasses;    /*MPIClasses can access private functions and variables of this class*/

};

typedef std::vector<Molecule*> MoleculeList;        /*container of molecules*/
};

#endif
