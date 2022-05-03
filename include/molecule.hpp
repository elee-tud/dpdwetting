#ifndef __MOLECULE__HPP
#define __MOLECULE__HPP
#include "types.hpp"

namespace dpd{

class Molecule{
private:
    std::string name;
    int num_mols;
    int num_atoms;
    int num_bonds;
    Ivec atoms;
    Rvec masses;
    Ivec2D bonds;
    Ivec atomtypes;
    Svec atomtypestrings;
    Ivec bondtypes;
    bool is_frozen;

public:
    Molecule(){}
    Molecule(std::string name, int num_mols);

    ~Molecule(){}
    void setName(std::string _name){ name=_name;}
    void setNMolecules(int _num_mols){ num_mols=_num_mols;}
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
    void printMolecule();
    void setFrozen();
    bool isFrozen(){ return is_frozen; }
    
    friend class MPIClasses;

};

typedef std::vector<Molecule*> MoleculeList;
};

#endif
