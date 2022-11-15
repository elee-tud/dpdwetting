#include "molecule.hpp"
#include <iomanip>
using namespace dpd;
Molecule::Molecule(std::string name, int num_mols):name(name), num_mols(num_mols){
    is_frozen=false;

}

void Molecule::addAtoms(int index, std::string atypestring, int atype, real mass){
    atoms.push_back(index);                     /*adding index of the particle*/
    atomtypestrings.push_back(atypestring);     /*adding particle type as string*/
    atomtypes.push_back(atype);                 /*adding particle type*/
    masses.push_back(mass);                     /*adding mass*/
    num_atoms=atoms.size();                     /*the number of atoms in a molecule*/
    return;
}

void Molecule::addBonds(int index1, int index2, int btype){
    Ivec bd{index1, index2};
    bonds.push_back(bd);                        /*adding the bond*/                    
    bondtypes.push_back(btype);                 /*adding the bond type*/
    num_bonds=bonds.size();                     /*the number of bonds in a molecule*/
    return;
}

void Molecule::setFrozen(){
    is_frozen=true;
}

void Molecule::printMolecule(){
    std::cout << "*****Information in Molecule " << name << "****" << std::endl;
    std::cout << "The number of molecules=" << num_mols << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "<Beads information>" << std::endl;
    std::cout << "#" << std::setw(9) << "Index" << std::setw(10) << "Typs"  << std::setw(10) << "Mass" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    for(int i=0;i<num_atoms;i++){
        std::cout << std::setw(10) << atoms[i] << std::setw(10) << atomtypes[i] << std::setw(10) << masses[i] << std::endl;
    }
    std::cout << "------------------------------------------------------------------" << std::endl;
    if(num_bonds>0){
        std::cout << "<Bonds information>" << std::endl;
        std::cout << "#" << std::setw(9) << "Index1" << std::setw(10) << "Index2"  << std::setw(10) << "Type" << std::endl;
        std::cout << "------------------------------------------------------------------" << std::endl;
        for(int i=0;i<num_bonds;i++){
            std::cout << std::setw(10) << bonds[i][0] << std::setw(10) << bonds[i][1] << std::setw(10) << bondtypes[i] << std::endl;
        }
    }
    return;
}


