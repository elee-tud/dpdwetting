/****************************************************************************
 *                             Class Particle
 *
 * This class contains information of each particle, e.g., position, velocity,
 * forces, and so on.
 ****************************************************************************/


#ifndef __PARTICLE__HPP
#define __PARTICLE__HPP

#include "real3d.hpp"
#include <string>
#include <vector>
#include <mpi.h>

namespace dpd{

class Particle;
typedef std::vector<Particle*> ParticleList;
typedef ParticleList::iterator ParticleIterator;
typedef std::vector<ParticleList> Particle2DList;
class Particle{
private:
    int particleindex;  /*particle index*/
    int molindex;       /*molecule index*/
    int atype;          /*particle type number*/
    std::string typestring;     /*particle type*/
    std::string molstring;      /*molecule type*/
    ParticleList bonded;        /*particles bonded to this*/
    ParticleList ss_bonded;     /*particles bonded by slip-springs to this*/
    Ivec bondtypes;             /*bond types*/
    real mass;                  /*mass of the particle*/
    bool is_frozen;             /*is this frozen in space?*/
    int exists;                 /*exists in a cell?*/
    bool is_pinned;             /*is pinned in space?*/
    int cellnum;                /*the cell number that this is placed*/
    
    


public:

    Real3D coord;           /*position*/
    Real3D veloc;           /*velocity*/
    Real3D force;           /*force*/
    Real3D prevcoord;       /*position in the previous step*/
    Real3D prevveloc;       /*velocity in the previous step*/
    real density;           /*local particle density*/
    Rvec stress;            /*virial*/

    Particle(){}            /*Constructor*/
    Particle(int particleindex, int molindex, int atype, real mass, bool is_frozen);

    ~Particle(){}           /*Destructor*/

    void inTheCell(int cellnum);        /*Setting up the cell number of this*/
    int inWhichCell(){ return cellnum; }        /*the cell number where this is placed*/
    void setTypeString(std::string _typestring);        /*assigning the particle type*/
    void setMoleculeString(std::string _molstring);     /*assigning molecule type*/
    void addBond(Particle* bondedparticle, int bondtype);       /*adding bonded particle*/
    
    void setCoord(Real3D newcoord);         /*assigning the position*/
    void setVeloc(Real3D newveloc);         /*assigning the velocity*/
    void setForce(Real3D newforce);         /*assigning the force*/

    /*Returning private variables*/
    Real3D getCoord(){ return coord; }      
    Real3D getVeloc(){ return veloc; }      
    Real3D getForce(){ return force; }      
    real getDensity(){ return density; }    
    Rvec getStress(){ return stress; }      
    int getParticleIndex(){ return particleindex; }     
    int getParticleType(){ return atype; }              
    real getMass(){ return mass; }                      
    bool isFrozen(){ return is_frozen; }                
    ParticleList& getBonds(){ return bonded;}
    Ivec& getBondTypes(){ return bondtypes; }
    int getMoleculeIndex(){ return molindex; }
    std::string getMoleculeName(){ return molstring; }
    std::string getTypeName(){ return typestring; }
    inline int existsHere(){ return exists; }
    bool isPinned(){return is_pinned; }     


    /*Setting up the existance of the particle*/
    inline void setTrue(){ exists=TRUEPTCL; }
    inline void setGhost(){ exists=GHOSTPTCL; }
    inline void setNone(){ exists=NOTEXIST; }

    
    void setPinned();   /*Setting up as pinned*/

    void reserveSSBond(int maxbnum);        /*reserving the slots for slip-springs*/
    void addSSBond(Particle* bondedparticle);       /*adding the slip-spring to this*/
    void removeSSBond(Particle* bondedparticle);    /*removing the slip-spring from this*/
    void removeAllSSBond();                         /*removing all slip-springs from this*/
    ParticleList& getSSBonds(){ return ss_bonded; }     /*Returning slip-springs to this*/


};

Rvec serializeNCoords(Ivec index, ParticleList particles);      /*serialization of the position for communication*/
Rvec serializeNVelocs(Ivec index, ParticleList particles);      /*serialization of the velocity for communication*/
Rvec serializeNForces(Ivec index, ParticleList particles);      /*serialization of the force for communication*/
Rvec serializeNStresses(Ivec index, ParticleList particles);    /*serialization of the Virial for communication*/
Rvec serializeNDensities(Ivec index, ParticleList particles);   /*serialization of the density for communication*/
void deserializeNCoords(Ivec index, ParticleList particles, Rvec serial);      /*de-serialization of the position for communication*/
void deserializeNVelocs(Ivec index, ParticleList particles, Rvec serial);      /*de-serialization of the velocity for communication*/
void deserializeNForces(Ivec index, ParticleList particles, Rvec serial);      /*de-serialization of the force for communication*/
void deserializeNStresses(Ivec index, ParticleList particles, Rvec serial);    /*de-serialization of the Virial for communication*/
void deserializeNDensities(Ivec index, ParticleList particles, Rvec serial);   /*de-serialization of the density for communication*/
};
#endif
