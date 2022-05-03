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
    int particleindex;
    int molindex;
    int atype;
    std::string typestring;
    std::string molstring;
    ParticleList bonded;
    ParticleList ss_bonded;
    Ivec bondtypes;
    real mass;
    bool is_frozen;
    int exists;
    bool is_pinned;
    int cellnum;
    
    


public:

    Real3D coord;
    Real3D veloc;
    Real3D force;
    Real3D prevcoord;
    Real3D prevveloc;
    real density;
    Rvec stress;

    Particle(){}
    Particle(int particleindex, int molindex, int atype, real mass, bool is_frozen);

    ~Particle(){}

    void inTheCell(int cellnum);
    int inWhichCell(){ return cellnum; }
    void setTypeString(std::string _typestring);
    void setMoleculeString(std::string _molstring);
    void addBond(Particle* bondedparticle, int bondtype);
    
    void setCoord(Real3D newcoord);
    void setVeloc(Real3D newveloc);
    void setForce(Real3D newforce);

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
    inline void setTrue(){ exists=TRUEPTCL; }
    inline void setGhost(){ exists=GHOSTPTCL; }
    inline void setNone(){ exists=NOTEXIST; }

    void setPinned();
    bool isPinned(){return is_pinned; }

    void reserveSSBond(int maxbnum);
    void addSSBond(Particle* bondedparticle);
    void removeSSBond(Particle* bondedparticle);
    void removeAllSSBond();
    ParticleList& getSSBonds(){ return ss_bonded; }


};

Rvec serializeNCoords(Ivec index, ParticleList particles);
Rvec serializeNVelocs(Ivec index, ParticleList particles);
Rvec serializeNForces(Ivec index, ParticleList particles);
Rvec serializeNStresses(Ivec index, ParticleList particles);
Rvec serializeNDensities(Ivec index, ParticleList particles);
void deserializeNCoords(Ivec index, ParticleList particles, Rvec serial);
void deserializeNVelocs(Ivec index, ParticleList particles, Rvec serial);
void deserializeNForces(Ivec index, ParticleList particles, Rvec serial);
void deserializeNStresses(Ivec index, ParticleList particles, Rvec serial);
void deserializeNDensities(Ivec index, ParticleList particles, Rvec serial);
};
#endif
