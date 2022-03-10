#ifndef __SELECTION__HPP
#define __SELECTION__HPP
#include "particle.hpp"
namespace dpd{


class Selection{
private:
    ParticleList particles;

public:
    Selection(){}
    Selection(ParticleList& particles);
    ~Selection(){}

    ParticleList selectByMoleculeName(std::string molname);
    ParticleList selectByParticleType(std::string ptclname);
    Particle2DList selectMoleculesByMoleculeName(std::string molname);
};

};


#endif
