#ifndef __PARTICLEGROUP__HPP
#define __PARTICLEGROUP__HPP
#include "../particle.hpp"
#include "../periodicboundary.hpp"

namespace dpd{
class ParticleGroup;
class ParticleGroup{
private:
    PeriodicBoundary pbc;
    Real3D ref0;
    Real3D ref1;
    ParticleList particles;
    Particle* ref;
    Real3D com;

public:
    ParticleGroup(){}
    ParticleGroup(PeriodicBoundary pbc);
    ParticleGroup(ParticleList particles, PeriodicBoundary dpc);
    ~ParticleGroup(){}
    

    void addParticle(Particle* ptcl);
    void setReference(Real3D ref0, Real3D ref1);
    void setReference(Particle* ref);
    Real3D calculateCenterOfMass(Particle* _ref);
    Real3D calculateCenterOfMass();
    Real3D calculateCenterOfMass(Real3D ref0, Real3D ref1);
    int size(){ return particles.size(); }
    void printParticleIndices();
    ParticleList& getParticles(){ return particles;}
    void clear(){ return particles.clear(); }
    Particle* at(int i){ return particles[i]; }
    real calculateRadiusOfGyration(Particle* _ref);
    real calculateRadiusOfGyration();
    Rvec calculateGyrationTensor(Particle* _ref);
    Rvec calculateGyrationTensor();



};
typedef std::vector<ParticleGroup*> MoleculeGroup;
};

#endif
