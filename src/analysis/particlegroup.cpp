#include "particlegroup.hpp"
using namespace dpd;
ParticleGroup::ParticleGroup(PeriodicBoundary pbc):pbc(pbc){
    particles=ParticleList{};
}

ParticleGroup::ParticleGroup(ParticleList particles, PeriodicBoundary pbc):particles(particles), pbc(pbc){
    ref=particles[0];
}



void ParticleGroup::addParticle(Particle* ptcl){
    particles.push_back(ptcl);
    return;
}

void ParticleGroup::setReference(Particle* _ref){
    ref=_ref;
    return;
}
void ParticleGroup::setReference(Real3D ref0, Real3D ref1){
    if(particles.size()!=0)
        ref=particles[0];

    for(int i=0;i<particles.size();i++){
        if(particles[i]->coord[0]>=ref0[0] and
                particles[i]->coord[0]<ref1[0] and
                particles[i]->coord[1]>=ref0[1] and
                particles[i]->coord[1]<ref1[1] and
                particles[i]->coord[2]>=ref0[2] and
                particles[i]->coord[2]<ref1[2]){
            ref=particles[i];
            break;
        }
    }
    return;
    
}


Real3D ParticleGroup::calculateCenterOfMass(Particle* _ref){
    ref=_ref;
    com=Real3D{0.0, 0.0, 0.0};
    for(int i=0;i<particles.size();i++){
        Real3D dr=pbc.getMinimumImageVector(particles[i]->coord, ref->coord);
        com+=dr;
    }
    com=com/particles.size()+ref->coord;
    return com;
}

Real3D ParticleGroup::calculateCenterOfMass(){
    com=Real3D{0.0, 0.0, 0.0};
    for(int i=0;i<particles.size()-1;i++){
        Real3D dr=pbc.getMinimumImageVector(particles[i+1]->coord, particles[i]->coord);
        com+=dr*(particles.size()-i-1);
    }
    com=com/particles.size()+particles[0]->coord;
    return com;
}

Real3D ParticleGroup::calculateCenterOfMass(Real3D ref0, Real3D ref1){
    com=Real3D{0.0, 0.0, 0.0};
    setReference(ref0, ref1);
    com=calculateCenterOfMass(ref);
    return com;
}

void ParticleGroup::printParticleIndices(){
    std::cout << "[" ;
    for(int i=0;i<particles.size()-1;i++){
        std::cout << particles[i]->getParticleIndex() << "," ;
    }
    std::cout << particles.back()->getParticleIndex() << "]" << std::endl;
    return;
}



real ParticleGroup::calculateRadiusOfGyration(Particle* _ref){
    ref=_ref;
    return calculateRadiusOfGyration();
}
        

real ParticleGroup::calculateRadiusOfGyration(){

//    calculateCenterOfMass();
    real rgsqr=0.;
    for(int i=0;i<particles.size();i++){
        Real3D rdiff(0.);
        for(int j=i+1;j<particles.size();j++){
            rdiff+=pbc.getMinimumImageVector(particles[j]->coord, particles[j-1]->coord);
            rgsqr+=rdiff.sqr();
        }
    }
    rgsqr/=particles.size()*particles.size();
    return rgsqr;
}
       
Rvec ParticleGroup::calculateGyrationTensor(Particle* _ref){
    ref=_ref;
    return calculateGyrationTensor();
}

Rvec ParticleGroup::calculateGyrationTensor(){
//    calculateCenterOfMass(ref);
    Rvec gt(6,0.0);
    for(int i=0;i<particles.size();i++){
        Real3D rdiff(0.);
        for(int j=i+1;j<particles.size();j++){
            rdiff+=pbc.getMinimumImageVector(particles[j]->coord, particles[j-1]->coord);
            gt[0]+=rdiff[0]*rdiff[0];
            gt[1]+=rdiff[1]*rdiff[1];
            gt[2]+=rdiff[2]*rdiff[2];
            gt[3]+=rdiff[0]*rdiff[1];
            gt[4]+=rdiff[1]*rdiff[2];
            gt[5]+=rdiff[2]*rdiff[3];
        }
    }
    for(int i=0;i<6;i++){
        gt[i]/=particles.size()*particles.size();
    }
    return gt;
}
        

