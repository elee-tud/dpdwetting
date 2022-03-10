#include "slipspringharmonic.hpp"
#include <algorithm>
using namespace dpd;
SlipSpringHarmonic::SlipSpringHarmonic(Topology* topol, Configuration* config, Decomposition* decomp):Bonded(topol, config, decomp){
    control=decomp->getControl();
    ssk=control->getSlipSpringForceConst();
    ssl=control->getSlipSpringLength();
}


Ivec SlipSpringHarmonic::getParticleIndexForCommunication(Particle* ptcl){
    Ivec needed;
    ParticleList ssbonded=ptcl->getSSBonds();
    ParticleList bonded=ptcl->getBonds();
    for(int i=0;i<ssbonded.size();i++){
        if(ptcl->getParticleIndex()<ssbonded[i]->getParticleIndex() && std::find(bonded.begin(), bonded.end(), ssbonded[i])==bonded.end())  {
            if(ssbonded[i]->existsHere()==NOTEXIST){
                ssbonded[i]->force=Real3D(0.0);
                ssbonded[i]->stress=Rvec(9,0.);
                needed.push_back(ssbonded[i]->getParticleIndex()-1);
            }

        }
    }
    return needed;
}





int SlipSpringHarmonic::calculateParticleForce(Particle* ptcl){
    ParticleList ssbonded=ptcl->getSSBonds();
    ParticleList bonded=ptcl->getBonds();
    int calc_nbonds=0;
    for(int i=0;i<ssbonded.size();i++){
        if(ptcl->getParticleIndex()<ssbonded[i]->getParticleIndex() && std::find(bonded.begin(), bonded.end(), ssbonded[i])==bonded.end())  {


            Real3D fij(0.0, 0.0, 0.0);
            Real3D rij=pbc.getMinimumImageVector(ptcl->coord, ssbonded[i]->coord);
            real dist=rij.abs();
            Real3D uij=rij.unit();
            fij=-ssk*(dist-ssl)*uij;
            ptcl->force+=fij;
            ssbonded[i]->force-=fij;

            calculateStress(ptcl, ssbonded[i], rij, fij);

            calc_nbonds++;
        }
    }
    return calc_nbonds;
}

            




int SlipSpringHarmonic::calculateParticleForce(Particle* ptcl, Real3D com, real** press, real dr){
    ParticleList ssbonded=ptcl->getSSBonds();
    int calc_nbonds=0;
    for(int i=0;i<ssbonded.size();i++){
        if(ptcl->getParticleIndex()<ssbonded[i]->getParticleIndex()){
            Real3D fij(0.0, 0.0, 0.0);
            Real3D rij=pbc.getMinimumImageVector(ptcl->coord, ssbonded[i]->coord);
            real dist=rij.abs();
            Real3D uij=rij.unit();
            fij=-ssk*(dist-ssl)*uij;
            ptcl->force+=fij;
            ssbonded[i]->force-=fij;

            calculateSphericalStress(ptcl, ssbonded[i], rij, fij, com, press, dr);

            calc_nbonds++;
        }
    }
    return calc_nbonds;
}

            




