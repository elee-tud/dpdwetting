#include "bondlenharmonic.hpp"
#include <algorithm>

using namespace dpd;
BondLenHarmonic::BondLenHarmonic(Topology* topol, Configuration* config, Decomposition* decomp):Bonded(topol, config, decomp){
    k=topol->getBondK();
    leq=topol->getBondL();


}
Ivec BondLenHarmonic::getParticleIndexForCommunication(Particle* ptcl){
    Ivec needed;
    ParticleList bonded=ptcl->getBonds();
    for(int i=0;i<bonded.size();i++){
        if(ptcl->getParticleIndex()<bonded[i]->getParticleIndex()){
            if(bonded[i]->existsHere()==NOTEXIST){
                bonded[i]->force=Real3D(0.0);
                bonded[i]->stress=Rvec(9,0.);
                needed.push_back(bonded[i]->getParticleIndex()-1);
//                std::cout <<  ptcl->getParticleIndex()-1 << "," << needed.back() << std::endl;
            }

        }
    }
    return needed;
}





int BondLenHarmonic::calculateParticleForce(Particle* ptcl){
    ParticleList bonded=ptcl->getBonds();
    Ivec type=ptcl->getBondTypes();
    int calc_nbonds=0;
    for(int i=0;i<bonded.size();i++){
        if(ptcl->getParticleIndex()<bonded[i]->getParticleIndex()){
//            std::cout << ptcl->getParticleIndex() << "," << bonded[i]->getParticleIndex() << std::endl;
//            if(bonded[i]->existsHere()==NOTEXIST)
//                err.missingBondedParticle(ptcl->getParticleIndex(), bonded[i]->getParticleIndex());
            Real3D fij(0.0, 0.0, 0.0);
            Real3D rij=pbc.getMinimumImageVector(ptcl->coord, bonded[i]->coord);
            real dist=rij.abs();
            Real3D uij=rij.unit();
            fij=-k[type[i]]*(dist-leq[type[i]])*uij;
//            if(mpi->rank()==0){
//          if(ptcl->getParticleIndex()==1){
//          std::cout << ptcl->getParticleIndex() << "," << bonded[i]->getParticleIndex() ;
//            std::cout << rij << dist << fij << std::endl;
//            }
            ptcl->force+=fij;
            bonded[i]->force-=fij;

            calculateStress(ptcl, bonded[i], rij, fij);

            calc_nbonds++;
        }
    }
    return calc_nbonds;
}

            




int BondLenHarmonic::calculateParticleForce(Particle* ptcl, Real3D com, real** press, real dr){
    ParticleList bonded=ptcl->getBonds();
    Ivec type=ptcl->getBondTypes();
    int calc_nbonds=0;
    for(int i=0;i<bonded.size();i++){
        if(ptcl->getParticleIndex()<bonded[i]->getParticleIndex()){
//            if(bonded[i]->existsHere()==NOTEXIST)
//                err.missingBondedParticle(ptcl->getParticleIndex(), bonded[i]->getParticleIndex());
            Real3D fij(0.0, 0.0, 0.0);
            Real3D rij=pbc.getMinimumImageVector(ptcl->coord, bonded[i]->coord);
            real dist=rij.abs();
            Real3D uij=rij.unit();
            fij=-k[type[i]]*(dist-leq[type[i]])*uij;
//            if(mpi->rank()==0){
 //       std::cout << ptcl->getParticleIndex() << "," << bonded[i]->getParticleIndex() ;
//            std::cout << rij << dist << fij << std::endl;
//            }
            ptcl->force+=fij;
            bonded[i]->force-=fij;

            calculateSphericalStress(ptcl, bonded[i], rij, fij, com, press, dr);

            calc_nbonds++;
        }
    }
    return calc_nbonds;
}

            




