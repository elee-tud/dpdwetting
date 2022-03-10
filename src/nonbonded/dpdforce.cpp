#include "dpdforce.hpp" 

using namespace dpd;
DPDForce::DPDForce(Topology* topol, Configuration* config, Decomposition* decomp):Nonbonded(topol, config, decomp){
    B=topol->getB();
    Brcut=topol->getBrcut();
    gamma=control->getGamma();
    invsqrtdt=1./sqrt(control->getTimeStep());
    temp=control->getTemperature();
    q=sqrt(2*gamma*temp);
    randseed=control->getRandomSeed();
    generator=std::default_random_engine(randseed*mpi->rank());
    normdist=std::normal_distribution<real>(0.0, 1.);
    

}

void DPDForce::calculatePairForce(int index1, int index2){
    Real3D rij=pbc.getMinimumImageVector(particles[index2]->coord, particles[index1]->coord);
    real dist=rij.abs();
    Real3D eij=rij.unit();
    Real3D vij=particles[index2]->veloc-particles[index1]->veloc;
    
    int bt1=particles[index1]->getParticleType();
    int bt2=particles[index2]->getParticleType();
    Real3D fij(0.0, 0.0, 0.0);
    if(dist<Brcut[bt1][bt2]){
        real wc=1.-dist/Brcut[bt1][bt2];
        real wd=wc*wc;
        real theta=normdist(generator);
        fij=B[bt1][bt2]*wc*eij;             //Repulsive
        if(integrator!=EMIN){
            fij+=-(gamma*wd*(vij*eij))*eij;     //Dissipative force
            fij+=q*theta*wc*invsqrtdt*eij;              //Random force
        }
    }
    
    if(fij.abs()>0.0){
        particles[index2]->force+=fij;
        particles[index1]->force-=fij;

        calculateStress(particles[index1], particles[index2], rij, fij);
    }
    return;
}

void DPDForce::calculatePairForce(int index1, int index2, Real3D com, real** press, real dr){
    Real3D rij=pbc.getMinimumImageVector(particles[index2]->coord, particles[index1]->coord);
    real dist=rij.abs();
    Real3D eij=rij.unit();
    Real3D vij=particles[index2]->veloc-particles[index1]->veloc;
    
    int bt1=particles[index1]->getParticleType();
    int bt2=particles[index2]->getParticleType();
    Real3D fij(0.0, 0.0, 0.0);
    if(dist<Brcut[bt1][bt2]){
        real wc=1.-dist/Brcut[bt1][bt2];
        real wd=wc*wc;
        real theta=normdist(generator);
        fij=B[bt1][bt2]*wc*eij;             //Repulsive
        if(integrator!=EMIN){
            fij+=-(gamma*wd*(vij*eij))*eij;     //Dissipative force
            fij+=q*theta*wc*invsqrtdt*eij;              //Random force
        }
    }

    if(fij.abs()>0.0){
        particles[index2]->force+=fij;
        particles[index1]->force-=fij;
        calculateSphericalStress(particles[index1], particles[index2], rij, fij, com, press, dr);
    }
    return;
}
