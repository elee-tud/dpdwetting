#include "particle.hpp"
#include <algorithm>

using namespace dpd;
Particle::Particle(int particleindex, int molindex, int atype, real mass, bool is_frozen):particleindex(particleindex), molindex(molindex), atype(atype), mass(mass), is_frozen(is_frozen){
    coord=Real3D(0.0);
    veloc=Real3D(0.0);
    prevcoord=Real3D(0.0);
    prevveloc=Real3D(0.0);
    force=Real3D(0.0);
    stress=Rvec(9, 0.0);
    bonded.reserve(4);
    bondtypes.reserve(4);
    density=0.0;
    is_pinned=false;

}
void Particle::inTheCell(int cidx){
    cellnum=cidx;
    return;
}

void Particle::setTypeString(std::string _typestring){ 
    typestring=_typestring; 
    return;
}

void Particle::setMoleculeString(std::string _molstring){ 
    molstring=_molstring; 
    return;
}

void Particle::addBond(Particle* bondedparticle, int bondtype){
    bonded.push_back(bondedparticle);
    bondtypes.push_back(bondtype);
    return;
}


void Particle::setCoord(Real3D newcoord){
    coord=newcoord;
    return;
}

void Particle::setVeloc(Real3D newveloc){
    veloc=newveloc;
    return;
}

void Particle::setForce(Real3D newforce){
    veloc=newforce;
    return;
}

void Particle::setPinned(){
    is_pinned=true;
    return;
}

void Particle::reserveSSBond(int maxbnum){
    ss_bonded.reserve(maxbnum);
    return;
}

void Particle::addSSBond(Particle* bondedparticle){
    ss_bonded.push_back(bondedparticle);
    return;
}
void Particle::removeSSBond(Particle* bondedparticle){
    ParticleIterator it=std::find(ss_bonded.begin(), ss_bonded.end(), bondedparticle);
    if(it==ss_bonded.end())
        MPI_Abort(MPI_COMM_WORLD, 20);
    else
        ss_bonded.erase(it);

    return;
}

void Particle::removeAllSSBond(){
    ss_bonded.clear();
    return;
}


namespace dpd{
Rvec serializeNCoords(Ivec index, ParticleList particles){
    int size=index.size();
    Rvec out(size*3);
    for(int i=0;i<size;i++){
        out[i*3]=particles[index[i]]->coord[0];
        out[i*3+1]=particles[index[i]]->coord[1];
        out[i*3+2]=particles[index[i]]->coord[2];
    }
    return out;
}

Rvec serializeNVelocs(Ivec index, ParticleList particles){
    int size=index.size();
    Rvec out(size*3);
    for(int i=0;i<size;i++){
        out[i*3]=particles[index[i]]->veloc[0];
        out[i*3+1]=particles[index[i]]->veloc[1];
        out[i*3+2]=particles[index[i]]->veloc[2];
    }
    return out;
}

Rvec serializeNForces(Ivec index, ParticleList particles){
    int size=index.size();
    Rvec out(size*3);
    for(int i=0;i<size;i++){
        out[i*3]=particles[index[i]]->force[0];
        out[i*3+1]=particles[index[i]]->force[1];
        out[i*3+2]=particles[index[i]]->force[2];
    }
    return out;
}

Rvec serializeNStresses(Ivec index, ParticleList particles){
    int size=index.size();
    Rvec out(size*9);
    for(int i=0;i<size;i++){
        out[i*9]=particles[index[i]]->stress[0];
        out[i*9+1]=particles[index[i]]->stress[1];
        out[i*9+2]=particles[index[i]]->stress[2];
        out[i*9+3]=particles[index[i]]->stress[3];
        out[i*9+4]=particles[index[i]]->stress[4];
        out[i*9+5]=particles[index[i]]->stress[5];
        out[i*9+5]=particles[index[i]]->stress[6];
        out[i*9+5]=particles[index[i]]->stress[7];
        out[i*9+5]=particles[index[i]]->stress[8];
    }
    return out;
}

Rvec serializeNDensities(Ivec index, ParticleList particles){
    int size=index.size();
    Rvec out(size);
    for(int i=0;i<size;i++){
        out[i]=particles[index[i]]->density;
    }
    return out;
}

void deserializeNCoords(Ivec index, ParticleList particles, Rvec serial){
    int size=index.size();
    for(int i=0;i<size;i++){
        particles[index[i]]->coord[0]=serial[i*3];
        particles[index[i]]->coord[1]=serial[i*3+1];
        particles[index[i]]->coord[2]=serial[i*3+2];
    }
    return;
}


void deserializeNVelocs(Ivec index, ParticleList particles, Rvec serial){
    int size=index.size();
    for(int i=0;i<size;i++){
        particles[index[i]]->veloc[0]=serial[i*3];
        particles[index[i]]->veloc[1]=serial[i*3+1];
        particles[index[i]]->veloc[2]=serial[i*3+2];
    }
    return;
}


void deserializeNForces(Ivec index, ParticleList particles, Rvec serial){
    int size=index.size();
    for(int i=0;i<size;i++){
        particles[index[i]]->force[0]=serial[i*3];
        particles[index[i]]->force[1]=serial[i*3+1];
        particles[index[i]]->force[2]=serial[i*3+2];
    }
    return;
}

void deserializeNStresses(Ivec index, ParticleList particles, Rvec serial){
    int size=index.size();
    for(int i=0;i<size;i++){
        particles[index[i]]->stress[0]=serial[i*9];
        particles[index[i]]->stress[1]=serial[i*9+1];
        particles[index[i]]->stress[2]=serial[i*9+2];
        particles[index[i]]->stress[3]=serial[i*9+3];
        particles[index[i]]->stress[4]=serial[i*9+4];
        particles[index[i]]->stress[5]=serial[i*9+5];
        particles[index[i]]->stress[6]=serial[i*9+6];
        particles[index[i]]->stress[7]=serial[i*9+7];
        particles[index[i]]->stress[8]=serial[i*9+8];
    }
    return;
}


void deserializeNDensities(Ivec index, ParticleList particles, Rvec serial){
    int size=index.size();
    for(int i=0;i<size;i++){
        particles[index[i]]->density=serial[i];
    }
    return;
}
}


