#include "selection.hpp"
#include <algorithm>
#include "parsing.hpp"

using namespace dpd;

Selection::Selection(ParticleList& particles):particles(particles){


}

ParticleList Selection::selectByMoleculeName(std::string molname){
    ParticleList ptclreturn;

//    std::cout << molname << "," << particles[0]->getMoleculeName() << std::endl;
    for(ParticleIterator ptcl=particles.begin();ptcl<particles.end();ptcl++){
        if((*ptcl)->getMoleculeName().compare(molname)==0){
            ptclreturn.push_back(*ptcl);
        }
    }
    return ptclreturn;
}

ParticleList Selection::selectByParticleType(std::string ptclname){
    ParticleList ptclreturn;
    for(ParticleIterator ptcl=particles.begin();ptcl<particles.end();ptcl++){
        if((*ptcl)->getTypeName().compare(ptclname)==0){
            ptclreturn.push_back(*ptcl);
        }
    }
    return ptclreturn;
}

Particle2DList Selection::selectMoleculesByMoleculeName(std::string molname){
    ParticleList selptcls=selectByMoleculeName(molname);
    Ivec molindex;
    for(int i=0;i<selptcls.size();i++){
        int idx=selptcls[i]->getMoleculeIndex();
        if(std::find(molindex.begin(), molindex.end(), idx)==molindex.end())
            molindex.push_back(idx);
    }

    int num_mols=molindex.size();
    Particle2DList molecules(num_mols, ParticleList{});
    for(int i=0;i<selptcls.size();i++){
        int idx=std::distance(molindex.begin(), std::find(molindex.begin(), molindex.end(), selptcls[i]->getMoleculeIndex()));
        molecules[idx].push_back(selptcls[i]);
    }



    return molecules;
}

