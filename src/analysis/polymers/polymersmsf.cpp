#include "polymersmsf.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>
#include "../../filecontrol.hpp"

using namespace dpd;
PolymerStructFactor::PolymerStructFactor(InitialSet initset):Property(initset){
    title="   Calculation of polymer structure factor";

    nprops=1;
    for_distrib=true;
    need_position=true;

    outfile="polymer_smsf.out";


    outtitle_dist="#Polymer structure factor";
    outheader_dist=Svec{"q", "S(q)"};

}





void PolymerStructFactor::initializeVariables(){
    Ivec polidx;

    for(int i=0;i<particles.size();i++){
        if(particles[i]->getMoleculeName().compare("POL")==0){
            int idx=particles[i]->getMoleculeIndex();
            if(std::find(polidx.begin(), polidx.end(), idx)==polidx.end()){
                polidx.push_back(idx);
            }
        }
    }
   
    num_pol=polidx.size();
    for(int i=0;i<num_pol;i++){
        polymers.push_back(new ParticleGroup(pbc));
    }

    for(int i=0;i<particles.size();i++){
        if(particles[i]->getMoleculeName().compare("POL")==0){
            int index=std::distance(polidx.begin(), std::find(polidx.begin(), polidx.end(), particles[i]->getMoleculeIndex()));
            polymers[index]->addParticle(particles[i]);
        }
    }

    for(int i=0;i<num_pol;i++){
        ref.push_back(polymers[i]->at(polymers[i]->size()/2));
    }
    dbin=2*PI/box.abs();
    q.push_back(0);
    q.push_back(dbin);
    real base=1.02;
    real qmax=2*PI/0.1;
    while(true){
        real qinst=q.back()*base;
        if(qinst>qmax)
            break;
        q.push_back(qinst);
    }
    ndbin=q.size();



            
    initializeResultArrays();
    return;
}

void PolymerStructFactor::calculateStep(int step){
    for(int i=0;i<NUMAVG;i++){
        int polindex=rand()%num_pol;
        int nmonpol=polymers[polindex]->size();
        for(int j=0;j<nmonpol;j++){
            for(int k=j+1;k<nmonpol;k++){
                real r=pbc.getMinimumImageVector(polymers[polindex]->at(k)->coord, polymers[polindex]->at(j)->coord).abs();
                dist[0][0]+=1.;
                for(int l=1;l<q.size();l++){
                    real qr=q[l]*r;
                    dist[0][l]+=sin(qr)/qr;
                }
            }
        }
    }

    return;
}

void PolymerStructFactor::normalizeResults(){
    for(int i=0;i<ndbin;i++){
        dist_sum[0][i]/=nsteps*NUMAVG*polymers[0]->size()*(polymers[0]->size()-1)/2;
    }
    return;
}

void PolymerStructFactor::writeOutput(){
    if(mpi->isMaster()){
        std::ofstream outstream;
        openFileWithBackup(outdistname, &outstream, false, false);
        outstream << outtitle_dist << std::endl;
        outstream << outdistheader << std::endl;
        int nrow=dist_sum.size();
        int ncol=dist_sum[0].size();
        for(int i=0;i<ncol;i++){
            outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << q[i];
            outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << dist_sum[0][i];
            outstream << std::endl;
        }
        outstream.close();
        std::cout << "Results are written in " << outdistname << "." << std::endl;
    }
    return;
}

