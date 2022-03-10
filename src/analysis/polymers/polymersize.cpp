#include "polymersize.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
PolymerSize::PolymerSize(InitialSet initset):Property(initset){
    title="   Calculation of polymer size ";

    nprops=12;
    for_timeevol=true;
    for_distrib=true;
    need_position=true;

    outfile="polymer_size.out";


    outtitle_tevol="#Polymer Radius of gyration and end-to-end vector";
    outheader_tevol=Svec{"Time", "Rg", "Re", "Rs", "Rgx", "Rgy", "Rgz", "Rex", "Rey", "Rez", "Rsx", "Rsy", "Rsz" };
    outtitle_dist="#Distribution of Polymer radius of gyration and end-to-end vectors";
    outheader_dist=Svec{"R^2", "P(Rg^2)", "P(Re^2)", "P(Rs^2)","P(Rgx^2)", "P(Rgy^2)", "P(Rgz^2)", "P(Rex^2)", "P(Rey^2)", "P(Rez^2)", "P(Rsx^2)", "P(Rsy^2)", "P(Psz^2)"};

}



void PolymerSize::getSpecificParameters(){
    molname="POL";
    command->getCommandSingleOption("-mol", molname, &molname);
    std::cout << "Caculation will be done for molecules with name=" << molname << std::endl;

    return;
}



void PolymerSize::initializeVariables(){
    Ivec polidx;

    for(int i=0;i<particles.size();i++){
        if(particles[i]->getMoleculeName().compare(molname)==0){
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
        if(particles[i]->getMoleculeName().compare(molname)==0){
            int index=std::distance(polidx.begin(), std::find(polidx.begin(), polidx.end(), particles[i]->getMoleculeIndex()));
            polymers[index]->addParticle(particles[i]);
        }
    }

    for(int i=0;i<num_pol;i++){
        ref.push_back(polymers[i]->at(static_cast<int>(polymers[i]->size()/2)));
    }

    dbin=0.1;
    dbegin=0.0;
    ndbin=int(10*(box.sqr()/2)/dbin)+1;
    numvarperstep=num_pol;

            
    initializeResultArrays();
    return;
}

void PolymerSize::calculateStep(int step){
    real rgsqr=0.;
    real resqr=0.;
    real rssqr=0.;
    real rgxsqr=0.;
    real rgysqr=0.;
    real rgzsqr=0.;
    real rexsqr=0.;
    real reysqr=0.;
    real rezsqr=0.;
    real rsxsqr=0.;
    real rsysqr=0.;
    real rszsqr=0.;
    for(int i=0;i<num_pol;i++){
        Rvec gt=polymers[i]->calculateGyrationTensor();
        real rgsqrmol=gt[0]+gt[1]+gt[2];
        rgsqr+=rgsqrmol;
        rgxsqr+=gt[0];
        rgysqr+=gt[1];
        rgzsqr+=gt[2];
        dist[0][int(rgsqrmol/dbin)]+=1.;
        dist[2][int(gt[0]/dbin)]+=1.;
        dist[3][int(gt[1]/dbin)]+=1.;
        dist[4][int(gt[2]/dbin)]+=1.;

        Real3D remol(0.0);
        int halfsize=static_cast<int>(polymers[i]->size()/2);
//        bool to_large=false;
        for(int j=0;j<polymers[i]->size()-1;j++){
            Real3D dr=pbc.getMinimumImageVector(polymers[i]->at(j+1)->coord, polymers[i]->at(j)->coord);
            remol+=dr;
//            if(dr.abs()>4){
//                std::cout << dr << std::endl;
//                to_large=true;
//            }
            if(j==halfsize-1){
                real rssqrmol=remol.sqr();
                rssqr+=rssqrmol;
                rsxsqr+=remol[0]*remol[0];
                rsysqr+=remol[1]*remol[1];
                rszsqr+=remol[2]*remol[2];
                dist[2][int(rssqrmol/dbin)]+=1.;
                dist[9][int((remol[0]*remol[0])/dbin)]+=1.;
                dist[10][int((remol[1]*remol[1])/dbin)]+=1.;
                dist[11][int((remol[2]*remol[2])/dbin)]+=1.;
            }
        }
        real resqrmol=remol.sqr();
//        if(to_large)
//            std::cout << remol << resqrmol << std::endl;
        resqr+=resqrmol;
        rexsqr+=remol[0]*remol[0];
        reysqr+=remol[1]*remol[1];
        rezsqr+=remol[2]*remol[2];
        dist[1][int(resqrmol/dbin)]+=1.;
        dist[6][int((remol[0]*remol[0])/dbin)]+=1.;
        dist[7][int((remol[1]*remol[1])/dbin)]+=1.;
        dist[8][int((remol[2]*remol[2])/dbin)]+=1.;

    }
    rgsqr/=num_pol;
    resqr/=num_pol;
    rssqr/=num_pol;
    rgxsqr/=num_pol;
    rgysqr/=num_pol;
    rgzsqr/=num_pol;
    rexsqr/=num_pol;
    reysqr/=num_pol;
    rezsqr/=num_pol;
    rsxsqr/=num_pol;
    rsysqr/=num_pol;
    rszsqr/=num_pol;
    tevol[0][step]=rgsqr;
    tevol[1][step]=resqr;
    tevol[2][step]=rssqr;
    tevol[3][step]=rgxsqr;
    tevol[4][step]=rgysqr;
    tevol[5][step]=rgzsqr;
    tevol[6][step]=rexsqr;
    tevol[7][step]=reysqr;
    tevol[8][step]=rezsqr;
    tevol[9][step]=rsxsqr;
    tevol[10][step]=rsysqr;
    tevol[11][step]=rszsqr;
    return;
}







