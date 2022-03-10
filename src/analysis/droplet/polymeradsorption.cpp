#include "polymeradsorption.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
PolymerAdsorption::PolymerAdsorption(InitialSet initset):Property(initset){
    title="   Calculation of # of polymer beads adsorbed on the surface";

    nprops=10;
    for_timeevol=true;
    need_position=true;

    outfile="polymer_adsorption.out";
    nbond_types=topol->getNbondTypes();
    if(nbond_types>0)
        bondl=topol->getBondL()[0];

    outtitle_tevol="#The number of polymer beads adsorbed on the surface";
    outheader_tevol=Svec{"Time", "n_ads", "m_ads", "<n_tail>", "<l_tail>", "<n_loop>", "<l_loop>", "<n_train>", "<l_train>", "<X_z>/<l>", "<m_l+ta>"};

}



void PolymerAdsorption::getSpecificParameters(){
    acrit=1.0;

    command->getCommandSingleOption("-ac", acrit, &acrit);
    surfaceb=control->getWallMinPosition()[2];
    surfacet=control->getWallMaxPosition()[2];
    return;
}


void PolymerAdsorption::initializeVariables(){
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

    initializeResultArrays();
    int lpol=polymers[0]->size();
    ads=Ivec(lpol, 0);
    return;
}

void PolymerAdsorption::calculateStep(int step){
    int mads=0;
    int ntrain=0, nloop=0, ntail=0;
    int iltrain=0, illoop=0, iltail=0;
    real ltrain=0., lloop=0., ltail=0.;
    real mtail=0.;
    bool istail=false;
    real stretch=0.;
    int nstretch=0;
    Rvec zcoords;
    zcoords.reserve(polymers[0]->size());
    for(int i=0;i<num_pol;i++){
        int lpol=polymers[i]->size();
        bool molads=false;
        for(int j=0;j<lpol;j++){
            if(polymers[i]->at(j)->coord[2]>=surfaceb && polymers[i]->at(j)->coord[2]<surfaceb+acrit){
                tevol[0][step]+=1.;
                ads[j]=1;
                molads=true;
            }
            else
                ads[j]=0;
        }
        if(molads){
            mads+=lpol;
            if(ads[0]==0){
                istail=true;
                iltail=0;
            }
            else
                iltrain=0;
        
            for(int j=0;j<lpol-1;j++){
                if(ads[j]==0 && ads[j+1]==0){
                    iltail++;
                    zcoords.push_back(polymers[i]->at(j)->coord[2]);
                    if(j==lpol-2){
                        ntail++;
                        iltail++;
                        zcoords.push_back(polymers[i]->at(j+1)->coord[2]);
                        ltail+=static_cast<real>(iltail);
                        mtail+=static_cast<real>(iltail*iltail);
                        if(iltail!=1){
                            nstretch+=iltail;
                            if(nbond_types>0)
                                stretch+=iltail*findStretchingFactor(zcoords);
                        }
                        zcoords.clear();
                        istail=false;

                    }
                }
                else if(ads[j]==0 && ads[j+1]==1){
                    if(istail){
                        ntail++;
                        iltail++;
                        zcoords.push_back(polymers[i]->at(j)->coord[2]);
                        ltail+=static_cast<real>(iltail);
                        mtail+=static_cast<real>(iltail*iltail);
                        if(iltail!=1){
                            nstretch+=iltail;
                            if(nbond_types>0)
                                stretch+=iltail*findStretchingFactor(zcoords);
                        }
                        zcoords.clear();

                        iltrain=0;
                        istail=false;
                        if(j==lpol-2){
                            ntrain++;
                            ltrain+=1.;
                        }
                            
                    }
                    else{
                        nloop++;
                        iltail++;
                        zcoords.push_back(polymers[i]->at(j)->coord[2]);
                        lloop+=static_cast<real>(iltail);
                        mtail+=static_cast<real>(iltail*iltail);
                        if(iltail!=1){
                            nstretch++;
                            if(nbond_types>0)
                                stretch+=iltail*findStretchingFactor(zcoords);
                        }
                        iltrain=0;
                        zcoords.clear();
                        if(j==lpol-2){
                            ntrain++;
                            ltrain+=1.;
                        }
                    }
                }
                else if(ads[j]==1 && ads[j+1]==1){
                    iltrain++;
                    if(j==lpol-2){
                        ntrain++;
                        ltrain+=static_cast<real>(iltrain);
                    }
                }
                else{
                    ntrain++;
                    iltrain++;
                    ltrain+=static_cast<real>(iltrain);
                    iltail=0;
                    if(j==lpol-2){
                        ntail++;
                        ltail+=1.0;
                    }
                }
            }
            
        }
        tevol[1][step]=mads;
        tevol[2][step]=ntail;
        tevol[3][step]=ltail;
        tevol[4][step]=nloop;
        tevol[5][step]=lloop;
        tevol[6][step]=ntrain;
        tevol[7][step]=ltrain;
        tevol[8][step]=stretch;
        tevol[9][step]=mtail;
    }
    tevol[9][step]/=tevol[3][step]+tevol[5][step];
    tevol[3][step]/=tevol[2][step];
    tevol[5][step]/=tevol[4][step];
    tevol[7][step]/=tevol[6][step];
    tevol[8][step]/=nstretch;


    return;
}

real PolymerAdsorption::findStretchingFactor(Rvec zcoords){
    real upmost=zcoords[0];
    real lowmost=zcoords[0];
    for(int i=1;i<zcoords.size();i++){
        if(zcoords[i]>upmost)
            upmost=zcoords[i];
        if(zcoords[i]<lowmost)
            lowmost=zcoords[i];
    }
    return (upmost-lowmost);
}



