#include "surfacerdf.hpp"
#include <algorithm>


using namespace dpd;

SurfaceRdf::SurfaceRdf(InitialSet initset):Property(initset){
    title= "  Calculation of RDF on the surface";
    
    nprops=4;
    for_distrib=true;
    need_position=true;

    outfile="rdf_surface.out";

    outtitle_dist="#RDF on surface";
    outheader_dist=Svec{"r", "g(r)", "g_SS(r)", "g_SP(r)", "g_PP(r)"};

}


void SurfaceRdf::getSpecificParameters(){
    dbin=0.02;
    dsz=1.0;
    command->getCommandSingleOption("-dr", dbin, &dbin);
    command->getCommandSingleOption("-dsz", dsz, &dsz);
    surfaceb=control->getWallMinPosition()[2];
    surfacet=control->getWallMaxPosition()[2];
    liquidgrps=control->getLiquidGroups();
    return;
}


void SurfaceRdf::initializeVariables(){
    ndbin=static_cast<int>(pow(box[0]*box[0]+box[1]*box[1], 0.5)/dbin)+1;
    if(liquidgrps.size()==0){
        for(int i=0;i<nbeads;i++){
            liquididx.push_back(i);
        }
    }
    else{
        for(int i=0;i<nbeads;i++){
            if(std::find(liquidgrps.begin(), liquidgrps.end(), particles[i]->getMoleculeName())!=liquidgrps.end())
                liquididx.push_back(i);
        }
    }
    nliqptcls=liquididx.size();

    nref=0;
    nref_s=0;
    nref_p=0;
    surfdens=0.;
    boxcx=Rvec{box[0]/2-0.5, box[0]/2+0.5};
    boxcy=Rvec{box[1]/2-0.5, box[1]/2+0.5};
    initializeResultArrays();
    return;

}

void SurfaceRdf::calculateStep(int step){
    for(int i=0;i<nliqptcls;i++){
        if(particles[liquididx[i]]->coord[2]>=surfaceb && particles[liquididx[i]]->coord[2]<surfaceb+dsz){
            if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0 )
                nref_s++;
            else if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0 )
                nref_p++;
            nref++;
            if(particles[liquididx[i]]->coord[0]>=boxcx[0] && particles[liquididx[i]]->coord[0]<boxcx[1]
                    && particles[liquididx[i]]->coord[1]>=boxcy[0] && particles[liquididx[i]]->coord[1]<boxcy[1])
                surfdens+=1.;
            for(int j=i+1;j<nliqptcls;j++){
                if(particles[liquididx[j]]->coord[2]>=surfaceb && particles[liquididx[j]]->coord[2]<surfaceb+dsz){
//                    Real3D r=pbc.getMinimumImageVector(particles[liquididx[i]]->coord, particles[liquididx[j]]->coord);
                    Real3D r=particles[liquididx[i]]->coord-particles[liquididx[j]]->coord;
                    int index=static_cast<int>(sqrt(r[0]*r[0]+r[1]*r[1])/dbin);
                    if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0 && particles[liquididx[j]]->getMoleculeName().compare("SOL")==0 )
                        dist[1][index]+=1;
                    else if(particles[liquididx[i]]->getMoleculeName().compare("SOL")==0 && particles[liquididx[j]]->getMoleculeName().compare("POL")==0 )
                        dist[2][index]+=1;
                    else if(particles[liquididx[i]]->getMoleculeName().compare("POL")==0 && particles[liquididx[j]]->getMoleculeName().compare("POL")==0 )
                        dist[3][index]+=1;
                dist[0][index]+=1.;
                }
            }
        }
    }


    return;
}

void SurfaceRdf::normalizeResults(){
    std::cout << "surface_density=" << surfdens/nsteps << std::endl;
    for(int j=0;j<ndbin;j++){
        dist_sum[0][j]/=nref*PI*(j+0.5)*dbin*dbin;
        dist_sum[1][j]/=nref_s*PI*(j+0.5)*dbin*dbin;
        dist_sum[2][j]/=nref_s*PI*(j+0.5)*dbin*dbin;
        dist_sum[3][j]/=nref_p*PI*(j+0.5)*dbin*dbin;
    }
    return;
}




