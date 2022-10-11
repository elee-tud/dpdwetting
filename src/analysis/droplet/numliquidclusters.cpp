#include "numliquidclusters.hpp"
#include <cmath>
#include <algorithm>

using namespace dpd;
NumberOfLiquidClusters::NumberOfLiquidClusters(InitialSet initset):Property(initset){
    title="  Calculation of the number of liquid clusters ";
    
    nprops=3;
    for_timeevol=true;
    need_position=true;

    outfile="clusters.out";
    
    outtitle_tevol="#The number and the size of liquid clusters";
    outheader_tevol=Svec{"Time", "N_cl", "<s_cl>", "s_large"};

}

void NumberOfLiquidClusters::getSpecificParameters(){
    distcrit=1.0;
    numminptcls=20;
    max_neighbors=100;
    command->getCommandSingleOption("-dc", distcrit, &distcrit);
    command->getCommandSingleOption("-mp", numminptcls, &numminptcls);
    command->getCommandSingleOption("-mb", max_neighbors, &max_neighbors);

    if( std::fmod(box[0],distcrit)>0.0001 || std::fmod(box[1],distcrit)>0.0001 || std::fmod(box[2],distcrit)>0.0001 ){
        std::cout << "Error. The criterion for distance of density calculation has to be a devisor of box sizes" << std::endl;
        exit(0);
    }



    liquidgrps=control->getLiquidGroups();
    

    return;
}

void NumberOfLiquidClusters::initializeVariables(){
    
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

    ncells=Ivec{static_cast<int>(box[0]/distcrit), static_cast<int>(box[1]/distcrit), static_cast<int>(box[2]/distcrit)};

    ntotcells=ncells[0]*ncells[1]*ncells[2];
    
    index=Indexing(ncells);
    nbsearch=index.buildNeighborSearchVector();
    nnbcells=nbsearch.size();
    neighborcells=Ivec2D(ntotcells, Ivec(nnbcells, -1));

    for(int i=0;i<ncells[0];i++){
        for(int j=0;j<ncells[1];j++){
            for(int k=0;k<ncells[2];k++){
                Ivec idx3d=Ivec{i,j,k};
                int idx=index.getIndexFrom3DIndex(idx3d);
                for(int l=0;l<nnbcells;l++){
                    int nbidx=index.addIndexToIndex(idx, nbsearch[l]);
                    nbidx=index.getIndexInBox(nbidx);
                    neighborcells[idx][l]=nbidx;
                }

            }
        }
    }



    initializeResultArrays();

    return;
}

void NumberOfLiquidClusters::findCellIndex(){
    cellidx=Ivec(nliqptcls, 0);
    ptcls_in_cell=Ivec2D(ntotcells, Ivec());
    for(int i=0;i<ntotcells;i++){
        ptcls_in_cell[i].reserve(max_neighbors);
    }
    for(int i=0;i<nliqptcls;i++){
        int pidx=liquididx[i];
        Real3D position=particles[pidx]->coord;
        Ivec cidx3d=Ivec{static_cast<int>(position[0]/distcrit), static_cast<int>(position[1]/distcrit), static_cast<int>(position[2]/distcrit)}; 
        int cidx=index.getIndexFrom3DIndex(cidx3d);
        cellidx[pidx]=cidx;
        ptcls_in_cell[cidx].push_back(pidx);

    }
    
}



void NumberOfLiquidClusters::calculateStep(int step){
}
