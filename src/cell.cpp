#include "cell.hpp" 
#include "indexing.hpp"
#include <algorithm>

using namespace dpd;

Cell::Cell(int myprocid, Ivec num_domains, int mycellid, Ivec num_cells):myprocid(myprocid), num_domains(num_domains), mycellid(mycellid), num_cells(num_cells){
    indexing=Indexing(num_cells);
    indexing_domain=Indexing(num_domains);
    totnum_cells=num_cells[0]*num_cells[1]*num_cells[2];
    num_original_cells=Ivec(3,0);
    num_original_cells[0]=num_cells[0]-2;
    num_original_cells[1]=num_cells[1]-2;
    num_original_cells[2]=num_cells[2]-2;
    mycell3did=indexing.get3DIndexFromIndex(mycellid);  //indexing class of the cell
    nbsearch=indexing.buildNeighborSearchVector();      //indexing for neighbor search

    myproc3did=indexing_domain.get3DIndexFromIndex(myprocid);
    mybeads.reserve(100);
    mybeads.clear();
    num_mybeads=0;
    checkCellType();
    if(is_ghost_cell){
        findRealCellOfGhost();
    }
    else{
        findNeighborCells();
    }
    if(is_real_cell)
        findGhostOfRealCell();


}



    
void Cell::clear(){
    mybeads.clear();
    num_mybeads=0;
    return;
}

void Cell::addBeads(int index){
    mybeads.push_back(index);
    num_mybeads++;
    return;
}

void Cell::addBeads(Ivec indices){
    for(int i=0;i<indices.size();i++){
        mybeads.push_back(indices[i]);
    }
    num_mybeads+=indices.size();
    return;
}

void Cell::removeBeads(int index){
    mybeads.erase(std::remove(mybeads.begin(), mybeads.end(), index), mybeads.end());
    num_mybeads--;
    return;
}

void Cell::removeBeads(Ivec indices){
    for(int i=0;i<indices.size();i++)
        mybeads.erase(std::remove(mybeads.begin(), mybeads.end(), indices[i]), mybeads.end());
    num_mybeads-=indices.size();;
    return;
}

void Cell::printBeads(){
    std::cout << "Procid:"<< myproc3did << "(" << myprocid << "), Cell ID" << mycell3did <<"(" << mycellid<< "): =>" ;
    if(isGhostCell())
        std::cout << "Ghost of Proc. " << realcell_index[0] << " Cell " << realcell_index[1] ;
    else
        std::cout << "Normal Cell" ;
    std::cout << "===>Beads= " << mybeads << std::endl;
    return;
}

void Cell::checkCellType(){
    checkCellType(num_cells);
    return;
}

void Cell::checkCellType(Ivec _num_cells){
    is_ghost_cell=false;
    is_real_cell=false;
    if(mycell3did[0]==0 || mycell3did[1]==0 || mycell3did[2]==0 || mycell3did[0]==_num_cells[0]-1 || mycell3did[1]==_num_cells[1]-1 || mycell3did[2]==_num_cells[2]-1){     //if the cell is at the first laye from the border of a domain
        is_ghost_cell=true;
    }
    else{
        if(mycell3did[0]==1 || mycell3did[1]==1 || mycell3did[2]==1 || mycell3did[0]==_num_cells[0]-2 || mycell3did[1]==_num_cells[1]-2 || mycell3did[2]==_num_cells[2]-2)  //if the cell is at the second layer from the border of a domain
            is_real_cell=true;
    }
    return;
}


void Cell::findNeighborCells(){
    neighbor_cells.clear();
    Ivec neighbor(3,0);
    for(int i=0;i<nbsearch.size();i++){
        neighbor=indexing.addIndexToIndex(mycell3did, nbsearch[i]);
        neighbor_cells.push_back(indexing.getIndexFrom3DIndex(neighbor));

    }
    return;
}

void Cell::findRealCellOfGhost(){
    realcell_index=Ivec(2,0);
    Ivec delta_proc(3,0);
    Ivec realcellid(3,0);
    for(int i=0;i<3;i++){
        if(mycell3did[i]==0){
            delta_proc[i]=-1;
            realcellid[i]=num_cells[i]-2;
        }
        else if(mycell3did[i]==num_cells[i]-1){
            delta_proc[i]=1;
            realcellid[i]=1;
        }
        else
            realcellid[i]=mycell3did[i];
    }
    realcell_index[0]=indexing_domain.addIndexToIndex(myprocid, delta_proc);    //Domain ID which contains its real cell
    realcell_index[1]=indexing.getIndexFrom3DIndex(realcellid);         //Cell ID of the real cell for this ghost cell
    return;
}


void Cell::findGhostOfRealCell(){
    Ivec dx{0};
    Ivec dy{0};
    Ivec dz{0};
    ghost_index.clear();
    Ivec2D nbsearch_back;
    if(mycell3did[0]==1)
        dx.push_back(-1);
    if(mycell3did[0]==num_cells[0]-2)
        dx.push_back(1);
    if(mycell3did[1]==1)
        dy.push_back(-1);
    if(mycell3did[1]==num_cells[1]-2) dy.push_back(1);
    if(mycell3did[2]==1)
        dz.push_back(-1);
    if(mycell3did[2]==num_cells[2]-2)
        dz.push_back(1);
    for(int ix=0;ix<dx.size();ix++){
        for(int iy=0;iy<dy.size();iy++){
            for(int iz=0;iz<dz.size();iz++){
                if(ix!=0 || iy!=0 || iz!=0){
                    nbsearch_back.push_back(Ivec{dx[ix], dy[iy], dz[iz]});
                }
            }
        }
    }
    for(int i=0;i<nbsearch_back.size();i++){
        int ghostproc=indexing_domain.addIndexToIndex(myprocid, nbsearch_back[i]);   //Domain index containing the ghost cell of the current real cell
        Ivec ghostid=mycell3did;
        for(int j=0;j<3;j++){
            if(nbsearch_back[i][j]==1)
                ghostid[j]=0;
            else if(nbsearch_back[i][j]==-1)
                ghostid[j]=num_cells[j]-1;
        }
        int ghostcell=indexing.getIndexFrom3DIndex(ghostid);  //Cell index of the ghost cell for the current real cell
        ghost_index.push_back(Ivec{ghostproc, ghostcell}); //A set of ghost domain/cell information
    }
    return;
}




    



void Cell::resetCell(int _mycellid, Ivec _num_cells){
    mycellid=_mycellid;
    num_cells=_num_cells;
    mybeads.reserve(100);
    clear();
    indexing.reset(num_cells);
    totnum_cells=num_cells[0]*num_cells[1]*num_cells[2];
    num_original_cells[0]=num_cells[0]-2;
    num_original_cells[1]=num_cells[1]-2;
    num_original_cells[2]=num_cells[2]-2;
    mycell3did=indexing.get3DIndexFromIndex(mycellid);
    nbsearch=indexing.buildNeighborSearchVector();
    checkCellType();
    if(is_ghost_cell){
        findRealCellOfGhost();
    }
    else{
        findNeighborCells();
    }
    if(is_real_cell)
        findGhostOfRealCell();

    return;

}

    
