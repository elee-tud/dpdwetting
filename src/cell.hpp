#ifndef __CELL__HPP
#define __CELL__HPP


/****************************************************************************
 * Class Cell
 *
 * The class to define cells in a domain
 ****************************************************************************/

#include "types.hpp"
#include "real3d.hpp"
#include "indexing.hpp"

namespace dpd{

class Cell{
private:
    Indexing indexing;      //Cell indexing
    Ivec num_cells;         //Number of original+ghost cells in domain,
    Indexing indexing_domain;       //Domain indexing
    Ivec num_domains;               //Number of domains
    Ivec num_original_cells;        //Number of original cells
    int totnum_cells;               //Number of all cells
    int myprocid;                   //My domain index
    Ivec myproc3did;                //My domain 3d index
    int mycellid;                   //My cell index
    Ivec mycell3did;                //My cell 3d index 

    bool is_ghost_cell;
    bool is_real_cell;


    int num_mybeads;                //Number of beads in a domain
    Ivec mybeads;                   //Indices of beads in a domain
    bool is_inner_cell;



    Ivec neighbor_cells;            //Indices of neighboring cells
    Ivec2D nbsearch;                //Neighbor searching index
    void findRealCellOfGhost();     //Finding a real cell of neighboring domain corresponding to the ghost cell of the current domain
    void findGhostOfRealCell();     //Finding a ghost cell of neighboring domain corresponding to the real cell of the current domain 
    void findNeighborCells();       //Finding neighbor cells
    
    Ivec realcell_index;            //Indices of the real cell
    Ivec2D ghost_index;             //Indices of ghost cells
    


public:
    Cell(){}
    Cell(int myprocid, Ivec num_domains, int mycellid, Ivec num_cells);
    ~Cell(){}
    void checkCellType();                   //Checking the cell type by 1D index, real or ghost?
    void checkCellType(Ivec _num_cells);    //Checking the cell type by 3D index
    bool isGhostCell() {return is_ghost_cell;}      //is ghost cell?
    bool isRealCell() {return is_real_cell; }       //is real cell (of a certain ghost cell) ?
    bool isTrueCell() {return !is_ghost_cell; }     //is a true cell (not a ghost cell)?
    Ivec& getBeads(){ return mybeads; }     //Returning indices in the cell
    void clear();                   //Clearing up all the beads the cell
    void addBeads(int index);       //Adding a bead to the cell
    void addBeads(Ivec indices);    //Adding beads to the cell
    void removeBeads(int index);    //Removing a bead from the cell
    void removeBeads(Ivec indices); //Removing beads from the cell


    Ivec& getMyCell3DIndex(){ return mycell3did; }  //Returning cell 3D index
    int& getMyCellIndex(){ return mycellid; }       //Returning cell 1D index
    void printBeads();      //Printing out the information of beads in the cell
    Ivec& getRealCellIndex() {return realcell_index; }      //Returning the real cell indices
    Ivec& getNeighborCells() {return neighbor_cells; }      //Returning the neighboring cell indices
    Ivec2D getGhostIndex() {return ghost_index; }           //Returning the ghost cell indices

    
    void resetCell(int mycellid, Ivec num_cells);       //Resetting the cell
};

typedef std::vector<Cell*> CellList;
};        
#endif
