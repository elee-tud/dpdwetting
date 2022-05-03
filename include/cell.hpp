#ifndef __CELL__HPP
#define __CELL__HPP
#include "types.hpp"
#include "real3d.hpp"
#include "indexing.hpp"

namespace dpd{

class Cell{
private:
    Indexing indexing;
    Ivec num_cells;         //Number of original+ghost cells in domain,
    Indexing indexing_domain;
    Ivec num_domains;
    Ivec num_original_cells;
    int totnum_cells;
    int myprocid;
    Ivec myproc3did;
    int mycellid;
    Ivec mycell3did;

    bool is_ghost_cell;
    bool is_real_cell;


    int num_mybeads;
    Ivec mybeads;
    bool is_inner_cell;



    Ivec neighbor_cells;
    Ivec2D nbsearch;
    void findRealCellOfGhost();
    void findGhostOfRealCell();
    void findNeighborCells();
    
    Ivec realcell_index;
    Ivec2D ghost_index;
    


public:
    Cell(){}
    Cell(int myprocid, Ivec num_domains, int mycellid, Ivec num_cells);
    ~Cell(){}
    void checkCellType();
    void checkCellType(Ivec _num_cells);
    bool isGhostCell() {return is_ghost_cell;}
    bool isRealCell() {return is_real_cell; }
    bool isTrueCell() {return !is_ghost_cell; }
    Ivec& getBeads(){ return mybeads; }
    void clear();
    void addBeads(int index);
    void addBeads(Ivec indices);
    void removeBeads(int index);
    void removeBeads(Ivec indices);


    Ivec& getMyCell3DIndex(){ return mycell3did; } 
    int& getMyCellIndex(){ return mycellid; }
    void printBeads();
    Ivec& getRealCellIndex() {return realcell_index; }
    Ivec& getNeighborCells() {return neighbor_cells; }
    Ivec2D getGhostIndex() {return ghost_index; }

    
    void resetCell(int mycellid, Ivec num_cells);
};

typedef std::vector<Cell*> CellList;
};        
#endif
