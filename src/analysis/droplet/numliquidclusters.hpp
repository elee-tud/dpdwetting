#ifndef __NUMLIQUIDCLUSTERS__HPP
#define __NUMLIQUIDCLUSTERS__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
#include "../../indexing.hpp"
namespace dpd{

class NumberOfLiquidClusters:public Property{
private:
    Ivec liquididx;
    Svec liquidgrps;
    int nliqptcls;
  
    real distcrit;
    int numminptcls;
    int max_neighbors;

    Indexing index;
    Ivec2D neighborcells;
    Ivec2D nbsearch;
    int nnbcells;

    Ivec2D ptcls_in_cell;
    Ivec cellidx;

    Ivec ncells;
    int ntotcells;






    


public:
    NumberOfLiquidClusters(){}
    NumberOfLiquidClusters(InitialSet initset);
    ~NumberOfLiquidClusters(){}

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);

    void findCellIndex();
};
};



#endif
