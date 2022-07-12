#ifndef __INDEXING__HPP
#define __INDEXING__HPP

#include "types.hpp"

namespace dpd{

class Indexing{
private:
    Ivec size;
public:
    Indexing(){}
    Indexing(Ivec& size):size(size){}
    ~Indexing(){}

    void reset(Ivec& _size);
    int getIndexFrom3DIndex(Ivec& index3d);
    Ivec get3DIndexFromIndex(int index);
    Ivec getIndexInBox(Ivec& index);
    Ivec addIndexToIndex(Ivec& index1, Ivec& index2);
    int addIndexToIndex(int index1, int index2);
    int addIndexToIndex(int index1, Ivec& index2);
    Ivec subractIndexFromIndex(Ivec& index1, Ivec& index2);
    Ivec2D buildNeighborSearchVector();
};
};

#endif
