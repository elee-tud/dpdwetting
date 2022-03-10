#include "indexing.hpp"

using namespace dpd;
void Indexing::reset(Ivec& _size){
    size=_size;
    return;
}
int Indexing::getIndexFrom3DIndex(Ivec& index3d){
    int idx=index3d[0]*size[1]*size[2]+index3d[1]*size[2]+index3d[2];
    return idx;
}

Ivec Indexing::get3DIndexFromIndex(int index){
    Ivec idx3d(3,0);
    idx3d[0]=index/(size[1]*size[2]);
    idx3d[1]=(index%(size[1]*size[2]))/size[2];
    idx3d[2]=(index%(size[1]*size[2]))%size[2];
    return idx3d;
}
Ivec Indexing::getIndexInBox(Ivec& index){
    Ivec result=index;
    if(result[0]<0)
        result[0]+=size[0];
    else if(result[0]>=size[0])
        result[0]-=size[0];
    if(result[1]<0)
        result[1]+=size[1];
    else if(result[1]>=size[1])
        result[1]-=size[1];
    if(result[2]<0)
        result[2]+=size[2];
    else if(result[2]>=size[2])
        result[2]-=size[2];
    return result;
}

    



Ivec Indexing::addIndexToIndex(Ivec& index1, Ivec& index2){
    Ivec result(3,0);
    result[0]=index1[0]+index2[0];
    result[1]=index1[1]+index2[1];
    result[2]=index1[2]+index2[2];
    result=getIndexInBox(result);
    return result;
}

int Indexing::addIndexToIndex(int idx1, int idx2){
    Ivec index1=get3DIndexFromIndex(idx1);
    Ivec index2=get3DIndexFromIndex(idx2);
    Ivec result=addIndexToIndex(index1, index2);
    return getIndexFrom3DIndex(result);
}
    

int Indexing::addIndexToIndex(int idx1, Ivec& index2){
    Ivec index1=get3DIndexFromIndex(idx1);
    Ivec result=addIndexToIndex(index1, index2);
    return getIndexFrom3DIndex(result);
}
    

    
    


Ivec Indexing::subractIndexFromIndex(Ivec& index1, Ivec& index2){ 
    Ivec result(3,0);
    result[0]=index1[0]-index2[0];
    result[1]=index1[1]-index2[1];
    result[2]=index1[2]-index2[2];
    result=getIndexInBox(result);
    return result;
}

Ivec2D Indexing::buildNeighborSearchVector(){
    Ivec2D nbsearch;
    nbsearch.push_back(Ivec{1,0,0});
    nbsearch.push_back(Ivec{0,1,0});
    nbsearch.push_back(Ivec{0,0,1});
    nbsearch.push_back(Ivec{1,1,0});
    nbsearch.push_back(Ivec{1,0,1});
    nbsearch.push_back(Ivec{0,1,1});
    nbsearch.push_back(Ivec{1,1,1});
    nbsearch.push_back(Ivec{1,0,-1});
    nbsearch.push_back(Ivec{1,-1,0});
    nbsearch.push_back(Ivec{0,1,-1});
    nbsearch.push_back(Ivec{1,1,-1});
    nbsearch.push_back(Ivec{1,-1,1});
    nbsearch.push_back(Ivec{1,-1,-1});
/*
    Ivec dx;
    Ivec dy;
    Ivec dz;
    if(size[0]>2)
        dx=Ivec{-1,0,1};
    else if(size[0]>1)
        dx=Ivec{0,1};
    else
        dx=Ivec{0};
    if(size[1]>2)
        dy=Ivec{-1, 0,1};
    else if(size[1]>1)
        dy=Ivec{0,1};
    else
        dy=Ivec{0};
    if(size[2]>2)
        dz=Ivec{-1,0,1};
    else if(size[2]>1)
        dz=Ivec{0,1};
    else
        dz=Ivec{0};
    nbsearch.clear();
    for(int ix=0;ix<dx.size();ix++){
        for(int iy=0;iy<dy.size();iy++){
            for(int iz=0;iz<dz.size();iz++){
                if(dx[ix]!=0 || dy[iy]!=0 || dz[iz]!=0){
                    nbsearch.push_back(Ivec{dx[ix], dy[iy], dz[iz]});
                }
            }
        }
    }
    */
    return nbsearch;
}



    


