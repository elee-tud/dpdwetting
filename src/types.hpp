#ifndef _TYPES_HPP
#define _TYPES_HPP
#include <vector>
#include <iostream>
#include <string>
#include <mpi.h>

namespace dpd{

#define PI 3.141592653589793
typedef double real;
typedef std::vector<real> Rvec;
typedef std::vector<Rvec> Rvec2D;
typedef std::vector<Rvec2D> Rvec3D;
typedef std::vector<int> Ivec;
typedef std::vector<Ivec> Ivec2D;
typedef std::vector<Ivec2D> Ivec3D;
typedef std::vector<int>::iterator IvecIterator;
typedef std::vector<bool> Bvec;
typedef std::vector<Bvec> Bvec2D;
typedef std::vector<std::string> Svec;





template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v){
    os << "[";
    for(int i=0;i<v.size();++i){
        os<<v[i];
        if(i!=v.size()-1)
            os << "," ;
    }
    os << "]";
    return os;

}

template <typename T> 
int sign(T val){
    return (T(0) < val)-(val < T(0));
}


#define RUN 0
#define DROPSIZE 1
#define SPHERICALSTRESS 2
#define VELOCITY 3
#define RADIALDENSITY 4
#define AVGSTRESS 5
#define POLADSORP 6
#define POLSIZE 7
#define POLEVRLX 8
#define POLORIENT 9
#define BONDLENGTH 10
#define POLSTRETCH 11
#define TRJTOGRO 13
#define POLSMSF 14
#define POLSUBSIZE 15
#define POLMSD 16
#define VELACF 17
#define SURFCOV 18
#define RDF 19
#define DEPORIENT 20
#define BRIDGESIZE 21
#define BRIDGEVEL 22
#define BRIDGEPCONC 23
#define BRIDGEADSCONC 24
#define BRIDGECLVEL 25
#define DROPZD 26
#define DROPSRDF 27
#define JUMPFREQ 28
#define BRIDGEJF 29
#define BRIDGESLVEL 30
#define BRIDGEVELX 31
#define BRIDGEVELXZ 32
#define BRIDGEALDENS 33
#define BRIDGEZD 34
#define BRIDGEDENSXZ 35

#define TRUEPTCL 0
#define GHOSTPTCL 1
#define NOTEXIST 2

#define MASTER 0
#define NECESSARY -2515312

#define VELVER 0
#define EMIN 1
#define SLLOD 2

#define NBDPD 1
#define NBMDPD 2

#define BLHARMONIC 1

#define NOWALL 0
#define SOLIDWALL 1

#define CRDVEL 0
#define STRESS 1
#define FORCE 2

#define NODEFORMATION 0
#define ELONGATION 1
#define DUMPPRECISION 5

#define LINEARSS 0
#define RINGSS 1
};
#endif


