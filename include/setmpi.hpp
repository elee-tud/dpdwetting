#ifndef __SETMPI_HPP
#define __SETMPI_HPP
#include <mpi.h>
#include "types.hpp"


namespace dpd{

class SetMPI{
private:
	int _size;
	int _rank;
public:
	SetMPI(){}
	SetMPI(int argc, char* argv[]);
	~SetMPI();
	void printRank();
	void printSize();
	int size(){	return _size;}
	int rank(){ return _rank;}
	void finalize();
    bool isMaster(){return _rank==MASTER;}

};
};
#endif
