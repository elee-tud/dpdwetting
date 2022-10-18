/*************************************************************
 *           Class "SetMPI"
 *************************************************************
 *  This class contains functions and variables for controlling
 *  MPI routines.
 *************************************************************/

#ifndef __SETMPI_HPP
#define __SETMPI_HPP
#include <mpi.h>
#include "types.hpp"


namespace dpd{

class SetMPI{
private:
	int _size;      /*Size of MPI*/
	int _rank;      /*Current process ID*/
public:
	SetMPI(){}      
	SetMPI(int argc, char* argv[]);
	~SetMPI();
	void printRank();       /*Printing out the current process ID*/
	void printSize();       /*Printing out the total number of processors*/
	int size(){	return _size;}  
	int rank(){ return _rank;}
	void finalize();
    bool isMaster(){return _rank==MASTER;}      /*Checking if it's master processor*/

};
};
#endif
