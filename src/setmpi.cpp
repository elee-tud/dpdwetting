#include "setmpi.hpp"
#include <iostream>
#include <fstream>


using namespace dpd;

SetMPI::SetMPI(int argc, char* argv[]){
	int rc=MPI_Init(&argc, &argv);
	rc=MPI_Comm_size(MPI_COMM_WORLD, &_size);
	rc=MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
}
SetMPI::~SetMPI(){
}

void SetMPI::printRank(){
    std::cout << "processor id=" << _rank << std::endl;
}

void SetMPI::printSize(){
    std::cout << "total processors=" << _size << std::endl;
}

void SetMPI::finalize(){
	MPI_Finalize();
}


