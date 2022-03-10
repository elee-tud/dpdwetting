#include "configuration.hpp"
#include <sstream>
#include <iomanip>
#include <iostream>

using namespace dpd;
Configuration::Configuration(Topology* topology, SetMPI* mpi):command(command), topology(topology), mpi(mpi){

    command=topology->getCommand();
    config_fname=command->configuration();
    num_beads=topology->getNbeads();
    num_bonds=topology->getNbonds();
    particles=topology->getParticles();
    mpiptcl=MPIClasses(mpi);
    err=Error(mpi);
    box=Real3D(0.);


}


void Configuration::openFile(){
	try{
		confstream.open(config_fname);
		if(!confstream.is_open())
			throw config_fname;
	}catch(std::string config_fname){
    std::cout << "File name " << config_fname << "does not exist." << std::endl;
		exit(0);
	}
    return;
}
int Configuration::readConfigurationWithoutErrorMismatch(){
    if(mpi->isMaster()){
        openFile();
        std::string line;
        std::getline(confstream, line);
        title_conf=line;
        std::getline(confstream, line);
        num_beads_conf=std::stoi(line);
        for(int i=0;i<num_beads;i++){
            std::getline(confstream, line);
            real x=std::stod(line.substr(20,8));
            real y=std::stod(line.substr(28,8));
            real z=std::stod(line.substr(36,8));
            Real3D newcoord(x, y, z);
            particles[i]->setCoord(newcoord);
            if(line.length()>44){
                real vx=std::stod(line.substr(44,8));
                real vy=std::stod(line.substr(52,8));
                real vz=std::stod(line.substr(60,8));
                Real3D newveloc(vx, vy, vz);
                particles[i]->setVeloc(newveloc);
            }
        }
        confstream >> box[0] >> box[1] >> box[2];
        closeFile();
    }
    return num_beads_conf;
}
void Configuration::readConfiguration(){
    if(mpi->isMaster()){
        openFile();
        std::string line;
        std::getline(confstream, line);
        title_conf=line;
        std::getline(confstream, line);
        num_beads_conf=std::stoi(line);
        err.mismatchedNumberOfParticles(num_beads_conf, num_beads);
        for(int i=0;i<num_beads;i++){
            std::getline(confstream, line);
            err.mismatchMoleculeName(line, particles[i]);
            real x=std::stod(line.substr(20,8));
            real y=std::stod(line.substr(28,8));
            real z=std::stod(line.substr(36,8));
            Real3D newcoord(x, y, z);
            particles[i]->setCoord(newcoord);
            if(line.length()>44){
                real vx=std::stod(line.substr(44,8));
                real vy=std::stod(line.substr(52,8));
                real vz=std::stod(line.substr(60,8));
                Real3D newveloc(vx, vy, vz);
                particles[i]->setVeloc(newveloc);
            }
        }
        confstream >> box[0] >> box[1] >> box[2];
        closeFile();
    }
    return;
}


void Configuration::bcastConfiguration(){
    Ivec allbeads(num_beads);
    for(int i=0;i<num_beads;i++){
        allbeads[i]=i;
    }
    Rvec commcoord(3*num_beads);
    Rvec commveloc(3*num_beads);
    if(mpi->isMaster()){
        commcoord=serializeNCoords(allbeads, particles);
        commveloc=serializeNVelocs(allbeads, particles);
    }
    MPI_Bcast(&commcoord[0], 3*num_beads, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&commveloc[0], 3*num_beads, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    if(!mpi->isMaster()){
        deserializeNCoords(allbeads, particles, commcoord);
        deserializeNVelocs(allbeads, particles, commveloc);
    }
    MPI_Bcast(&box[0], 3, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    return;
}

void Configuration::writeConfiguration(std::string fname, int step, real t){
    if(mpi->isMaster()){
        std::ofstream outconfstream;
        outconfstream.open(fname);
        std::stringstream title;
        title << title_conf << " at step=" << step <<  ", time=" << t;
        outconfstream << title.str() << std::endl;
        outconfstream << particles.size() << std::endl;
        for(int i=0;i<particles.size();i++){
            outconfstream << std::setw(5) << particles[i]->getMoleculeIndex();
            outconfstream << std::left << std::setw(5) << particles[i]->getMoleculeName();
            outconfstream << std::right << std::setw(5) << particles[i]->getTypeName();
            outconfstream << std::setw(5) << particles[i]->getParticleIndex();
            outconfstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getCoord()[0];
            outconfstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getCoord()[1];
            outconfstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getCoord()[2];
            outconfstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getVeloc()[0];
            outconfstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getVeloc()[1];
            outconfstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getVeloc()[2];
            outconfstream << std::endl;
        }
        outconfstream << std::setw(10) << std::fixed << std::setprecision(5) << box[0];
        outconfstream << std::setw(10) << std::fixed << std::setprecision(5) << box[1];
        outconfstream << std::setw(10) << std::fixed << std::setprecision(5) << box[2] << std::endl;
        outconfstream.close();
    }
    return;
}


    
