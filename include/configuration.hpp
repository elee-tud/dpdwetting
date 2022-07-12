#ifndef __CONFIGURATION__HPP
#define __CONFIGURATION__HPP


#include <string>
#include "types.hpp"
#include "topology.hpp"
#include "setmpi.hpp"
#include "command.hpp"
#include "real3d.hpp"
#include "errorhandling.hpp"
#include "mpiclasses.hpp"


namespace dpd{

class Configuration{
private:
    Command* command;
    Topology* topology;
    SetMPI* mpi;
    Error err;
    MPIClasses mpiptcl;

    std::string config_fname;
    std::ifstream confstream;
    std::string title_conf;
    int num_beads_conf;
   
    int num_beads;
    int num_bonds;

    Svec molname;
    Svec beadtype;
    ParticleList particles;
    Real3D box;
    
    void openFile();
    void closeFile(){ confstream.close(); }


public:
    Configuration(){}
    Configuration(Topology* topology, SetMPI* mpi);
    ~Configuration(){}
    void readConfiguration();
    int readConfigurationWithoutErrorMismatch();
    void bcastConfiguration();

    Real3D& getBox(){ return box; }
    void setBox(Real3D _box){ box=_box; }
    Topology* getTopology(){ return topology; }
    ParticleList& getParticles(){ return particles; }
    std::string getTitle(){ return title_conf; }

    void writeConfiguration(std::string fname, int step, real t);


};
};
#endif

