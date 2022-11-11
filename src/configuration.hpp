#ifndef __CONFIGURATION__HPP
#define __CONFIGURATION__HPP

/****************************************************************************
 * Class Configuration
 *
 * This class is to define the configuration of the system and to read and 
 * to write it from/to a file. You can also broadcast configuration to all 
 * other processors.
****************************************************************************/

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

    std::string config_fname;           /*File name of the configuration*/
    std::ifstream confstream;           /*Input stream of the configuration*/
    std::string title_conf;             /*Title of the gro file*/
    int num_beads_conf;                 /*The number of beads in a gro file*/
   
    int num_beads;                      /*The number of total beads*/
    int num_bonds;                      /*the number of bonds*/

    Svec molname;                       /*List of molecule names*/
    Svec beadtype;                      /*List of particle types*/
    ParticleList particles;             /*List of particle pointers*/
    Real3D box;                         /*Box size*/
    
    void openFile();                    /*Opening a gro file*/
    void closeFile(){ confstream.close(); }     /*Closing the gro file*/


public:
    Configuration(){}
    Configuration(Topology* topology, SetMPI* mpi);
    ~Configuration(){}
    void readConfiguration();           /*Reading configuration*/
//    int readConfigurationWithoutErrorMismatch();        /*Reading  configuration ignoring mismatch of the particle numbers*/
    void bcastConfiguration();          /*Broadcasting configuration*/

    Real3D& getBox(){ return box; }     /*Returning box size*/
    void setBox(Real3D _box){ box=_box; }       /*Setting the box size*/
    Topology* getTopology(){ return topology; }     /*Returning topology pointer*/
    ParticleList& getParticles(){ return particles; }       /*Returning particle list*/
    std::string getTitle(){ return title_conf; }            /*Returning the title of the gro file*/

    void writeConfiguration(std::string fname, int step, real t);       /*Writing configuration to a certain file*/


};
};
#endif

