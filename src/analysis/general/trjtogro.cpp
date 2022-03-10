#include "trjtogro.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>
#include "../../filecontrol.hpp"

using namespace dpd;
TrjToGro::TrjToGro(InitialSet initset):Property(initset){
    property=TRJTOGRO;
    title="  Converting <.trj>(BINARY) trajectory file to <.gro>(ASCII) file ";
    nprops=5;

    for_timeevol=true;
    need_position=true;
    need_velocity=true;

    if(mpi->isMaster()){
        outname=command->output_analysis();
        if(outname.compare("default")==0){
            std::string trjname=command->trajectory();
            outfile=trjname;
            outfile.replace(outfile.find(".trj"), std::string(".trj").length(), ".gro");
        }
        else
            outfile=outname;
    }
    dumpfrozen=control->getDumpFrozen();
    dpd::openFileWithBackup(outfile, &trajstream, false, false);
    nufbeads=topol->getNumUnfrozenBeads();


}

void TrjToGro::dumpGro(int step, real t){

    std::stringstream title;
    title << config->getTitle() << " at step=" << step <<  ", time=" << t;
    trajstream << title.str() << std::endl;
    if(!dumpfrozen)
        trajstream << nufbeads << std::endl;
    else
        trajstream << particles.size() << std::endl;
    for(int i=0;i<particles.size();i++){
        if(particles[i]->isFrozen() && !dumpfrozen){
        }
        else{
            trajstream << std::setw(5) << particles[i]->getMoleculeIndex()%100000;
            trajstream << std::left << std::setw(5) << particles[i]->getMoleculeName();
            trajstream << std::right << std::setw(5) << particles[i]->getTypeName();
            trajstream << std::setw(5) << particles[i]->getParticleIndex()%100000;
            trajstream << std::setprecision(3) << std::setw(8);
            trajstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getCoord()[0];
            trajstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getCoord()[1];
            trajstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getCoord()[2];
            trajstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getVeloc()[0];
            trajstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getVeloc()[1];
            trajstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getVeloc()[2];
            trajstream << std::endl;
        }
    }
    trajstream << std::setw(10) << std::fixed << std::setprecision(5) << box[0];
    trajstream << std::setw(10) << std::fixed << std::setprecision(5) << box[1];
    trajstream << std::setw(10) << std::fixed << std::setprecision(5) << box[2] << std::endl;
    return;
}

void TrjToGro::writeOutput(){
    std::cout << "Trajectory file, " << outfile << " is generated." << std::endl;
    trajstream.close();
}
