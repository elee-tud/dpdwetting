#include "restart.hpp"
#include "filecontrol.hpp"

using namespace dpd;

Restart::Restart(Command* command, Configuration* config, Topology* topol, Decomposition* decomp, SetMPI* mpi):command(command), config(config), topol(topol), decomp(decomp), mpi(mpi) {
    ckp_fname=command->checkpoint();
    if(mpi->isMaster()){
        openFileWithBackup(ckp_fname, &ckpstream, true, false);
        dumpmult=pow(10,DUMPPRECISION);
        ckpstream.close();
    }
    particles=config->getParticles();


}

void Restart::writeCheckPoint(int step, real t){
    if(mpi->isMaster()){
        ckpstream.open(ckp_fname, std::ios::out | std::ios::binary);
        box=config->getBox();
        ckpstream.write((char*)&step, sizeof(int));
        ckpstream.write((char*)&t, sizeof(real));
        int ntotparticles=particles.size();
        ckpstream.write((char*)&ntotparticles, sizeof(int));
        for(int i=0;i<particles.size();i++){
            int modval=particles[i]->getCoord()[0]*dumpmult;
            ckpstream.write((char*)&modval, sizeof(int));
            modval=particles[i]->getCoord()[1]*dumpmult;
            ckpstream.write((char*)&modval, sizeof(int));
            modval=particles[i]->getCoord()[2]*dumpmult;
            ckpstream.write((char*)&modval, sizeof(int));
            modval=particles[i]->getVeloc()[0]*dumpmult;
            ckpstream.write((char*)&modval, sizeof(int));
            modval=particles[i]->getVeloc()[1]*dumpmult;
            ckpstream.write((char*)&modval, sizeof(int));
            modval=particles[i]->getVeloc()[2]*dumpmult;
            ckpstream.write((char*)&modval, sizeof(int));
        }
        int modval=box[0]*dumpmult;
        ckpstream.write((char*)&modval, sizeof(int));
        modval=box[1]*dumpmult;
        ckpstream.write((char*)&modval, sizeof(int));
        modval=box[2]*dumpmult;
        ckpstream.write((char*)&modval, sizeof(int));
        ckpstream.close();
    }
    return;

}
    


void Restart::readCheckPoint(int* step, real* t){
    if(mpi->isMaster()){
        rst_fname=command->restart();
        try{
            rststream.open(rst_fname, std::ios::in | std::ios::binary);
            if(!rststream.is_open())
                throw rst_fname;
        }catch(std::string traj_fname){
            std::cout << "File name " << rst_fname << " does not exist." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 11);
            exit(0);
        }

        int readval;
        int nbeads;
        rststream.read((char*)step, sizeof(int));
        rststream.read((char*)t, sizeof(real));
        rststream.read((char*)&nbeads, sizeof(int));
        for(int i=0;i<nbeads;i++){
            Real3D newcoord{0,0,0};
            rststream.read((char*)&readval, sizeof(int));
            newcoord[0]=(real)readval/dumpmult;
            rststream.read((char*)&readval, sizeof(int));
            newcoord[1]=(real)readval/dumpmult;
            rststream.read((char*)&readval, sizeof(int));
            newcoord[2]=(real)readval/dumpmult;
            particles[i]->setCoord(newcoord);
            Real3D newveloc{0,0,0};
            rststream.read((char*)&readval, sizeof(int));
            newveloc[0]=(real)readval/dumpmult;
            rststream.read((char*)&readval, sizeof(int));
            newveloc[1]=(real)readval/dumpmult;
            rststream.read((char*)&readval, sizeof(int));
            newveloc[2]=(real)readval/dumpmult;
            particles[i]->setVeloc(newveloc);
        }
        rststream.read((char*)&readval, sizeof(int));
        box[0]=(real)readval/dumpmult;
        rststream.read((char*)&readval, sizeof(int));
        box[1]=(real)readval/dumpmult;
        rststream.read((char*)&readval, sizeof(int));
        box[2]=(real)readval/dumpmult;
    }
    config->setBox(box);
    return;
}
