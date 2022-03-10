#include "trajectory.hpp"
using namespace dpd;

Trajectory::Trajectory(InitialSet initset, int type):initset(initset), type(type){
    command=initset.command;
    control=initset.control;
    topol=initset.topol;
    config=initset.config;
    mpi=initset.mpi;
    decomp=initset.decomp;
    particles=config->getParticles();
    err=Error(mpi);
    dumpbinary=control->getDumpBinary();
    dumpmult=pow(10,DUMPPRECISION);
    if(type==CRDVEL)
        traj_fname=command->trajectory();
    else if(type==STRESS)
        traj_fname=command->stress();


    readFirstStep();




}

void Trajectory::readFirstStep(){
    openTrajectory();
    readTrajectoryStep();
    closeTrajectory();
    return;
}


void Trajectory::openTrajectory(){
    if(mpi->isMaster()){
        try{
            if(dumpbinary)
                trajstream.open(traj_fname, std::ios::in | std::ios::binary);
            else
                trajstream.open(traj_fname);
            if(!trajstream.is_open())
                throw traj_fname;
        }catch(std::string traj_fname){
            std::cout << "File name " << traj_fname << " does not exist." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 11);
            exit(0);
        }
    }
    return;
}

bool Trajectory::readTrajectoryStep(){
    if(mpi->isMaster()){
        if(dumpbinary){
            int readval;
            int step;
            float t;
            int nbeads;
            if(type==CRDVEL){
                trajstream.read((char*)&step, sizeof(int));
                trajstream.read((char*)&t, sizeof(float));
                trajstream.read((char*)&nbeads, sizeof(int));
                for(int i=0;i<nbeads;i++){
                    Real3D newcoord{0,0,0};
                    trajstream.read((char*)&readval, sizeof(int));
                    newcoord[0]=(real)readval/dumpmult;
                    trajstream.read((char*)&readval, sizeof(int));
                    newcoord[1]=(real)readval/dumpmult;
                    trajstream.read((char*)&readval, sizeof(int));
                    newcoord[2]=(real)readval/dumpmult;
                    particles[i]->setCoord(newcoord);
                    Real3D newveloc{0,0,0};
                    trajstream.read((char*)&readval, sizeof(int));
                    newveloc[0]=(real)readval/dumpmult;
                    trajstream.read((char*)&readval, sizeof(int));
                    newveloc[1]=(real)readval/dumpmult;
                    trajstream.read((char*)&readval, sizeof(int));
                    newveloc[2]=(real)readval/dumpmult;
                    particles[i]->setVeloc(newveloc);
                }
                trajstream.read((char*)&readval, sizeof(int));
                box[0]=(real)readval/dumpmult;
                trajstream.read((char*)&readval, sizeof(int));
                box[1]=(real)readval/dumpmult;
                trajstream.read((char*)&readval, sizeof(int));
                box[2]=(real)readval/dumpmult;
                config->setBox(box);
                decomp->resetBox(box);
            }
            else if(type==STRESS){
                trajstream.read((char*)&step, sizeof(int));
                trajstream.read((char*)&t, sizeof(float));
                trajstream.read((char*)&nbeads, sizeof(int));
                for(int i=0;i<nbeads;i++){
                    Rvec newstress(6, 0);
                    trajstream.read((char*)&readval, sizeof(int));
                    newstress[0]=(real)readval/dumpmult;
                    trajstream.read((char*)&readval, sizeof(int));
                    newstress[1]=(real)readval/dumpmult;
                    trajstream.read((char*)&readval, sizeof(int));
                    newstress[2]=(real)readval/dumpmult;
                    trajstream.read((char*)&readval, sizeof(int));
                    newstress[3]=(real)readval/dumpmult;
                    trajstream.read((char*)&readval, sizeof(int));
                    newstress[4]=(real)readval/dumpmult;
                    trajstream.read((char*)&readval, sizeof(int));
                    newstress[5]=(real)readval/dumpmult;
                    particles[i]->stress=newstress;
                }
            }

        }

        else{
            std::string line;
            if(type==CRDVEL){
                if(!std::getline(trajstream, line)) return false;
                title_traj=line;
                if(!std::getline(trajstream, line)) return false;
                nbeads=std::stoi(line);

                for(int i=0;i<nbeads;i++){
                    if(!std::getline(trajstream, line)) return false;
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
                if(!std::getline(trajstream, line)) return false;
                box[0]=std::stod(line.substr(0,10));
                box[1]=std::stod(line.substr(10,10));
                box[2]=std::stod(line.substr(20,10));
                config->setBox(box);
                decomp->resetBox(box);
            }
            else if(type==STRESS){
                if(!std::getline(trajstream, line)) return false;
                title_traj=line;
                if(!std::getline(trajstream, line)) return false;
                nbeads=std::stoi(line);

                for(int i=0;i<nbeads;i++){
                    if(!std::getline(trajstream, line)) return false;
                    real sxx=std::stod(line.substr(5,15));
                    real syy=std::stod(line.substr(20,15));
                    real szz=std::stod(line.substr(35,15));
                    real sxy=std::stod(line.substr(50,15));
                    real syz=std::stod(line.substr(65,15));
                    real szx=std::stod(line.substr(80,15));
                    particles[i]->stress=Rvec{sxx, syy, szz, sxy, syz, szx};
                }
            }
        }

    }
    if(trajstream.eof())
        return false;
    else
        return true;
}
bool Trajectory::skipTrajectoryStep(){
    if(mpi->isMaster()){
        if(dumpbinary){
            std::size_t skiplength;
            int readval;
            int step;
            float t;
            int nbeads;

            if(type==CRDVEL){
                trajstream.read((char*)&step, sizeof(int));
                trajstream.read((char*)&t, sizeof(float));
                trajstream.read((char*)&nbeads, sizeof(int));
                skiplength=(nbeads*6+3)*sizeof(int);
                trajstream.ignore(skiplength);
            }
            else if(type==STRESS){
                trajstream.read((char*)&step, sizeof(int));
                trajstream.read((char*)&t, sizeof(float));
                trajstream.read((char*)&nbeads, sizeof(int));
                skiplength=(nbeads*6)*sizeof(int);
                trajstream.ignore(skiplength);
            }

        }
        else{

            std::string line;
            int nline;
            if(type==CRDVEL)
                nline=nbeads+3;
            else if(type=STRESS)
                nline=nbeads+2;

            for(int i=0;i<nbeads+3;i++){
                if(!std::getline(trajstream, line)) return false;
            }
        }
    }
    if(trajstream.eof())
        return false;
    else
        return true;
}


void Trajectory::closeTrajectory(){ 
    if(mpi->isMaster()){
        trajstream.close(); 
    }
}

