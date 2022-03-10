#include "dump.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
#include "filecontrol.hpp"

using namespace dpd;
Dump::Dump(Initialization* init, Analysis* analysis):init(init), analysis(analysis){
    
    config=init->getConfiguration();

    box=config->getBox();
    particles=config->getParticles();
    decomp=init->getDecomposition();
    mpi=decomp->getMPI();
    control=init->getControl();
    mpiptcl=MPIClasses(mpi);
    command=init->getCommand();
    dumpfrozen=control->getDumpFrozen();
    nufbeads=init->getTopology()->getNumUnfrozenBeads();
    dumpbinary=control->getDumpBinary();
    dumpmult=pow(10,DUMPPRECISION);
    restart=control->getRestart();
    if(mpi->isMaster()){
        traj_fname=command->trajectory();
        log_fname=command->log();
        force_fname=command->force();
        stress_fname=command->stress();
        if(!dumpbinary){
            if(traj_fname.find(".trj")!=std::string::npos)
                traj_fname.replace(traj_fname.find(".trj"), std::string(".trj").length(), ".gro");
            if(force_fname.find(".frc")!=std::string::npos)
                force_fname.replace(force_fname.find(".frc"), std::string(".frc").length(), ".fou");
            if(stress_fname.find(".str")!=std::string::npos)
                stress_fname.replace(stress_fname.find(".str"), std::string(".str").length(), ".sou");
        }
        dpd::openFileWithBackup(traj_fname, &trajstream, dumpbinary, restart);
        dpd::openFileWithBackup(log_fname, &logstream, false, restart);

        if(control->getStressFrequency()>0)
            dpd::openFileWithBackup(stress_fname, &stressstream, dumpbinary, restart);
        if(control->getForceFrequency()>0)
            dpd::openFileWithBackup(force_fname, &forcestream, dumpbinary, restart);
        if(!command->doRestart())
            writeLogHeader();
    }
}

Dump::~Dump(){
    if(mpi->isMaster()){
        trajstream.close();
        logstream.close();
        if(control->getStressFrequency()>0)
            stressstream.close();
    }
}






void Dump::dumpConf(bool on_screen, int step, real t){
    if(on_screen){
        decomp->gatherCoords();
        decomp->gatherVelocs();
        box=config->getBox();
        if(mpi->isMaster()){
            std::stringstream title;
            title << config->getTitle() << " at step=" << step <<  ", time=" << t;
            std::cout << title.str() << std::endl;
            if(!dumpfrozen)
                std::cout << nufbeads << std::endl;
            else
                std::cout << particles.size() << std::endl;
            for(int i=0;i<particles.size();i++){
                if(particles[i]->isFrozen() && !dumpfrozen){
                }
                else{
                    std::cout << std::setw(5) << particles[i]->getMoleculeIndex();
                    std::cout << std::left << std::setw(5) << particles[i]->getMoleculeName();
                    std::cout << std::right << std::setw(5) << particles[i]->getTypeName();
                    std::cout << std::setw(5) << particles[i]->getParticleIndex();
                    std::cout << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getCoord()[0];
                    std::cout << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getCoord()[1];
                    std::cout << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getCoord()[2];
                    std::cout << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getVeloc()[0];
                    std::cout << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getVeloc()[1];
                    std::cout << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getVeloc()[2];
                    std::cout << std::endl;
                }
            }
            std::cout << std::setw(10) << std::fixed << std::setprecision(5) << box[0];
            std::cout << std::setw(10) << std::fixed << std::setprecision(5) << box[1];
            std::cout << std::setw(10) << std::fixed << std::setprecision(5) << box[2] << std::endl;
        }
    }
    else
        dumpConf(step, t);
    return;
}


/*Dump final conformation(ASCII)*/
void Dump::dumpConf(int step, real t, std::string filename){
    decomp->gatherCoords();
    decomp->gatherVelocs();
    box=config->getBox();
    if(mpi->isMaster()){
        std::ofstream outstream;
        dpd::openFileWithBackup(command->finalconf(), &outstream, false, false);
        dumpfrozen=true;
        dumpbinary=false;
        std::stringstream title;
        title << config->getTitle() << " at step=" << step <<  ", time=" << t;
        outstream << title.str() << std::endl;
        if(!dumpfrozen)
            outstream << nufbeads << std::endl;
        else
            outstream << particles.size() << std::endl;
        for(int i=0;i<particles.size();i++){
            if(particles[i]->isFrozen() && !dumpfrozen){
            }
            else{
                outstream << std::setw(5) << particles[i]->getMoleculeIndex()%100000;
                outstream << std::left << std::setw(5) << particles[i]->getMoleculeName();
                outstream << std::right << std::setw(5) << particles[i]->getTypeName();
                outstream << std::setw(5) << particles[i]->getParticleIndex()%100000;
                outstream << std::setprecision(3) << std::setw(8);
                outstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getCoord()[0];
                outstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getCoord()[1];
                outstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getCoord()[2];
                outstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getVeloc()[0];
                outstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getVeloc()[1];
                outstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getVeloc()[2];
                outstream << std::endl;
            }
        }
        outstream << std::setw(10) << std::fixed << std::setprecision(5) << box[0];
        outstream << std::setw(10) << std::fixed << std::setprecision(5) << box[1];
        outstream << std::setw(10) << std::fixed << std::setprecision(5) << box[2] << std::endl;
        outstream.close();
    }




    return;
}


void Dump::dumpConf(int step, real t){
    decomp->gatherCoords();
    decomp->gatherVelocs();
    box=config->getBox();
    if(mpi->isMaster()){
        if(dumpbinary){
            trajstream.write((char*)&step, sizeof(int));
            trajstream.write((char*)&t, sizeof(float));
            if(!dumpfrozen)
                trajstream.write((char*)&nufbeads, sizeof(int));
            else{
                int ntotparticles=particles.size();
                trajstream.write((char*)&ntotparticles, sizeof(int));
            }
            for(int i=0;i<particles.size();i++){
                if(particles[i]->isFrozen() && !dumpfrozen){}
                else{
                    int modval=particles[i]->getCoord()[0]*dumpmult;
                    trajstream.write((char*)&modval, sizeof(int));
                    modval=particles[i]->getCoord()[1]*dumpmult;
                    trajstream.write((char*)&modval, sizeof(int));
                    modval=particles[i]->getCoord()[2]*dumpmult;
                    trajstream.write((char*)&modval, sizeof(int));
                    modval=particles[i]->getVeloc()[0]*dumpmult;
                    trajstream.write((char*)&modval, sizeof(int));
                    modval=particles[i]->getVeloc()[1]*dumpmult;
                    trajstream.write((char*)&modval, sizeof(int));
                    modval=particles[i]->getVeloc()[2]*dumpmult;
                    trajstream.write((char*)&modval, sizeof(int));
                }
            }
            int modval=box[0]*dumpmult;
            trajstream.write((char*)&modval, sizeof(int));
            modval=box[1]*dumpmult;
            trajstream.write((char*)&modval, sizeof(int));
            modval=box[2]*dumpmult;
            trajstream.write((char*)&modval, sizeof(int));

        }
        else{
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
        }
    }
    return;
}

void Dump::dumpForce(int step, real t){
    decomp->gatherForces();
    if(mpi->isMaster()){
        if(dumpbinary){
            forcestream.write((char*)&step, sizeof(int));
            forcestream.write((char*)&t, sizeof(float));
            if(!dumpfrozen)
                forcestream.write((char*)&nufbeads, sizeof(int));
            else{
                int ntotparticles=particles.size();
                forcestream.write((char*)&ntotparticles, sizeof(int));
            }
            for(int i=0;i<particles.size();i++){
                if(particles[i]->isFrozen() && !dumpfrozen){
                }
                else{
                    int modval=particles[i]->getForce()[0]*dumpmult;
                    forcestream.write((char*)&modval, sizeof(int));
                    modval=particles[i]->getForce()[1]*dumpmult;
                    forcestream.write((char*)&modval, sizeof(int));
                    modval=particles[i]->getForce()[2]*dumpmult;
                    forcestream.write((char*)&modval, sizeof(int));
                }
            }
        }
        else{
            forcestream << "Force at step "  << step << ", time " << t << std::endl;
            if(!dumpfrozen)
                forcestream << nufbeads << std::endl;
            else
                forcestream << particles.size() << std::endl;
            for(int i=0;i<particles.size();i++){
                if(particles[i]->isFrozen() && !dumpfrozen){
                }
                else{
                    forcestream << std::setw(5) << particles[i]->getMoleculeIndex();
                    forcestream << std::left << std::setw(5) << particles[i]->getMoleculeName();
                    forcestream << std::right << std::setw(5) << particles[i]->getTypeName();
                    forcestream << std::setw(5) << particles[i]->getParticleIndex();
                    forcestream << std::setprecision(3) << std::setw(8);
                    forcestream << std::setw(12) << std::fixed << std::setprecision(3) << particles[i]->getForce()[0];
                    forcestream << std::setw(12) << std::fixed << std::setprecision(3) << particles[i]->getForce()[1];
                    forcestream << std::setw(12) << std::fixed << std::setprecision(3) << particles[i]->getForce()[2];
                    forcestream << std::endl;
                }
            }
        }
    }
    return;
}

void Dump::dumpLocalDensity(){
    decomp->gatherDensities();
    if(mpi->isMaster()){
        std::cout << "*** Local Density ***" << std::endl;
        for(int i=0;i<particles.size();i++){
            if(particles[i]->isFrozen() && !dumpfrozen){
            }
            else{
                std::cout << std::setw(5) << particles[i]->getMoleculeIndex();
                std::cout << std::left << std::setw(5) << particles[i]->getMoleculeName();
                std::cout << std::right << std::setw(5) << particles[i]->getTypeName();
                std::cout << std::setw(5) << particles[i]->getParticleIndex();
                std::cout << std::setprecision(3) << std::setw(8);
                std::cout << std::setw(8) << std::fixed << std::setprecision(3) << particles[i]->getDensity();
                std::cout << std::endl;
            }
        }
    }
    return;
}

void Dump::writeLogHeader(){
    logstream << "#" << std::setw(7) << "Step";
    logstream << std::setw(13) << "Time";
    logstream << std::setw(13) << "Temperature";
    logstream << std::setw(13) << "Pressure";
    logstream << std::setw(13) << "Stress_xx";
    logstream << std::setw(13) << "Stress_yy";
    logstream << std::setw(13) << "Stress_zz";
    logstream << std::setw(13) << "Stress_xy";
    logstream << std::setw(13) << "Stress_xz";
    logstream << std::setw(13) << "Stress_yx";
    logstream << std::setw(13) << "Stress_yz";
    logstream << std::setw(13) << "Stress_zx";
    logstream << std::setw(13) << "Stress_zy";
    logstream << std::setw(15) << "SimTime_left";
    logstream << std::endl;
    return;
}


void Dump::dumpLog(int step, real time, std::string remaint){
    if(mpi->isMaster()){
        Rvec virial=analysis->getStressTensor();
        logstream << std::setw(8) << step;
        logstream << std::setw(13) << std::fixed << std::setprecision(4) << time;
        logstream << std::setw(13) << std::fixed << std::setprecision(4) << std::scientific << analysis->getTemperature();
        logstream << std::setw(13) << std::fixed << std::setprecision(4) << std::scientific << analysis->getPressure();
        logstream << std::setw(13) << std::fixed << std::setprecision(4) << std::scientific << virial[0];
        logstream << std::setw(13) << std::fixed << std::setprecision(4) << std::scientific << virial[1];
        logstream << std::setw(13) << std::fixed << std::setprecision(4) << std::scientific << virial[2];
        logstream << std::setw(13) << std::fixed << std::setprecision(4) << std::scientific << virial[3];
        logstream << std::setw(13) << std::fixed << std::setprecision(4) << std::scientific << virial[4];
        logstream << std::setw(13) << std::fixed << std::setprecision(4) << std::scientific << virial[5];
        logstream << std::setw(13) << std::fixed << std::setprecision(4) << std::scientific << virial[6];
        logstream << std::setw(13) << std::fixed << std::setprecision(4) << std::scientific << virial[7];
        logstream << std::setw(13) << std::fixed << std::setprecision(4) << std::scientific << virial[8];
        logstream << std::setw(15) << std::fixed << remaint;
        logstream << std::endl;
    }
    return;
}
    

void Dump::dumpStress(int step, real t){
    decomp->gatherStresses();
    if(mpi->isMaster()){
        if(dumpbinary){
            stressstream.write((char*)&step, sizeof(int));
            stressstream.write((char*)&t, sizeof(float));
            if(!dumpfrozen)
                stressstream.write((char*)&nufbeads, sizeof(int));
            else{
                int ntotparticles=particles.size();
                stressstream.write((char*)&ntotparticles, sizeof(int));
            }
            for(int i=0;i<particles.size();i++){
                if(particles[i]->isFrozen() && !dumpfrozen){
                }
                else{
                    int modval=particles[i]->getStress()[0]*dumpmult;
                    stressstream.write((char*)&modval, sizeof(int));
                    modval=particles[i]->getStress()[1]*dumpmult;
                    stressstream.write((char*)&modval, sizeof(int));
                    modval=particles[i]->getStress()[2]*dumpmult;
                    stressstream.write((char*)&modval, sizeof(int));
                    modval=particles[i]->getStress()[3]*dumpmult;
                    stressstream.write((char*)&modval, sizeof(int));
                    modval=particles[i]->getStress()[4]*dumpmult;
                    stressstream.write((char*)&modval, sizeof(int));
                    modval=particles[i]->getStress()[5]*dumpmult;
                    stressstream.write((char*)&modval, sizeof(int));
                    modval=particles[i]->getStress()[6]*dumpmult;
                    stressstream.write((char*)&modval, sizeof(int));
                    modval=particles[i]->getStress()[7]*dumpmult;
                    stressstream.write((char*)&modval, sizeof(int));
                    modval=particles[i]->getStress()[8]*dumpmult;
                    stressstream.write((char*)&modval, sizeof(int));
                }
            }
        }

        else{
            stressstream << "Step=" << step << ", Time=" << t << std::endl;
            if(!dumpfrozen)
                stressstream << nufbeads << std::endl;
            else
                stressstream << particles.size() << std::endl;
            for(int i=0;i<particles.size();i++){
                if(particles[i]->isFrozen() && !dumpfrozen){
                }
                else{
                    stressstream << std::setw(5) << particles[i]->getParticleIndex();
                    stressstream << std::setw(15) << std::fixed << std::setprecision(5) << std::scientific  << particles[i]->stress[0];
                    stressstream << std::setw(15) << std::fixed << std::setprecision(5) << std::scientific  << particles[i]->stress[1];
                    stressstream << std::setw(15) << std::fixed << std::setprecision(5) << std::scientific  << particles[i]->stress[2];
                    stressstream << std::setw(15) << std::fixed << std::setprecision(5) << std::scientific  << particles[i]->stress[3];
                    stressstream << std::setw(15) << std::fixed << std::setprecision(5) << std::scientific  << particles[i]->stress[4];
                    stressstream << std::setw(15) << std::fixed << std::setprecision(5) << std::scientific  << particles[i]->stress[5];
                    stressstream << std::setw(15) << std::fixed << std::setprecision(5) << std::scientific  << particles[i]->stress[6];
                    stressstream << std::setw(15) << std::fixed << std::setprecision(5) << std::scientific  << particles[i]->stress[7];
                    stressstream << std::setw(15) << std::fixed << std::setprecision(5) << std::scientific  << particles[i]->stress[8];
                    stressstream << std::endl;
                }
            }
        }
    }
    return;
}


