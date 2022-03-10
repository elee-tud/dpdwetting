#include "property.hpp"
#include "../filecontrol.hpp"

using namespace dpd;

Property::Property(InitialSet initset):initset(initset){
    command=initset.command;
    control=initset.control;
    topol=initset.topol;
    config=initset.config;
    inter=initset.interactions;
    mpi=initset.mpi;
    
    for_distrib=false;
    for_timeevol=false;
    for_dynamics=false;
    need_multicore=false;

    need_position=false;
    need_velocity=false;
    need_stress=false;
    need_force=false;

    dbegin=0.0;

    box=config->getBox();
    pbc=PeriodicBoundary(box);
    particles=config->getParticles();

}

void Property::printPropertyTitle(){
    if(mpi->isMaster()){
        
        std::cout << "***************************************************************************"<< std::endl;
        std::cout << title << std::endl;
        std::cout << "***************************************************************************"<< std::endl;
    }
    return;
}


void Property::getGeneralParameters(){
    if(mpi->isMaster()){
        command->getCommandSingleOption("-bs", 0, &begstep);
        command->getCommandSingleOption("-es", -1, &endstep);
        command->getCommandSingleOption("-ss", 1, &skipstep);
        
        if(endstep==-1){
            endstep=control->getTotalSteps()/control->getTrajFrequency();
        }
        nsteps=static_cast<int>((endstep-begstep)/skipstep)+1;

    }

    MPI_Bcast(&begstep, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&endstep, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&skipstep, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&nsteps, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

    if(control->getDumpFrozen())
        nbeads=topol->getNbeads();
    else
        nbeads=topol->getNumUnfrozenBeads();
    

    return;
}


void Property::initializeResultArrays(){
    if(for_timeevol){
        tevol=Rvec2D(nprops, Rvec(nsteps, 0.0));
        tevol_sum=Rvec2D(nprops, Rvec(nsteps, 0.0));
        for(int i=0;i<nprops;i++){
            for(int j=0;j<nsteps;j++){
                tevol[i][j]=0.0;
                tevol_sum[i][j]=0.0;
            }
        }
    }
    if(for_distrib){
        dist=Rvec2D(nprops, Rvec(ndbin, 0.0));
        dist_sum=Rvec2D(nprops, Rvec(ndbin, 0.0));
        for(int i=0;i<nprops;i++){
            for(int j=0;j<ndbin;j++){
                dist[i][j]=0.0;
                dist_sum[i][j]=0.0;
            }
        }
    }
    if(for_dynamics){
        generateDeltaSteps();
        tevol=Rvec2D(nprops, Rvec(ndsteps, 0.0));
        tevol_sum=Rvec2D(nprops, Rvec(ndsteps, 0.0));
    }

    return;
}

void Property::generateDeltaSteps(){
    real base=1.01;
    real nds=1.0;
    dsteps.push_back(0);
    dsteps.push_back(int(nds));
    while(1){
        nds=nds*base;
        if(int(nds)>nsteps-1)
            break;
        if(int(nds)!=dsteps.back()){
            dsteps.push_back(int(nds));
        }
    }

    ndsteps=dsteps.size();
    numavgs=Ivec(ndsteps, 0);
    for(int i=0;i<ndsteps;i++){
        if(nsteps>=dsteps[i]+MAXAVG)
            numavgs[i]=MAXAVG;
        else
            numavgs[i]=nsteps-dsteps[i];
    }
    return;
}

       
int Property::getTrajectoryType(){
    if(need_position|| need_velocity)
        return CRDVEL;
    else if(need_stress)
        return STRESS;
    else if(need_force)
        return FORCE;
    else
        return -1;
}


void Property::reduceResults(){
    if(need_multicore){
        if(for_timeevol)
            reduceTimeEvolution();
        if(for_distrib)
            reduceDistribution();
    }
    else{
        if(for_timeevol){
            for(int i=0;i<nprops;i++){
                for(int j=0;j<nsteps;j++){
                    tevol_sum[i][j]=tevol[i][j];
                }
            }
        }
        if(for_distrib){
            for(int i=0;i<nprops;i++){
                for(int j=0;j<ndbin;j++){
                    dist_sum[i][j]=dist[i][j];
                }
            }
        }
        if(for_dynamics){
            calculateDynamicProperty();
        }
    }


    return;
}



void Property::reduceTimeEvolution(){
    /*
    MPI_Barrier(MPI_COMM_WORLD);

    if(mpi->isMaster()){
        for(int i=0;i<nprops;i++){
            for(int j=0;j<nsteps_proc;j++){
                results[i][j*mpi->size()]=valvst[i][j];
            }
        }
        int nval;
        for(int i=1;i<mpi->size();i++){
            MPI_Recv(&nval, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            Rvec recv(nval);
            for(int j=0;j<nprops;j++){
                MPI_Recv(&recv[0], nval, MPI_DOUBLE, i, 1+j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for(int k=0;k<nval;k++)
                    results[j][k*mpi->size()+i]=recv[k];
            }
        }
    }
    else{
        MPI_Send(&nsteps_proc, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
        for(int j=0;j<nprops;j++){
            MPI_Send(&valvst[j][0], nsteps_proc, MPI_DOUBLE, MASTER, 1+j, MPI_COMM_WORLD);
        }
    }
    */
    return;
}

void Property::reduceDistribution(){
    for(int i=0;i<nprops;i++){
        MPI_Reduce(&dist[i][0], &dist_sum[i][0], ndbin, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
    }
    return;
                
}
void Property::normalizeResults(){
    if(mpi->isMaster() && for_distrib){
        for(int i=0;i<nprops;i++){
            for(int j=0;j<ndbin;j++){
                dist_sum[i][j]/=nsteps*dbin*numvarperstep;
            }
        }
    }
    if(mpi->isMaster() && for_dynamics){
        for(int i=0;i<nprops;i++){
            for(int j=0;j<ndsteps;j++){
//                tevol_sum[i][j]/=numavgs[j];
            }
                
        }
    }
    return;
}


void Property::resetBox(Real3D newbox){
    box=newbox;
    pbc.resetBoxSize(box);
    return;
}



void Property::generateOutFileName(){
    if(mpi->isMaster()){
        outname=command->output_analysis();
        if(outname.compare("default")==0){
            if(for_distrib && for_timeevol){
                size_t idx=outfile.find(".out");
                if(idx!=std::string::npos){
                    outtimename=outfile;
                    outdistname=outfile.replace(idx, 5, ".dist");
                }
                else{
                    idx=outfile.find(".dist");
                    if(idx!=std::string::npos){
                        outdistname=outfile;
                        outtimename=outfile.replace(idx, 5, ".out");
                    }
                    else{
                        outtimename=outfile+".out";
                        outdistname=outfile+".dist";
                    }
                }
            }
            else{
                outtimename=outfile;
                outdistname=outfile;
                outdynmname=outfile;
            }

        }
        else{
            outtimename=outname;
            outdistname=outname;
            outdynmname=outname;
        }

        generateOutputTitleHeader();
    }
    return;
}

void Property::generateOutputTitleHeader(){
    if(getTrajectoryType()==CRDVEL)
        dt=control->getTimeStep()*control->getTrajFrequency();
    else if(getTrajectoryType()==STRESS)
        dt=control->getTimeStep()*control->getStressFrequency();
    else if(getTrajectoryType()==FORCE)
        dt=control->getTimeStep()*control->getForceFrequency();

    if(for_timeevol){
        outtimeheader=generateOutputHeader(outheader_tevol);
    }
    if(for_distrib){
        outdistheader=generateOutputHeader(outheader_dist);
    }
    if(for_dynamics){
        outdynmheader=generateOutputHeader(outheader_dynm);
    }

}

std::string Property::generateOutputHeader(Svec hd){
    std::stringstream header;
    header << "#";
    header << std::setw(15) << hd[0];
    for(int i=1;i<hd.size();i++){
        header << std::setw(16) << hd[i];
    }
    return header.str();
}

void Property::writeOutput(){
    if(mpi->isMaster()){
        if(for_timeevol)
            writeTimeEvolution();
        if(for_distrib)
            writeDistribution();
        if(for_dynamics)
            writeDynamics();
    }
    return;

}

void Property::writeTimeEvolution(){
    std::ofstream outstream;
    openFileWithBackup(outtimename, &outstream, false, false);
    outstream << outtitle_tevol << std::endl;
    outstream << outtimeheader << std::endl;
    int nrow=tevol_sum.size();
    int ncol=tevol_sum[0].size();
    for(int i=0;i<ncol;i++){
        outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << (begstep+(i*skipstep))*dt ;
        for(int j=0;j<nrow;j++){
            outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << tevol_sum[j][i];
        }
        outstream << std::endl;
    }
    outstream.close();
    std::cout << "Results are written in " << outtimename << "." << std::endl;
    return;
}

void Property::writeDistribution(){
    std::ofstream outstream;
    openFileWithBackup(outdistname, &outstream, false, false);
    outstream << outtitle_dist << std::endl;
    outstream << outdistheader << std::endl;
    int nrow=dist_sum.size();
    int ncol=dist_sum[0].size();
    for(int i=0;i<ncol;i++){
        outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << dbegin+(i+0.5)*dbin;
        for(int j=0;j<nrow;j++){
            outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << dist_sum[j][i];
        }
        outstream << std::endl;
    }
    outstream.close();
    std::cout << "Results are written in " << outdistname << "." << std::endl;
    return;
}


void Property::writeDynamics(){
    std::ofstream outstream;
    openFileWithBackup(outdynmname, &outstream, false, false);
    outstream << outtitle_dynm << std::endl;
    outstream << outdynmheader << std::endl;
    int nrow=tevol_sum.size();
    int ncol=tevol_sum[0].size();
    for(int i=0;i<ncol;i++){
        outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << dsteps[i]*dt;
        for(int j=0;j<nrow;j++){
            outstream << std::setw(16) << std::fixed << std::setprecision(5) << std::scientific << tevol_sum[j][i];
        }
        outstream << std::endl;
    }
    outstream.close();
    std::cout << "Results are written in " << outdynmname << "." << std::endl;
    return;
}

