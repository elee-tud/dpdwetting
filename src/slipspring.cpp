#include "slipspring.hpp" 
#include <algorithm>
#include <chrono>
#include <thread>
#include "filecontrol.hpp"
#include <sstream>
#include <iostream>
#include <iomanip>
using namespace dpd;

SlipSpring::SlipSpring(Control* control, Topology* topol, Configuration* config, Decomposition* decomp):control(control), topol(topol), config(config), decomp(decomp){
    mpi=decomp->getMPI();
    particles=config->getParticles();
    ntotptcls=particles.size();
    for(int i=0;i<ntotptcls;i++){
        particles[i]->reserveSSBond(MAXSS);
    }
    totnumss=control->getNumSlipSprings();
    nmcsteps=control->getNumMCSeqSteps();
    ndpdsteps=control->getNumDPDSeqSteps();
    ssl0=control->getSlipSpringLength();
    ssk=control->getSlipSpringForceConst();
    ssmaxd=control->getSlipSpringMaxCutoff();
    ssmind=control->getSlipSpringMinCutoff();
    temp=control->getTemperature();
    sspoltype=control->getSlipSpringPolymerType();
    intrassbias=control->getIntraSlipSpringBias();

    box=config->getBox();
    pbc=PeriodicBoundary(box);
    ssdecomp=new Decomposition(control, config, mpi, 1, ssmaxd);
    cells=ssdecomp->getCells();

    nssproc=Ivec(mpi->size(),0);
    startidx=Ivec(mpi->size(),0);
    for(int i=0;i<mpi->size();i++){
        nssproc[i]=static_cast<int>(totnumss/mpi->size());
        if(i<totnumss%mpi->size()){
            nssproc[i]++;
        }
        for(int j=0;j<i;j++){
            startidx[i]+=nssproc[j];
        }
    }
    numss=nssproc[mpi->rank()];
    springs=Ivec2D(numss, Ivec(2,0));

    endindex.reserve(particles.size());
    for(int i=0;i<ntotptcls;i++){
        if(particles[i]->getBonds().size()==1)
            endindex.push_back(i);
    }

    numends=endindex.size();
    getPolymerLength(); 


    srand(control->getRandomSeed()*mpi->rank());
    restart=control->getRestart();
    if(mpi->isMaster()){
        ssfname=control->getCommand()->slipspring();
        dpd::openFileWithBackup(ssfname, &ssstream, false, restart);
        sstrjfname=control->getCommand()->slipspringtrj();
        dpd::openFileWithBackup(sstrjfname, &sstrjstream, false, restart);
        if(!restart)
            writeLogHeader();
    }

}

SlipSpring::~SlipSpring(){
    ssstream.close();
    sstrjstream.close();
    return;
}

void SlipSpring::getPolymerLength(){
    if(mpi->isMaster()){
        pollength=0;

        bool first=true;
        int previdx=-1;
        for(int i=0;i<particles.size();i++){
            if(particles[i]->getMoleculeName().compare("POL")==0){
                if(first){
                    previdx=particles[i]->getMoleculeIndex();
                    first=false;
                }
                if(particles[i]->getMoleculeIndex()==previdx)
                    pollength++;
                else
                    break;
            }
        }
    }
    MPI_Bcast(&pollength, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    return;
}

void SlipSpring::initializeSlipSpring(){
    Ivec initbeads;
    saveCellIndexOfParticles();
    allsprings=Ivec2D(totnumss, Ivec(2, 0));
    if(mpi->isMaster()){
        for(int i=0;i<totnumss;i++){
            while(true){
                int pick=rand()%ntotptcls;
                int pair=choosePairRandom(pick);
                if(pair>-1){
                    allsprings[i]=Ivec{pick, pair};
                    break;
                }
            }

        }
        

    }

    bcastSprings();
    /*
    if(mpi->isMaster()){
        for(int i=0;i<totnumss;i++){
            std::cout << "(" << allsprings[i][0] << "," << allsprings[i][1] << ")" ;
            std::cout << ":[" << allsprings[i][0]/25+1 << "," << allsprings[i][1]/25+1 << "]" ;
        }
        std::cout << std::endl;
    }
    */
    
    calculateTotalEnergy();
    doMonteCarloSequence(0);

    accepted=0;


    return;
}

void SlipSpring::saveCellIndexOfParticles(){
    ssdecomp->clearAllBeadInformation();
    ssdecomp->allocateBeadsToDomain();
    ssdecomp->allocateBeadsToCells();
    ssdecomp->communicateGhostBeads();
//    if(mpi->isMaster()) ssdecomp->printCellInformation();

    cells=ssdecomp->getCells();
    for(int i=0;i<cells.size();i++){
        if(cells[i]->isTrueCell()){
//            if(mpi->isMaster()){
//                std::cout << "HERE";
//                cells[i]->printBeads();
//            }
            Ivec beads=cells[i]->getBeads();
            for(int j=0;j<beads.size();j++){
                particles[beads[j]]->inTheCell(i);
            }
        }
    }
    /*
    if(mpi->isMaster()){
        for(int i=0;i<particles.size();i++){
            std::cout << i << ":" << particles[i]->inWhichCell() << std::endl;
        }
    }
    */
    return;
}




int SlipSpring::choosePairRandom(int index1){
    int cellidx=particles[index1]->inWhichCell();
    Ivec beads=cells[cellidx]->getBeads();
    Ivec pospair;
    pospair.reserve(100);
    for(int i=0;i<beads.size();i++){
        if(beads[i]!=index1){
            real dist=pbc.getMinimumDistance(particles[index1]->coord, particles[beads[i]]->coord);
            if(dist>=ssmind && dist < ssmaxd){
                pospair.push_back(beads[i]);
            }
        }
    }
    Ivec nbcells=cells[cellidx]->getNeighborCells();
    for(int j=0;j<nbcells.size();j++){
        beads=cells[nbcells[j]]->getBeads();
        for(int i=0;i<beads.size();i++){
            real dist=pbc.getMinimumDistance(particles[index1]->coord, particles[beads[i]]->coord);
            if(dist>=ssmind && dist < ssmind){
                pospair.push_back(beads[i]);
            }
        }
    }
    int npospair=pospair.size();
    int pick=rand()%npospair;
//    std::cout << pospair << std::endl;
    if(pospair.size()>0)
        return pospair[pick];
    else
        return -1;
}


void SlipSpring::bcastSprings(){
    for(int i=0;i<totnumss;i++){
        MPI_Bcast(&allsprings[i][0], 2, MPI_INT, MASTER, MPI_COMM_WORLD);
    }
    for(int j=0;j<numss;j++){
        springs[j]=allsprings[startidx[mpi->rank()]+j];
    }
    return;


}

void SlipSpring::calculateTotalEnergy(){
    ssenergy=0;
    for(int i=0;i<numss;i++){
        real dist=pbc.getMinimumDistance(particles[springs[i][0]]->coord, particles[springs[i][1]]->coord);
        ssenergy+=springEnergy(dist);
    }
    sum_ssenergy=0;
    MPI_Allreduce(&ssenergy, &sum_ssenergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return;
}

void SlipSpring::doMonteCarloSequence(int dpdstep){
    /*
    decomp->gatherCoords();
    config->bcastConfiguration();
    saveCellIndexOfParticles();
    calculateTotalEnergy();
    */
    for(int i=0;i<nmcsteps;i++){
        doMonteCarloStep();
    }
    reduceSprings();
//    if(mpi->isMaster())
//        std::cout << allsprings << std::endl;
    addSpringsToConfiguration();
    sum_ssenergy=0;
    MPI_Allreduce(&ssenergy, &sum_ssenergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(mpi->isMaster()){
        dumpLog(dpdstep);
        dumpTrajectory(dpdstep);
    }
    return;
}


void SlipSpring::doMonteCarloStep(){
    accepted=0;
    for(int i=0;i<numss;i++){
        /*
        int spidx=rand()%numss;
        int spdir=rand()%2;
        int pick=springs[spidx][spdir];
        ParticleList bonded=particles[pick]->getBonds();
        int nbonded=bonded.size();
        int move=rand()%nbonded;
        */
        int bead1=springs[i][0];
        int bead2=springs[i][1];
        ParticleList b1bonded=particles[bead1]->getBonds();
        ParticleList b2bonded=particles[bead2]->getBonds();
        int move1=rand()%b1bonded.size();
        int move2=rand()%b2bonded.size();
        Particle* newbead1=b1bonded[move1];
        Particle* newbead2=b2bonded[move2];

//        std::cout << spidx << "," << spdir << "," << pick << "," << move << "," <<  std::endl;
        real newdist=pbc.getMinimumDistance(newbead1->coord, newbead2->coord);
//        std::cout << newdist << std::endl;
        if(newbead1!=newbead2 && newdist<ssmaxd && newdist >=ssmind){
            real olddist=pbc.getMinimumDistance(particles[bead1]->coord, particles[bead2]->coord);
//            std::cout << olddist << std::endl;
            real oldener=springEnergy(olddist);
            real newener=springEnergy(newdist); 
            real de=newener-oldener;
//            std::cout << de << std::endl;
            if(static_cast<real>(rand()/static_cast<real>(RAND_MAX)) < exp(-de/temp)){
                springs[i][0]=newbead1->getParticleIndex()-1;
                springs[i][1]=newbead2->getParticleIndex()-1;
                ssenergy+=de;
                accepted++;
//                  std::cout << "accepted" << std::endl;
            }
//            else
//                std::cout << "rejected" << std::endl;
        }
    }
    sum_accepted=0;
    MPI_Reduce(&accepted, &sum_accepted, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
    
    return;
}

void SlipSpring::reduceSprings(){
    if(mpi->isMaster()){
        for(int j=0;j<numss;j++){
           allsprings[j]=springs[j];
        }
        for(int i=1;i<mpi->size();i++){
            for(int j=0;j<nssproc[i];j++){
                MPI_Recv(&allsprings[startidx[i]+j][0], 2, MPI_INT, i, j,  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
    else{
        for(int j=0;j<numss;j++){
            MPI_Send(&springs[j][0], 2, MPI_INT, MASTER, j, MPI_COMM_WORLD);
        }
    }
    for(int i=0;i<totnumss;i++){
        MPI_Bcast(&allsprings[i][0], 2, MPI_INT, MASTER, MPI_COMM_WORLD);
    }
    return;
}




void SlipSpring::addSpringsToConfiguration(){
    removeAllSpringsFromConfiguration();
    for(int i=0;i<totnumss;i++){
        particles[allsprings[i][0]]->addSSBond(particles[allsprings[i][1]]);
        particles[allsprings[i][1]]->addSSBond(particles[allsprings[i][0]]);
    }
//    std::cout << mpi->rank() <<" :" <<  totnumss  << std::endl;
    return;
}

void SlipSpring::removeAllSpringsFromConfiguration(){
    for(int i=0;i<ntotptcls;i++){
        particles[i]->removeAllSSBond();
    }
    return;
}

void SlipSpring::relocation(){
    decomp->gatherCoords();
    config->bcastConfiguration();
    saveCellIndexOfParticles();
    calculateTotalEnergy();
    relocated=0;
    numintrass=0.;

    if(sspoltype==LINEARSS)
        relocationForLinearPolymer();
    else if(sspoltype==RINGSS)
        relocationForRingPolymer();

    return;
}

void SlipSpring::relocationForLinearPolymer(){
//    if(mpi->isMaster())
//        std::cout << "blc" << springs << std::endl;
    for(int i=0;i<numss;i++){
        if(particles[springs[i][0]]->getBonds().size()==1 || particles[springs[i][1]]->getBonds().size()==1){
            int pick=endindex[rand()%numends];
            int pair=choosePairRandom(pick);
            if(pair>-1){
                real newdist=pbc.getMinimumDistance(particles[pick]->coord, particles[pair]->coord);
                if(newdist<ssmaxd && newdist >=ssmind){
                    real olddist=pbc.getMinimumDistance(particles[springs[i][0]]->coord, particles[springs[i][1]]->coord);
                    real oldener=springEnergy(olddist);
                    real newener=springEnergy(newdist); 
                    real de=newener-oldener;
                    if(static_cast<real>(rand()/static_cast<real>(RAND_MAX)) < exp(-de/temp)){
                        springs[i][0]=pick;
                        springs[i][1]=pair;
                        ssenergy+=de;
                        relocated++;
                    }
                }
            }
        }
    }
    for(int i=0;i<numss;i++){
        if(particles[springs[i][0]]->getMoleculeIndex()==particles[springs[i][1]]->getMoleculeIndex())
            numintrass++;
    }
//        if(mpi->isMaster()) std::cout << "RELOCATIION" << springs[0] << std::endl;
    sum_ssenergy=0;
    MPI_Allreduce(&ssenergy, &sum_ssenergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sum_relocated=0;
    MPI_Reduce(&relocated, &sum_relocated, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
    sum_numintrass=0;
    MPI_Reduce(&numintrass, &sum_numintrass, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
//    if(mpi->isMaster())
//        std::cout << "after relocation" << springs << std::endl;
//    if(mpi->isMaster()) std::cout << sum_relocated << std::endl;
//    if(mpi->isMaster())
//        std::cout << "alc" << springs << std::endl;
    
    return;
}


void SlipSpring::findPseudoEnds(){
    reduceSprings();
    endindex.clear();
    if(mpi->isMaster()){
        for(int i=0;i<totnumss;i++){
            if(particles[allsprings[i][0]]->getMoleculeIndex()==particles[allsprings[i][1]]->getMoleculeIndex()){
                int molnum=particles[allsprings[i][0]]->getMoleculeIndex()-1;
                int idx1=allsprings[i][0]%pollength;
                int idx2=allsprings[i][1]%pollength;
//                std::cout << "SPRINGS=" << allsprings[i] << ":";
//                std::cout << "MONIDX=" << idx1 << "," << idx2 << "->" ;
                int psend=idx1+idx2;

                if(psend%2==0)
                    psend/=2;
                else
                    psend=psend/2+rand()%2;
                psend+=molnum*pollength;
//                std::cout << psend << "," ;
                if(std::find(endindex.begin(), endindex.end(), psend)==endindex.end()){
//                    std::cout << psend << std::endl;
                    endindex.push_back(psend);
                }

                psend=idx1+idx2+pollength;
                if(psend%2==0)
                    psend/=2;
                else
                    psend=psend/2+rand()%2;
                psend=psend%pollength+molnum*pollength;
//                std::cout << psend <<  std::endl;
                if(std::find(endindex.begin(), endindex.end(), psend)==endindex.end())
                    endindex.push_back(psend);
            }
            
        }
        numends=endindex.size();
    }
    MPI_Bcast(&numends, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    if(!mpi->isMaster())
        endindex=Ivec(numends, 0);
    MPI_Bcast(&endindex[0], numends, MPI_INT, MASTER, MPI_COMM_WORLD);
    return;
}





void SlipSpring::relocationForRingPolymer(){
    findPseudoEnds();
    for(int i=0;i<numss;i++){
        if(std::find(endindex.begin(), endindex.end(), springs[i][0])!=endindex.end() || std::find(endindex.begin(), endindex.end(), springs[i][1])!=endindex.end() ){
            int pick=endindex[rand()%numends];
            int pair=choosePairRandom(pick);
            int oldsstype=INTERSS;
            if(particles[springs[i][0]]->getMoleculeIndex()==particles[springs[i][1]]->getMoleculeIndex())
                oldsstype=INTRASS;

            if(pair>-1){
                real newdist=pbc.getMinimumDistance(particles[pick]->coord, particles[pair]->coord);
                if(newdist<ssmaxd && newdist >=ssmind){
                    real olddist=pbc.getMinimumDistance(particles[springs[i][0]]->coord, particles[springs[i][1]]->coord);
                    real oldener=springEnergy(olddist);
                    real newener=springEnergy(newdist); 
                    real de=newener-oldener;
                    real bias;
                    int newsstype=INTERSS;
                    if(particles[pick]->getMoleculeIndex()==particles[pair]->getMoleculeIndex())
                        newsstype=INTRASS;

                    if(oldsstype==newsstype)
                        bias=0.;
                    else{
                        if(oldsstype==INTERSS && newsstype==INTRASS)
                            bias=-intrassbias;
                        else
                            bias=intrassbias;
                    }
//                    std::cout << springs[i] << "->" << pick << "," << pair << ":" << bias ;
                    if(static_cast<real>(rand()/static_cast<real>(RAND_MAX)) < exp(-(de+bias)/temp)){
                        springs[i][0]=pick;
                        springs[i][1]=pair;
                        ssenergy+=de;
                        relocated++;
//                        std::cout << "A" << std::endl;
                    }
//                    else
//                        std::cout << "R" << std::endl;
                }
            }
        }
    }
//        if(mpi->isMaster()) std::cout << "RELOCATIION" << springs[0] << std::endl;
    for(int i=0;i<numss;i++){
        if(particles[springs[i][0]]->getMoleculeIndex()==particles[springs[i][1]]->getMoleculeIndex())
            numintrass++;
    }
    sum_ssenergy=0;
    MPI_Allreduce(&ssenergy, &sum_ssenergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sum_relocated=0;
    MPI_Reduce(&relocated, &sum_relocated, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
    sum_numintrass=0;
    MPI_Reduce(&numintrass, &sum_numintrass, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
//    if(mpi->isMaster())
//        std::cout << "after relocation" << springs << std::endl;
//    if(mpi->isMaster()) std::cout << sum_relocated << std::endl;
//    if(mpi->isMaster())
//        std::cout << "alc" << springs << std::endl;
    return;
}
                

void SlipSpring::writeLogHeader(){
    ssstream << "#" << std::setw(14) << "Step";
    ssstream << std::setw(15) << "Time";
    ssstream << std::setw(15) << "E_ss";
    ssstream << std::setw(15) << "Accept_Ratio";
    ssstream << std::setw(15) << "Reloc_Ratio";
    ssstream << std::setw(15) << "Frac_intrass";
    ssstream << std::endl;
}

void SlipSpring::dumpLog(int step){
    ssstream << std::setw(15) << step;
    ssstream << std::setw(15) << std::fixed << std::setprecision(4) << step*control->getTimeStep();
    ssstream << std::setw(15) << std::fixed << std::setprecision(6) << sum_ssenergy;
    ssstream << std::setw(15) << std::fixed << std::setprecision(6) << static_cast<real>(accepted)/totnumss/nmcsteps;
    ssstream << std::setw(15) << std::fixed << std::setprecision(6) << static_cast<real>(sum_relocated)/totnumss;
    ssstream << std::setw(15) << std::fixed << std::setprecision(6) << static_cast<real>(sum_numintrass)/totnumss;
    ssstream << std::endl;
    return;
}

void SlipSpring::dumpTrajectory(int step){
    box=config->getBox();
    std::stringstream title;
    title << "Slip spring positions at step=" << step ;
    sstrjstream << title.str() << std::endl;
    sstrjstream << totnumss << std::endl;


    for(int i=0;i<totnumss;i++){
        for(int j=0;j<2;j++){
            sstrjstream << std::setw(5) << i+1;
            sstrjstream << std::left << std::setw(5) << "SS";
            sstrjstream << std::right << std::setw(5) << "S";
            sstrjstream << std::setw(5) << i*2+j+1;
            sstrjstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[allsprings[i][j]]->getCoord()[0];
            sstrjstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[allsprings[i][j]]->getCoord()[1];
            sstrjstream << std::setw(8) << std::fixed << std::setprecision(3) << particles[allsprings[i][j]]->getCoord()[2];
            sstrjstream << std::endl;
        }
    }
    sstrjstream << std::setw(10) << std::fixed << std::setprecision(5) << box[0];
    sstrjstream << std::setw(10) << std::fixed << std::setprecision(5) << box[1];
    sstrjstream << std::setw(10) << std::fixed << std::setprecision(5) << box[2] << std::endl;

    return;
}










