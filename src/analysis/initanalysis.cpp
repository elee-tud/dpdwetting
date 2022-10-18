#include "initanalysis.hpp"

using namespace dpd;
InitAnalysis::InitAnalysis(Command* command):Initialization(command){

    initset=InitialSet{.command=command, .control=control, .topol=topol, .config=config, .mpi=mpi, .decomp=decomp, .interactions=interactions};
    program=command->getProgram();
    loadProperty();
    checkNumberOfCoresForProgram();
    prop->printPropertyTitle();
    prop->getGeneralParameters();
    prop->getSpecificParameters();
    prop->initializeVariables();
    traj=new Trajectory(initset, prop->getTrajectoryType()); 
    nbeads=prop->getNbeads();
    nsteps=prop->getNsteps();
    begstep=prop->getBeginStep();
    endstep=prop->getEndStep();
    skipstep=prop->getSkipStep();
    timer->reset(nsteps);


}


void InitAnalysis::loadProperty(){
    if(program==DROPSIZE){
        prop=new DropSize(initset);
    }
    else if(program==VELOCITY){
        prop=new DropVelocity(initset);
    }
    else if(program==SPHERICALSTRESS){
        prop=new DropSphericalStress(initset);
    }
    else if(program==RADIALDENSITY){
        prop=new RadialDensity(initset);
    }
    else if(program==AVGSTRESS){
        prop=new AverageStress(initset);
    }
    else if(program==POLADSORP){
        prop=new PolymerAdsorption(initset);
    }
    else if(program==POLSIZE){
        prop=new PolymerSize(initset);
    }
    else if(program==POLEVRLX){
        prop=new PolymerEevecRelax(initset);
    }
    else if(program==POLORIENT){
        prop=new PolymerOrientation(initset);
    }
    else if(program==BONDLENGTH){
        prop=new PolymerBondlength(initset);
    }
    else if(program==POLSTRETCH){
        prop=new PolymerStretch(initset);
    }
    else if(program==TRJTOGRO){
        prop=new TrjToGro(initset);
    }
    else if(program==POLSMSF){
        prop=new PolymerStructFactor(initset);
    }
    else if(program==POLSUBSIZE){
        prop=new PolymerInternalSize(initset);
    }
    else if(program==POLMSD){
        prop=new MeanSquareDisplacement(initset);
    }
    else if(program==VELACF){
        prop=new VelocityAutoCorrelation(initset);
    }
    else if(program==SURFCOV){
        prop=new SurfaceCoverage(initset);
    }
    else if(program==RDF){
        prop=new RadialDistributionFunction(initset);
    }
    else if(program==DEPORIENT){
        prop=new PolymerDepositOrient(initset);
    }
    else if(program==BRIDGESIZE){
        prop=new BridgeSize(initset);
    }
    else if(program==BRIDGEVEL){
        prop=new BridgeVelocity(initset);
    }
    else if(program==BRIDGEPCONC){
        prop=new BridgePolymerConc(initset);
    }
    else if(program==BRIDGEADSCONC){
        prop=new BridgeAdsorptionConc(initset);
    }
    else if(program==BRIDGECLVEL){
        prop=new BridgeContLineVelocity(initset);
    }
    else if(program==DROPZD){
        prop=new DropZDensity(initset);
    }
    else if(program==DROPSRDF){
        prop=new SurfaceRdf(initset);
    }
    else if(program==JUMPFREQ){
        prop=new JumpingFrequency(initset);
    }
    else if(program==BRIDGEJF){
        prop=new BridgeJumpingFrequency(initset);
    }
    else if(program==BRIDGESLVEL){
        prop=new BridgeSlipVelocity(initset);
    }
    else if(program==BRIDGEVELX){
        prop=new BridgeVelocityX(initset);
    }
    else if(program==BRIDGEVELXZ){
        prop=new BridgeVelocityXZ(initset);
    }
    else if(program==BRIDGEALDENS){
        prop=new BridgeAdsorblayerDensity(initset);
    }
    else if(program==BRIDGEZD){
        prop=new BridgeZdensity(initset);
    }
    else if(program==BRIDGEDENSXZ){
        prop=new BridgeDensityXZ(initset);
    }
    else if(program==BRIDGECLINE){
        prop=new BridgeCline(initset);
    }
    else if(program==BRIDGEINTERF){
        prop=new BridgeInterface(initset);
    }
    else if(program==BRIDGEGPDEN){
        prop=new BridgeGrooveParticleDensity(initset);
    }
    else if(program==NUMCLUSTER){
        prop=new NumberOfLiquidClusters(initset);
    }
    return;
}

void InitAnalysis::checkNumberOfCoresForProgram(){
    if(mpi->size()>1 && !prop->needMultiCore()){
        if(mpi->isMaster())
            std::cout << " [Error] This program has no advantage of using multicore processor.\n         Use single core for this property" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 0);
        exit(0);
    }

    return;
}

InitAnalysis::~InitAnalysis(){
    delete traj;
//    delete prop;
}
