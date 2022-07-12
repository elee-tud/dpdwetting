#include "control.hpp"
#include "parsing.hpp"
#include <algorithm>
#include <sstream>
#include <cstring>

using namespace dpd;
Control::Control(Command* command, SetMPI* mpi):command(command), mpi(mpi){
    control_fname=command->control();
}

Control::Control(std::string control_fname, SetMPI* mpi):control_fname(control_fname), mpi(mpi){
}

void Control::setDefaults(){
    restart=command->doRestart();
    randomseed=time(NULL);
    dt=NECESSARY;
    temperature=NECESSARY;
    gamma=4.5;
    lambda=0.5;
    integrator=NECESSARY;
    trajfreq=100;
    stressfreq=-1;
    forcefreq=-1;
    logfreq=10;
    dumpbinary=true;
    saferatio=10.0;
    pull_springk=0.;
    pull_coord=Real3D(-256.0);
    pull_direct=0;

    gravity_field=0.;
    gravity_direct=-1;

    nonbonded=-1;
    bondlength=-1;
    walltype=NOWALL;
    walldirect=Ivec(3,0);
    wallmin=Rvec(3,-10.0);
    wallmax=Rvec(3,-10.0);
    dumpfrozen=false;
    rmcomvdir=Ivec{1,1,1};
    targetforce=100;
    rmcomvfreq=-1;

    deformation=NODEFORMATION;
    deformfreq=1;
    deformfactor=1.;

    pindirect=Ivec{0,0,0};
    pinindex=Ivec{};
    pinrandom=false;
    npinpermol=0;

    uniflow=Ivec{0, 0, 0};
    flowrate=0.;
    sheartensor=Rvec(9,0.);
    liquidgrps=Svec{};
    sheardir=-1;
    shearrate=0.;

    wsheardir=-1;
    wshearrate=0.;


    slipspring=false;
    numss=0;
    nsteps_ssmc=0;
    nsteps_ssdpd=0;
    ssl0=0.;
    ssk=0.;
    sspoltype=LINEARSS;
    intrassbias=0.;


    tempannrate=0.;
    return;
}

void Control::readControl(){
    if(mpi->isMaster()){
        try{
            constream.open(control_fname);
            if(!constream.is_open())
                throw control_fname;
        }catch(std::string control_fname){
            std::cout << "File name " << control_fname << "does not exist." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 12);
            exit(0);
        }
        constream.seekg(0, std::ios::beg);
        std::string line;
        while(std::getline(constream, line)){
            if(line[0]!=';' && line.length()!=0){
                tokens=dpd::parsing(line);

                if(tokens[0].compare("temperature")==0)
                    temperature=std::stod(tokens[1]);

                else if(tokens[0].compare("tempanneal")==0)
                    tempannrate=std::stod(tokens[1]);

                else if(tokens[0].compare("timestep")==0)
                    dt=std::stod(tokens[1]);

                else if(tokens[0].compare("gamma")==0)
                    gamma=std::stod(tokens[1]);

                else if(tokens[0].compare("saferatio")==0)
                    saferatio=std::stod(tokens[1]);

                else if(tokens[0].compare("maxforce")==0)
                    targetforce=std::stod(tokens[1]);

                else if(tokens[0].compare("lambda")==0)
                    lambda=std::stod(tokens[1]);

                else if(tokens[0].compare("cellcutoff")==0)
                    cellcut=std::stod(tokens[1]);

                else if(tokens[0].compare("totalsteps")==0)
                    totalsteps=std::stoi(tokens[1]);

                else if(tokens[0].compare("dumpbinary")==0){
                    if(tokens[1].compare("N")==0||
                            tokens[1].compare("n")==0||
                            tokens[1].compare("No")==0||
                            tokens[1].compare("no")==0)
                        dumpbinary=false;
                }

                else if(tokens[0].compare("xtrjfreq")==0)
                    trajfreq=std::stoi(tokens[1]);

                else if(tokens[0].compare("xlogfreq")==0)
                    logfreq=std::stoi(tokens[1]);

                else if(tokens[0].compare("xstrfreq")==0)
                    stressfreq=std::stoi(tokens[1]);

                else if(tokens[0].compare("xfrcfreq")==0)
                    forcefreq=std::stoi(tokens[1]);

                else if(tokens[0].compare("rmcomvfreq")==0)
                    rmcomvfreq=std::stoi(tokens[1]);

                else if(tokens[0].compare("rmcomvdir")==0){
                    rmcomvdir=Ivec{0,0,0};
                    for(int i=1;i<tokens.size();i++){
                        if(tokens[i].compare("x")==0){
                            rmcomvdir[0]=1;
                        }
                        if(tokens[i].compare("y")==0){
                            rmcomvdir[1]=1;
                        }
                        if(tokens[i].compare("z")==0){
                            rmcomvdir[2]=1;
                        }
                    }
                }
                else if(tokens[0].compare("randseed")==0)
                    randomseed=std::stoi(tokens[1]);
                else if(tokens[0].compare("dumpfrozen")==0){
                    if(tokens[1].compare("y")==0 ||
                            tokens[1].compare("Y")==0 ||
                            tokens[1].compare("Yes")==0 ||
                            tokens[1].compare("yes")==0 ||
                            tokens[1].compare("YES")==0 )
                        dumpfrozen=true;
                }
                else if(tokens[0].compare("integrator")==0){
                    if(tokens[1].compare("vv")==0)
                        integrator=VELVER;
                    else if(tokens[1].compare("emin")==0)
                        integrator=EMIN;
                    else if(tokens[1].compare("vv-sllod")==0)
                        integrator=SLLOD;
                }
                else if(tokens[0].compare("pullspringk")==0){
                    pull_springk=std::stod(tokens[1]);
                    if(tokens.size()==2)
                        pull_direct=0;
                    else{
                        if(tokens[2].compare("x")==0)
                            pull_direct=1;
                        else if(tokens[2].compare("y")==0)
                            pull_direct=2;
                        else if(tokens[2].compare("z")==0)
                            pull_direct=3;
                        else if(tokens[2].compare("xy")==0 || tokens[2].compare("yx")==0)
                            pull_direct=4;
                        else if(tokens[2].compare("xz")==0 || tokens[2].compare("zx")==0)
                            pull_direct=5;
                        else if(tokens[2].compare("yz")==0 || tokens[2].compare("zy")==0)
                            pull_direct=6;
                    }
                }

                else if(tokens[0].compare("pullcenter")==0){
                    pull_coord[0]=std::stod(tokens[1]);
                    pull_coord[1]=std::stod(tokens[2]);
                    pull_coord[2]=std::stod(tokens[3]);
                }

                else if(tokens[0].compare("gravity")==0){
                    gravity_field=std::stod(tokens[1]);
                    if(tokens[2].compare("x")==0)
                        gravity_direct=0;
                    if(tokens[2].compare("y")==0)
                        gravity_direct=1;
                    if(tokens[2].compare("z")==0)
                        gravity_direct=2;
                }


                else if(tokens[0].compare("nonbonded")==0){
                    if(tokens[1].compare("mdpd")==0)
                        nonbonded=NBMDPD;
                    else if(tokens[1].compare("dpd")==0)
                        nonbonded=NBDPD;
                    else{
                        std::cout << "The given value " << tokens[1] << " of the argument \"" << tokens[0] << "\" in the control file cannot be understood." << std::endl;
                    exit(0);
                    }
                }
                else if(tokens[0].compare("bondlength")==0){
                    if(tokens[1].compare("harmonic")==0)
                        bondlength=BLHARMONIC;
                    else{
                        std::cout << "The given value " << tokens[1] << " of the argument \"" << tokens[0] << "\" in the control file cannot be understood." << std::endl;
                    exit(0);
                    }
                }
                else if(tokens[0].compare("wall")==0){
                    if(tokens[1].compare("solid")==0)
                        walltype=SOLIDWALL;
                    for(int i=2;i<tokens.size();i++){
                        if(tokens[i].compare("x")==0){
                            walldirect[0]=1;
                        }
                        if(tokens[i].compare("y")==0){
                            walldirect[1]=1;
                        }
                        if(tokens[i].compare("z")==0){
                            walldirect[2]=1;
                        }
                    }
                }
                else if(tokens[0].compare("wallposx")==0){
                    wallmin[0]=std::stod(tokens[1]);
                    wallmax[0]=std::stod(tokens[2]);
                }

                else if(tokens[0].compare("wallposy")==0){
                    wallmin[1]=std::stod(tokens[1]);
                    wallmax[1]=std::stod(tokens[2]);
                }

                else if(tokens[0].compare("wallposz")==0){
                    wallmin[2]=std::stod(tokens[1]);
                    wallmax[2]=std::stod(tokens[2]);
                }

                else if(tokens[0].compare("deform")==0){
                    if(tokens[1].compare("elongation")==0)
                        deformation=ELONGATION;
                    else{
                        std::cout << "The given value " << tokens[1] << " of the argument \"" << tokens[0] << "\" in the control file cannot be understood." << std::endl;
                    exit(0);
                    }
                }
                else if(tokens[0].compare("deformfreq")==0){
                    deformfreq=std::stoi(tokens[1]);
                }
                else if(tokens[0].compare("deformfact")==0){
                    deformfactor=std::stod(tokens[1]);
                }

                else if(tokens[0].compare("pinning")==0){
                    pinmolname=tokens[1];
                }
                else if(tokens[0].compare("pindirect")==0){
                    for(int i=1;i<tokens.size();i++){
                        if(tokens[i].compare("x")==0){
                            pindirect[0]=1;
                        }
                        if(tokens[i].compare("y")==0){
                            pindirect[1]=1;
                        }
                        if(tokens[i].compare("z")==0){
                            pindirect[2]=1;
                        }
                    }
                }
                else if(tokens[0].compare("pinatom")==0){
                    if(tokens[1].compare("random")==0){
                        pinrandom=true;
                    }
                    else{
                        pinindex.clear();
                        for(int i=1;i<tokens.size();i++){
                            pinindex.push_back(std::stoi(tokens[i]));
                        }
                    }
                }
                else if(tokens[0].compare("npinpermol")==0){
                        npinpermol=std::stoi(tokens[1]);
                }
                /*
                else if(tokens[0].compare("pullgroup")==0){
                        pullgroup=tokens[1];
                }
                else if(tokens[0].compare("pulldirect")==0){
                    for(int i=1;i<tokens.size();i++){
                        if(tokens[i].compare("x")==0){
                            pulldirect[0]=1;
                        }
                        if(tokens[i].compare("y")==0){
                            pulldirect[1]=1;
                        }
                        if(tokens[i].compare("z")==0){
                            pulldirect[2]=1;
                        }
                    }
                }
                else if(tokens[0].compare("pullvelocity")==0){
                        pullvel=std::stof(tokens[1]);
                }
                else if(tokens[0].compare("pull")==0){
                        pullvel=std::stof(tokens[1]);
                }
                */


                else if(tokens[0].compare("uniflow")==0){
                    for(int i=1;i<tokens.size();i++){
                        if(tokens[i].compare("x")==0){
                            uniflow[0]=1;
                        }
                        else if(tokens[i].compare("y")==0){
                            uniflow[1]=1;
                        }
                        else if(tokens[i].compare("z")==0){
                            uniflow[2]=1;
                        }
                    }
                }
                
                else if(tokens[0].compare("flowrate")==0){
                    flowrate=std::stod(tokens[1]);
                    buildShearTensor();
                }
                else if(tokens[0].compare("shearflow")==0){
                    if(tokens[1].compare("xy")==0)
                        sheardir=0;
                    else if(tokens[1].compare("xz")==0)
                        sheardir=1;
                    else if(tokens[1].compare("yx")==0)
                        sheardir=2;
                    else if(tokens[1].compare("yz")==0)
                        sheardir=3;
                    else if(tokens[1].compare("zx")==0)
                        sheardir=4;
                    else if(tokens[1].compare("zy")==0)
                        sheardir=5;
                    shearrate=std::stod(tokens[2]);
                    buildShearTensor();
                }

                else if(tokens[0].compare("wallshear")==0){
                    if(tokens[1].compare("xy")==0)
                        wsheardir=0;
                    else if(tokens[1].compare("xz")==0)
                        wsheardir=1;
                    else if(tokens[1].compare("yx")==0)
                        wsheardir=2;
                    else if(tokens[1].compare("yz")==0)
                        wsheardir=3;
                    else if(tokens[1].compare("zx")==0)
                        wsheardir=4;
                    else if(tokens[1].compare("zy")==0)
                        wsheardir=5;
                    wshearrate=std::stod(tokens[2]);
                }


                else if(tokens[0].compare("wallshrgrp")==0){
                    shrwallgrp1=tokens[1];
                    shrwallgrp2=tokens[2];
                }


        

                else if(tokens[0].compare("liquidgrps")==0){
                    for(int i=1;i<tokens.size();i++){
                        liquidgrps.push_back(tokens[i]);
                    }
                }
                
                else if(tokens[0].compare("slipspring")==0){
                    numss=std::stoi(tokens[1]);
                    slipspring=true;

                }
                else if(tokens[0].compare("sscutoff")==0){
                    sscutmin=std::stod(tokens[1]);
                    sscutmax=std::stod(tokens[2]);
                }

                else if(tokens[0].compare("seqnmcsteps")==0){
                    nsteps_ssmc=std::stoi(tokens[1]);
                }
                else if(tokens[0].compare("seqndpdsteps")==0){
                    nsteps_ssdpd=std::stoi(tokens[1]);
                }
                else if(tokens[0].compare("ssparam")==0){
                    ssk=std::stod(tokens[1]);
                    ssl0=std::stod(tokens[2]);
                }
                else if(tokens[0].compare("sspoltype")==0){
                    if(tokens[1].compare("linear")==0)
                        sspoltype=LINEARSS;
                    else if(tokens[1].compare("ring")==0)
                        sspoltype=RINGSS;
                }
                else if(tokens[0].compare("intrassbias")==0){
                    intrassbias=std::stod(tokens[1]);
                }

                else{
                    std::cout << "Argument \"" << tokens[0] << "\" in the control file cannot be understood." << std::endl;
                    exit(0);
                }
            }
        }
        constream.close();
    }

            
    return;
}

void Control::setDependency(){
    if(mpi->isMaster()){
        if(isUniflow() || sheardir>-1){
            if(integrator!=SLLOD){
                std::cout << "The integrator SLLOD is used for the external flow." << std::endl;
                integrator=SLLOD;
            }
        }
    }
    return;
}


void Control::bcastControl(){
    MPI_Bcast(&restart, 1, MPI_CXX_BOOL, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&temperature, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&tempannrate, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&dt, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&gamma, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&lambda, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&saferatio, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&randomseed, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&integrator, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&totalsteps, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&dumpbinary, 1, MPI_CXX_BOOL, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&trajfreq, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&stressfreq, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&logfreq, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&forcefreq, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&cellcut, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&nonbonded, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&bondlength, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&dumpfrozen, 1, MPI_CXX_BOOL, MASTER, MPI_COMM_WORLD);

    MPI_Bcast(&pull_springk, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&pull_coord[0], 3, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&pull_direct, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

    MPI_Bcast(&gravity_field, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&gravity_direct, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

    MPI_Bcast(&targetforce, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&walltype, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    if(walltype==SOLIDWALL){
        MPI_Bcast(&walldirect[0], 3, MPI_INT, MASTER, MPI_COMM_WORLD);
        MPI_Bcast(&wallmin[0], 3, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
        MPI_Bcast(&wallmax[0], 3, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    }

    MPI_Bcast(&rmcomvfreq, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&rmcomvdir[0], 3, MPI_INT, MASTER, MPI_COMM_WORLD);
    
    MPI_Bcast(&deformation, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&deformfreq, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&deformfactor, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    MPI_Bcast(&pindirect[0], 3, MPI_INT, MASTER, MPI_COMM_WORLD);
    int size=pinmolname.size();
    char buf[8];
    std::strcpy(buf, pinmolname.c_str());
    MPI_Bcast(&size, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&buf[0], size, MPI_CHAR, MASTER, MPI_COMM_WORLD);
    pinmolname=std::string(buf).substr(0,size);
    size=pinindex.size();
    MPI_Bcast(&size, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&pinindex[0], size, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&pinrandom, 1, MPI_CXX_BOOL, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&npinpermol, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

   
    size=shrwallgrp1.size();
    std::strcpy(buf, shrwallgrp1.c_str());
    MPI_Bcast(&size, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&buf[0], size, MPI_CHAR, MASTER, MPI_COMM_WORLD);
    shrwallgrp1=std::string(buf).substr(0,size);

    size=shrwallgrp2.size();
    std::strcpy(buf, shrwallgrp2.c_str());
    MPI_Bcast(&size, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&buf[0], size, MPI_CHAR, MASTER, MPI_COMM_WORLD);
    shrwallgrp2=std::string(buf).substr(0,size);

    MPI_Bcast(&wshearrate, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&wsheardir, 1, MPI_INT, MASTER, MPI_COMM_WORLD);


    MPI_Bcast(&uniflow[0], 3, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&flowrate, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&sheardir, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&shearrate, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    if(isUniflow() || sheardir!=-1)
        MPI_Bcast(&sheartensor[0], 9, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);


    MPI_Bcast(&slipspring, 1, MPI_CXX_BOOL, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&numss, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&nsteps_ssmc, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&nsteps_ssdpd, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&ssl0, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&ssk, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&sscutmax, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&sscutmin, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&sspoltype, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&intrassbias, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);





    MPI_Barrier(MPI_COMM_WORLD);
    return;
}


bool Control::getPinning(){
    for(int i=0;i<3;i++){
        if(pindirect[i]!=0){
            return true;
        }
    }
    return false;
}

bool Control::isUniflow(){
    for(int i=0;i<3;i++){
        if(uniflow[i]!=0){
            return true;
        }
    }
    return false;
}

void Control::buildShearTensor(){
    sheartensor=Rvec(9,0.);
    if(uniflow[0]==1){
        sheartensor[0]=flowrate;
        sheartensor[4]=-flowrate/2;
        sheartensor[8]=-flowrate/2;
    }
    else if(uniflow[1]==1){
        sheartensor[0]=-flowrate/2;
        sheartensor[4]=flowrate;
        sheartensor[8]=flowrate;
    }
    else if(uniflow[2]==1){
        sheartensor[0]=flowrate;
        sheartensor[4]=flowrate;
        sheartensor[8]=-flowrate/2;
    }
    if(sheardir=0)
        sheartensor[1]=shearrate;
    else if(sheardir=1)
        sheartensor[2]=shearrate;
    else if(sheardir=2)
        sheartensor[3]=shearrate;
    else if(sheardir=3)
        sheartensor[5]=shearrate;
    else if(sheardir=4)
        sheartensor[6]=shearrate;
    else if(sheardir=5)
        sheartensor[7]=shearrate;


    return;
}
bool Control::isWallSheared(){
    if(wsheardir==-1)
        return false;
    else
        return true;
}

void Control::annealingTemperature(){
    temperature+=tempannrate*dt;
    return;
}
