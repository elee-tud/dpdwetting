#include "topology.hpp"
#include <sstream>
#include <algorithm>
#include "parsing.hpp"
#include "mpiclasses.hpp"

using namespace dpd;
Topology::Topology(Command* command, Control *control, SetMPI* mpi):command(command), control(control), mpi(mpi){
    nmol_types=0;
    nbead_types=0;
    nbond_types=0;
    nonbonded=control->getNonbondedInteraction();
    bondlength=control->getBondLengthInteraction();
    topol_fname=command->topology();


}

void Topology::openFile(){
	try{
		topstream.open(topol_fname);
		if(!topstream.is_open())
			throw topol_fname;
	}catch(std::string topol_fname){
        std::cout << "File name " << topol_fname << " does not exist." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 15);
		exit(0);
	}
    return;
}

void Topology::readTopology(){
    if(mpi->isMaster()){
        openFile();
        readMoleculeTypes();
        readBeadTypes();
        readBondTypes();
        readAtomsInMolecule();
        readNonbondedParameters();
        if(nbond_types>0){
            readBondsInMolecule();
            readBondsBetweenMolecules();
            readBondLengthParameters();
        }
        closeFile();
    }
    return;
}



std::streampos Topology::seekArgument(std::string argument){
    topstream.clear();
    topstream.seekg(0, std::ios::beg);
    std::string line;
    std::streampos current_position=-1;
    while(std::getline(topstream, line)){
        tokens=dpd::parsing(line);
        if(line[0]!=';' && line.length()!=0 && !tokens.empty()){
            if(line.compare(argument)==0){
                current_position=topstream.tellg();
                break;
            }
        }

    }
    return current_position;
}

void Topology::readMoleculeTypes(){
    std::string argument="[ System ]";
    std::streampos position=seekArgument(argument);
    exitErrorMissingArgument(position, argument);
    topstream.seekg(position, std::ios::beg);
    while(std::getline(topstream, line)){
        tokens=dpd::parsing(line);
        if(line[0]=='[' || line.length()==0 || tokens.empty())
            break;
        if(line[0]!=';' && line.length()!=0){
            molecules.push_back(new Molecule(tokens[0], std::stoi(tokens[1])));
            if(tokens.size()>2){
                if(tokens[2].compare("frozen")==0)
                    molecules.back()->setFrozen();
            }
        }
    }

    nmol_types=molecules.size();

    return;
}

void Topology::readBeadTypes(){

    for(int i=0;i<nmol_types;i++){
        std::string argument="[ "+molecules[i]->getName()+" Atoms ]";
        std::streampos position=seekArgument(argument);
        exitErrorMissingArgument(position, argument);
        topstream.seekg(position, std::ios::beg);
        while(std::getline(topstream, line)){
//            if(line[0]=='[' || line.length()==0)
            tokens=dpd::parsing(line);
            if(line[0]=='[' || line.length()==0 || tokens.empty())
                break;
            if(line[0]!=';' && line.length()!=0){
                if(find(bead_types.begin(), bead_types.end(), tokens[1])==bead_types.end())
                    bead_types.push_back(tokens[1]);
            }
        }
    }
    nbead_types=bead_types.size();
    return;
}

void Topology::readBondTypes(){
    std::string argument="[ Bondlength ]";
    std::streampos position=seekArgument(argument);
    topstream.seekg(position, std::ios::beg);
    if(position==-1){
        if(control->getBondLengthInteraction()==-1){
            nbond_types=0;
            return;
        }
        else
            exitErrorMissingArgument(position, argument);
    }
    while(std::getline(topstream, line)){
        tokens=dpd::parsing(line);
        if(line[0]=='[' || line.length()==0 || tokens.empty())
//        if(line[0]=='[' || line.length()==0)
            break;
        if(line[0]!=';' && line.length()!=0){
            if(find(bond_types.begin(), bond_types.end(), tokens[0])==bond_types.end())
                bond_types.push_back(tokens[0]);
        }
    }
    nbond_types=bond_types.size();
    /*
    for(int i=0;i<nbond_types;i++)
        std::cout << bond_types[i] << std::endl;
        */

    return;
}

    




void Topology::readAtomsInMolecule(){

    for(int i=0;i<nmol_types;i++){
        std::string argument="[ "+molecules[i]->getName()+" Atoms ]";
        std::streampos position=seekArgument(argument);
        exitErrorMissingArgument(position, argument);
        topstream.seekg(position, std::ios::beg);
        while(std::getline(topstream, line)){
            tokens=dpd::parsing(line);
            if(line[0]=='[' || line.length()==0 || tokens.empty())
//            if(line[0]=='[' || line.length()==0)
                break;
            if(line[0]!=';' && line.length()!=0){
                molecules[i]->addAtoms(std::stoi(tokens[0]), tokens[1], distance(bead_types.begin(), find(bead_types.begin(), bead_types.end(), tokens[1])), std::stod(tokens[2]));
            }
        }
        /*
        std::cout << molecules[i].getName() << std::endl;
        for(int j=0;j<molecules[i].getNatoms();j++){
            std::cout << molecules[i].getAtomTypes()[j] << std::endl;
        }
        */
    }

    return;
}


void Topology::readBondsInMolecule(){
    for(int i=0;i<nmol_types;i++){
        if(molecules[i]->getAtoms().size()>1){
            std::string argument="[ "+molecules[i]->getName()+" Bonds ]";
            std::streampos position=seekArgument(argument);
            exitErrorMissingArgument(position, argument);
            if(position!=-1){
                topstream.seekg(position, std::ios::beg);
                while(std::getline(topstream, line)){
                    tokens=dpd::parsing(line);
                    if(line[0]=='[' || line.length()==0 || tokens.empty())
    //                if(line[0]=='[' || line.length()==0)
                        break;
                    if(line[0]!=';' && line.length()!=0){
                        if(bond_types.size()==0)
                            exitErrorMissingArgument(-1, "[ Bondlength ]");
                        molecules[i]->addBonds(std::stoi(tokens[0]), std::stoi(tokens[1]), distance(bond_types.begin(), find(bond_types.begin(), bond_types.end(), tokens[2])));
                    }
                }
            }
        }
    }
    /*
    for(int i=0;i<molecules.size();i++){
        std::cout << molecules[i]->getBonds() << std::endl;
    }
    */


}

void Topology::readBondsBetweenMolecules(){
    std::string argument="[ Interbonds ]";
    std::streampos position=seekArgument(argument);
    topstream.seekg(position, std::ios::beg);
    while(std::getline(topstream, line)){
//            if(line[0]=='[' || line.length()==0)
        tokens=dpd::parsing(line);
        if(line[0]=='[' || line.length()==0 || tokens.empty())
            break;
        if(line[0]!=';' && line.length()!=0){
            interbonds.push_back(Ivec{std::stoi(tokens[0]), std::stoi(tokens[1]), static_cast<int>(distance(bond_types.begin(), find(bond_types.begin(), bond_types.end(), tokens[2])))});
        }
    }

    return;
}

void Topology::readNonbondedParameters(){
    if(nonbonded==NBMDPD){
        A=Rvec2D(nbead_types, Rvec(nbead_types, 0.0));
        B=Rvec2D(nbead_types, Rvec(nbead_types, 0.0));
        Arcut=Rvec2D(nbead_types, Rvec(nbead_types, 0.0));
        Brcut=Rvec2D(nbead_types, Rvec(nbead_types, 0.0));
        std::string argument="[ Nonbonded ]";
        std::streampos position=seekArgument(argument);
        topstream.seekg(position, std::ios::beg);
        exitErrorMissingArgument(position, argument);
        while(std::getline(topstream, line)){
//            if(line[0]=='[' || line.length()==0)
            tokens=dpd::parsing(line);
            if(line[0]=='[' || line.length()==0 || tokens.empty())
                break;
            if(line[0]!=';' && line.length()!=0){
                int index1=distance(bead_types.begin(), find(bead_types.begin(), bead_types.end(), tokens[0]));
                int index2=distance(bead_types.begin(), find(bead_types.begin(), bead_types.end(), tokens[1]));
                B[index1][index2]=std::stod(tokens[2]);
                B[index2][index1]=std::stod(tokens[2]);
                Brcut[index1][index2]=std::stod(tokens[3]);
                Brcut[index2][index1]=std::stod(tokens[3]);
                A[index1][index2]=std::stod(tokens[4]);
                A[index2][index1]=std::stod(tokens[4]);
                Arcut[index1][index2]=std::stod(tokens[5]);
                Arcut[index2][index1]=std::stod(tokens[5]);
            }
        }
    }
    else if(nonbonded==NBDPD){
        B=Rvec2D(nbead_types, Rvec(nbead_types, 0.0));
        Brcut=Rvec2D(nbead_types, Rvec(nbead_types, 0.0));
        std::string argument="[ Nonbonded ]";
        std::streampos position=seekArgument(argument);
        topstream.seekg(position, std::ios::beg);
        exitErrorMissingArgument(position, argument);
        while(std::getline(topstream, line)){
//            if(line[0]=='[' || line.length()==0)
            tokens=dpd::parsing(line);
            if(line[0]=='[' || line.length()==0 || tokens.empty())
                break;
            if(line[0]!=';' && line.length()!=0){
                int index1=distance(bead_types.begin(), find(bead_types.begin(), bead_types.end(), tokens[0]));
                int index2=distance(bead_types.begin(), find(bead_types.begin(), bead_types.end(), tokens[1]));
                B[index1][index2]=std::stod(tokens[2]);
                B[index2][index1]=std::stod(tokens[2]);
                Brcut[index1][index2]=std::stod(tokens[3]);
                Brcut[index2][index1]=std::stod(tokens[3]);
            }
        }
    }


    return;
}

void Topology::readBondLengthParameters(){
    if(bondlength==BLHARMONIC){
        bondk=Rvec(nbond_types, 0.0);
        bondl=Rvec(nbond_types, 0.0);
        std::string argument="[ Bondlength ]";
        std::streampos position=seekArgument(argument);
        topstream.seekg(position, std::ios::beg);
        exitErrorMissingArgument(position, argument);
        while(std::getline(topstream, line)){
//            if(line[0]=='[' || line.length()==0)
            tokens=dpd::parsing(line);
            if(line[0]=='[' || line.length()==0 || tokens.empty())
                break;
            if(line[0]!=';' && line.length()!=0){
                int index=distance(bond_types.begin(), find(bond_types.begin(), bond_types.end(), tokens[0]));
                bondk[index]=std::stod(tokens[1]);
                if(tokens.size()>2 && tokens[2][0]!=';'){
                    bondl[index]=std::stod(tokens[2]);
                }
            }
        }
    }

    return;
}

void Topology::exitErrorMissingArgument(std::streampos position, std::string argument){
    if(position==-1){
        std::cout << "An argument "<< argument << " should be given in the topology file." << std::endl;
        exit(0);
    }
    return;
}


void Topology::bcastTopology(){
    MPI_Bcast(&nbead_types, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&nbond_types, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    if(nonbonded==NBMDPD){
        if(!mpi->isMaster()){
            A=Rvec2D(nbead_types, Rvec(nbead_types, 0.));
            Arcut=Rvec2D(nbead_types, Rvec(nbead_types, 0.));
            B=Rvec2D(nbead_types, Rvec(nbead_types, 0.));
            Brcut=Rvec2D(nbead_types, Rvec(nbead_types, 0.));

        }
        for(int i=0;i<nbead_types;i++){
            MPI_Bcast(&A[i][0], nbead_types, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
            MPI_Bcast(&Arcut[i][0], nbead_types, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
            MPI_Bcast(&B[i][0], nbead_types, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
            MPI_Bcast(&Brcut[i][0], nbead_types, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
        }
    }
    else if(nonbonded==NBDPD){
        if(!mpi->isMaster()){
            B=Rvec2D(nbead_types, Rvec(nbead_types, 0.));
            Brcut=Rvec2D(nbead_types, Rvec(nbead_types, 0.));

        }
        for(int i=0;i<nbead_types;i++){
            MPI_Bcast(&B[i][0], nbead_types, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
            MPI_Bcast(&Brcut[i][0], nbead_types, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
        }
    }

    if(bondlength==BLHARMONIC){
        if(!mpi->isMaster()){
            bondk=Rvec(nbond_types, 0.);
            bondl=Rvec(nbond_types, 0.);
        }
        MPI_Bcast(&bondk[0], nbond_types, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
        MPI_Bcast(&bondl[0], nbond_types, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
        int numinterbonds=interbonds.size();
        MPI_Bcast(&numinterbonds, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
        if(mpi->rank()!=MASTER)
            interbonds=Ivec2D(numinterbonds, Ivec(3, 0));
        for(int i=0;i<numinterbonds;i++){
            MPI_Bcast(&interbonds[i][0], 3, MPI_INT, MASTER, MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    return;
}

void Topology::bcastMolecules(){
    MPI_Bcast(&nmol_types, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    molecules.reserve(nmol_types);
    if(!mpi->isMaster()){
        for(int i=0;i<nmol_types;i++){
            molecules.push_back(new Molecule());
        }
    }
    MPIClasses mpimolecule(mpi);
    for(int i=0;i<nmol_types;i++){
        mpimolecule.bcastMolecule(MASTER, molecules[i]);
    }
    /*
    if(mpi->rank()==10){
        std::cout << molecules[1]->getBonds() << std::endl;
    }
    */
    MPI_Barrier(MPI_COMM_WORLD);
    return;
}




void Topology::buildTopology(){

    nbeads=0;
    nbeads_unfrozen=0;
    nbonds=0;
    totmass=0.;
    int beadnum=1;
    int molnum=1;
    int index=0;
    for(int i=0;i<nmol_types;i++)
        nbeads+=molecules[i]->getNmolecules()*molecules[i]->getNatoms();

    for(int i=0;i<nmol_types;i++){
        for(int j=0;j<molecules[i]->getNmolecules();j++){
            for(int k=0;k<molecules[i]->getNatoms();k++){
                particles.push_back(new Particle(beadnum, molnum, molecules[i]->getAtomTypes()[k], molecules[i]->getMasses()[k], molecules[i]->isFrozen()));
                particles.back()->setNone();
                if(mpi->isMaster()){
                    particles.back()->setMoleculeString(molecules[i]->getName());
                    particles.back()->setTypeString(molecules[i]->getAtomTypeStrings()[k]);
                }

                beadnum++;
                totmass+=molecules[i]->getMasses()[k];
            }
            Ivec2D bonds=molecules[i]->getBonds();
            Ivec types=molecules[i]->getBondTypes();
            for(int k=0;k<molecules[i]->getNbonds();k++){
                particles[index+bonds[k][0]]->addBond(particles[index+bonds[k][1]], types[k]);
                particles[index+bonds[k][1]]->addBond(particles[index+bonds[k][0]], types[k]);
            }
            molnum++;
            nbonds+=molecules[i]->getNbonds();
            if(!molecules[i]->isFrozen())
                nbeads_unfrozen+=molecules[i]->getNatoms();
            index+=molecules[i]->getNatoms();
        }
    }

    /*Add intermolecular bonds*/
    for(int i=0;i<interbonds.size();i++){
        particles[interbonds[i][0]]->addBond(particles[interbonds[i][1]], interbonds[i][2]);
        particles[interbonds[i][1]]->addBond(particles[interbonds[i][0]], interbonds[i][2]);
        nbonds++;
        /*
        ParticleList bds=particles[interbonds[i][0]]->getBonds();
        for(int j=0;j<bds.size();j++){
            std::cout << particles[interbonds[i][0]]->getParticleIndex() <<":" << bds[j]->getParticleIndex()<< std::endl;
        }
        */

    }
    /*
    if(mpi->rank()==0){
        for(int i=0;i<particles.size();i++){
            std::cout << "myid=" << particles[i]->getParticleIndex() << ":";
            ParticleList bonded=particles[i]->getBonds();
            for(int j=0;j<bonded.size();j++){
                std::cout << bonded[j]->getParticleIndex() << " " ;
            }
            std::cout << std::endl;
        }
    }
    */

    /*
    if(mpi->rank()==1)

        std::cout << &particles[0]->coord[0] << std::endl;
        */

    return;
}

void Topology::printSystemInformation(){
    if(mpi->isMaster()){
        std::cout << " " << nmol_types << " molecule types " << std::endl;
        std::cout << " " << nbead_types << " bead types " << std::endl;;
        std::cout << " " << nbond_types << " bond types "  << std::endl;;
        std::cout << " " << nbeads << " total particles "  << std::endl;;
        if(nbonds!=0)
            std::cout << " " << nbonds << " total bonds " << std::endl;;
        if(nbeads_unfrozen!=nbeads)
            std::cout << " " << nbeads-nbeads_unfrozen << " frozen particless ."  << std::endl;;
        std::cout << std:: endl;
    }

    return; 
}

void Topology::checkDependency(){
    if(bondlength>0 && nbond_types==0){
        if(mpi->isMaster())
            std::cout << "[Error] Bond length parameters has to be given in the topology file\n" << std::endl;
        mpi->finalize();
        exit(0);
    }

    if(bondlength<0 && nbond_types>0){
        if(mpi->isMaster())
            std::cout << "[Error] A type of bond length energy has to be given in the control file with \"bondlength\" option.\n" << std::endl;
        mpi->finalize();
        exit(0);
    }
    return;
}



