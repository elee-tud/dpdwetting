#include "index.hpp"
#include "parsing.hpp"

using namespace dpd;
Index::Index(std::string ndx_fname):ndx_fname(ndx_fname){

}

void Index::openFile(){
	try{
		ndxstream.open(ndx_fname);
		if(!ndxstream.is_open())
			throw ndx_fname;
	}catch(std::string ndx_fname){
        std::cout << "File name " << ndx_fname << " does not exist." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 15);
		exit(0);
	}
    return;
}


std::streampos Index::seekArgument(std::string argument){
    Svec tokens;
    ndxstream.clear();
    ndxstream.seekg(0, std::ios::beg);
    std::string line;
    std::streampos current_position=-1;
    while(std::getline(ndxstream, line)){
        tokens=dpd::parsing(line);
        if(line[0]!=';' && line.length()!=0 && !tokens.empty()){
            if(line.compare(argument)==0){
                current_position=ndxstream.tellg();
                break;
            }
        }

    }
    return current_position;
}

Ivec Index::getParticleIndexOfGroup(std::string group){
    Ivec index;
    std::string line;
    Svec tokens;
    openFile();
    std::string argument="[ "+group+" ]"; 
    std::streampos pos=seekArgument(argument);
    exitErrorMissingArgument(pos, argument);
    ndxstream.seekg(pos, std::ios::beg);
    while(std::getline(ndxstream, line)){
        tokens=dpd::parsing(line);
        if(line[0]=='[' || line.length()==0 || tokens.empty())
            break;
        if(line[0]!=';' && line.length()!=0){
            for(int i=0;i<tokens.size();i++){
                index.push_back(std::stoi(tokens[i])-1);
            }
        }
    }
    return index;
    ndxstream.close();
}




void Index::exitErrorMissingArgument(std::streampos position, std::string argument){
    if(position==-1){
        std::cout << "An group "<< argument << " should be given in the index file." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 16);
        exit(0);
    }
    return;
}

