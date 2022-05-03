#ifndef __INDEX__HPP
#define __INDEX__HPP

#include <fstream>
#include "types.hpp"
namespace dpd{
class Index{
private:
    std::string ndx_fname;
    std::ifstream ndxstream;


public:
    Index(){}
    Index(std::string ndx_fname);
    ~Index(){}

    void openFile();
    std::streampos seekArgument(std::string argument);
    Ivec getParticleIndexOfGroup(std::string argument);
    void exitErrorMissingArgument(std::streampos position, std::string argument);
};
};



#endif
