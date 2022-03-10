#include "parsing.hpp"
#include <sstream>
using namespace dpd;
Svec dpd::parsing(std::string instring){
    std::string buf;
    Svec tokens;
    std::stringstream ss(instring);
    while(ss>>buf){
        tokens.push_back(buf);
    }
    return tokens;
}

