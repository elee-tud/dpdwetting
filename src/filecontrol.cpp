#include "filecontrol.hpp"
#include <sys/stat.h>
#include <sstream>
#include <iostream>
namespace dpd{

inline bool fileExists (const std::string& fileName) {
    struct stat buff;  
    return (stat (fileName.c_str(), &buff) == 0);
}

void copyFile(std::string srce_file, std::string dest_file){
    std::ifstream srce(srce_file, std::ios::binary);
    std::ofstream dest(dest_file, std::ios::binary);
    dest << srce.rdbuf() ;
}


void openFileWithBackup(std::string filename, std::ofstream *stream, bool binary, bool append){
    if(append){
        if(binary){
            stream->open(filename, std::ios::out | std::ios::binary | std::ios::app);
        }
        else{
            stream->open(filename, std::ios::out | std::ios::app);
        }
    }
    else{
        if(fileExists(filename)){
            std::stringstream newfile;
            int index=1;
            while(true){
                newfile.str("");
                newfile << "#" << filename << "." << index << "#";
                if(!fileExists(newfile.str()))
                    break;
                index++;
            }
            std::rename(filename.c_str(), newfile.str().c_str());
        }
        if(binary){
            stream->open(filename, std::ios::out | std::ios::binary);
        }
        else{
            stream->open(filename);
        }
    }
    return;
}

}
