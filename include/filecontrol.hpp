#ifndef __FILECONTROL__HPP
#define __FILECONTROL__HPP

#include <string>
#include <fstream>
namespace dpd{
    inline bool fileExists (const std::string& fileName); 
    void copyFile(std::string srce_file, std::string dest_file);
    void openFileWithBackup(std::string filename, std::ofstream *stream, bool binary, bool append);
};

#endif
