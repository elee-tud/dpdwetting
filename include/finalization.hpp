#ifndef __FINALIZATION__HPP
#define __FINALIZATION__HPP
#include "initialization.hpp"
#include "system.hpp"

namespace dpd{

class Finalization{
private:
    Initialization* init;
    System* system;
    Configuration* config;
    Control* control;
    Command* command;
    Dump* dump;
    Timer* timer;
    void printFinalInfo();
public:
    Finalization(){}
    Finalization(Initialization* init, System* sys);
    ~Finalization(){}
};
};

#endif
