/********************************************************
 *                  class NEWPROPERTY
 *
 *  This class contains an example code for implementing
 * new analysis algorithm. One should rename the class
 *********************************************************/


#ifndef __NEWPROPERTY__HPP
#define __NEWPROPERTY__HPP

#include "property.hpp"
#include "particlegroup.hpp"


namespace dpd{

class NewProperty:public Property{
private:
    /*example parameter from command line
    int new_param;
    */
    std::string molname;
    int molidx;

    /*Container of particle pointers to be used*/
    ParticleGroup particles_to_calculate;



    


public:
    
    NewProperty(){}
    NewProperty(InitialSet initset);
    ~NewProperty();

    /*Obtaining the parameters necessary for the calculation from the command line*/
    virtual void getSpecificParameters();
    /*Initializing the simulation variables*/
    virtual void initializeVariables();
    /*Calculating properties at the current step*/
    virtual void calculateStep(int step);
};
};



#endif
