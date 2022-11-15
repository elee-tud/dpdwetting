#include "newproperty.hpp"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace dpd;
/*General variables for property calcualtion*/
NewProperty::NewProperty(InitialSet initset):Property(initset){
    /*Title of the calculation*/
    title="   Tutorial for Calculation of new properties  ";

    nprops=2;               /*The number of properties to be calculated*/
    for_timeevol=true;      /*Is the property calculated for each time?*/
    for_distrib=false;      /*Is the property calculated by averaing whole time?*/
    need_position=true;     /*Do you need positions for the property?*/
    need_velocity=false;    /*Do you need velocities for the property?*/
    need_force=false;       /*Do you need forces for the property*/

    outfile="new_property.out";     /*Default output file name*/

    /*Header in the output file name(first line) */
    outtitle_tevol="#Polymer Radius of gyration and end-to-end vector"; 
    /*Keys in the output file (second line)*/
    outheader_tevol=Svec{"Time", "prop1", "prop2"};

}


/*Getting parameters specific for calculation from command line*/
void NewProperty::getSpecificParameters(){
    /*Obtaining a parameter from the command line, 1 is default
    command->getCommandSingleOption("-np", 1, &new_param);  

    Example)
    */
    command->getCommandSingleOption("-mol", "POL", &molname);
    command->getCommandSingleOption("-molidx", 1, &molidx);
    return;
}



/*Initializing the simulation variables*/
void NewProperty::initializeVariables(){
    /*Container of particle pointers to be used to calculate the property*/
    particles_to_calculate=ParticleGroup(pbc);


    for(int i=0;i<particles.size();i++){
        /*If molecule name is 'molname' and molecule index is 'molidx'*/
        if(particles[i]->getMoleculeName().compare(molname)==0 && particles[i]->getMoleculeIndex()==molidx){
            /*Add the particle to the container*/
            particles_to_calculate.addParticle(particles[i]);
        }
    }


    initializeResultArrays();
    return;
}


/*Calculating properties at the current step*/
void NewProperty::calculateStep(int step){
    /*
    tevol[0][step]=First property to calculate at the current step
    tevol[1][step]=Second property to calculate at the current step

    particle[i]->coord : 3d coordinate of the particle with the index i
    particle[i]->veloc : 3d velocity of the particle with the index i
    particle[i]->force : 3d force of the particle with the index i

    Example)
    */
    tevol[0][step]=particles_to_calculate.calculateRadiusOfGyration();



    return;
}







