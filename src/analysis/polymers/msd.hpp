#ifndef __POLMSD__HPP
#define __POLMSD__HPP

#include "../property.hpp"
#include "../particlegroup.hpp"
namespace dpd{

class MeanSquareDisplacement:public Property{
private:
    MoleculeGroup polymers;
    int num_pol;
    ParticleList ref;
    R3vec2D mon_pos;
    R3vec2D com_pos;
    R3vec2D mon_com;
    std::string molname;



    


public:
    MeanSquareDisplacement(){}
    MeanSquareDisplacement(InitialSet initset);
    ~MeanSquareDisplacement();

    virtual void getSpecificParameters();
    virtual void initializeVariables();
    virtual void calculateStep(int step);
    void calculateDynamicProperty();
    void unfoldTrajectory();
};
};



#endif
