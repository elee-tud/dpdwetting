#ifndef __ELLIPSE_FIT_HPP
#define __ELLIPSE_FIT_HPP



#include "../../dep/Eigen/Dense"
#include "../../dep/Eigen/Core"
#include <iostream>
#include <vector>
using namespace std;
using namespace Eigen;
namespace dpd{

class ellipse_fit{


public:

    void set(vector<vector<double> > input);
    void fit(double& result_center_x, double& result_center_y, double& result_phi, double& result_width, double& result_hight);

private:

    vector<vector<double> >input_matrix;
    void take_input(vector<vector<double> > input);

};
};
#endif


