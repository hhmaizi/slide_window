#include "parameters.h"
#include <boost/concept_check.hpp>

// cam to imu extrinsic params
Eigen::Matrix3d ric;
Eigen::Vector3d tic(-2.16401462e-02, -6.46769851e-02, 9.81073081e-03); // weirdo tic didn't init in setParameter
Eigen::Quaterniond qic;
Eigen::Vector3d G(5.75561345e-01, -9.08304977e+00, -3.66101575e+00);// V10302

// left cam intrinsic param, after undistort
double fx(4.35204681e+02);
double fy(4.35204681e+02);
double cx(3.67451721e+02);
double cy(2.52200851e+02);
double bf(4.7906394e+01);