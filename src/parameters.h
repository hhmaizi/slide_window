#ifndef 	PARAMETERS_H_
#define PARAMETERS_H_

// #include <ros/ros.h>
#include <vector>
// #include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Dense>
#include "utility/utility.h"
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
// #include <fstream>

const bool ESTIMATE_EXTRINSIC = false;

const double FOCAL_LENGTH = 460.0;
const int WINDOW_SIZE = 5;
const int NUM_OF_CAM = 1;
const int NUM_OF_F = 1000;
// 0 FOR UNIT_SPHERE_ERROR, 1 FOR MONO CAM, 2 FOR STEREO CAM
#define PROJECTION_ERROR_TYPE 2
// 9.08365697e-02, 9.25736523e+00, 3.24484825e+00
// const Eigen::Vector3d G(-9.08365697e-02, -9.25736523e+00, -3.24484825e+00);
// const Eigen::Vector3d G(9.08365697e-02, 9.25736523e+00, 3.24484825e+00);
// const Eigen::Vector3d G(6.07326210e-01, -9.08026218e+00, -3.66279817e+00);// V103
// -5.75561345e-01, 9.08304977e+00, 3.66101575e+00
// const Eigen::Vector3d G(5.75561345e-01, -9.08304977e+00, -3.66101575e+00);// V10302
extern Eigen::Vector3d G;// V10302

const double SOLVER_TIME = 30.0;

enum SIZE_PARAMETERIZATION
{
    SIZE_POSE = 7,
    SIZE_SPEEDBIAS = 9
//     SIZE_FEATURE = 3
};

const int NUM_ITERATIONS = 200;// 20 originally
// extrinsic cam params
extern Eigen::Matrix3d ric;
extern Eigen::Vector3d tic;
extern Eigen::Quaterniond qic;
// intrinsic cam params
extern double fx;
extern double fy;
extern double cx;
extern double cy;
extern double bf;

#endif