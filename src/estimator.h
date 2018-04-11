#ifndef ESTIMATOR_H_
#define ESTIMATOR_H_

#include<fstream>
#include "parameters.h"

#include <iostream>
#include <string>
#include <vector>

#include "utility/utility.h"
#include "utility/tic_toc.h"
#include "utility/Converter.h"

#include <ceres/ceres.h>
// #include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/Dense>
#include "factor/imu_factor.h"
#include "factor/pose_local_parameterization.h"
#include "factor/stereo_projection_factor.hpp"

#include <opencv2/core/core.hpp>
#include <unordered_map>
#include <queue>
#include <opencv2/core/eigen.hpp>

// using namespace Eigen;

class Estimator
{
  public:
    Estimator::Estimator();

    void setParameter();

    // interface
	bool readIMUPreint();
	// try to use txt instead of yaml to read info
	bool readIMUPreint(char* imuFileName);
	bool readFeatures();
	bool readFeatures(char* prjFileName);

	// internal
    void clearState();
	void optimization_simple();
    void vector2double();
    void double2vector();

    enum SolverFlag
    {
        INITIAL,
        NON_LINEAR
    };

// 	std::stirng imuFileName = "";
// 	std::string featureFileName = "";
    SolverFlag solver_flag;
    Eigen::Vector3d g;

    Eigen::Matrix3d ric;
    Eigen::Vector3d tic;

	// howto assign?
    Eigen::Vector3d acc_0, gyr_0;
//     Eigen::Vector3d Ps[(WINDOW_SIZE + 1)];
//     Eigen::Vector3d Vs[(WINDOW_SIZE + 1)];
//     Eigen::Matrix3d Rs[(WINDOW_SIZE + 1)];
    std::vector<Eigen::Vector3d> Ps;
    std::vector<Eigen::Vector3d> Vs;
	std::vector<Eigen::Quaterniond> Qs;
    std::vector<Eigen::Matrix3d> Rs;
	std::vector<Eigen::Vector3d> Bas;
	std::vector<Eigen::Vector3d> Bgs;

    std::vector<double> dt_buf[(WINDOW_SIZE + 1)];
    std::vector<Eigen::Vector3d> linear_acceleration_buf[(WINDOW_SIZE + 1)];
    std::vector<Eigen::Vector3d> angular_velocity_buf[(WINDOW_SIZE + 1)];

	std::vector<int> frame_ids;
    std::vector<Eigen::Vector3d> point_cloud;// using map instead
//     std::map<int, Eigen::Vector3d> point_cloud;
	std::vector<Eigen::Vector3d> features; // using map instead
// 	std::map<int, Eigen::Vector3d> features;
	std::vector<Eigen::Matrix3d> stereo_infos;
// 	std::map<int, Eigen::Matrix3d> stereo_infos;
	double initial_timestamp;
	
// using eigen matrix to store all your data.
#if 1
	// for imu factor
	std::vector<double> dTs;
    std::vector<Eigen::Vector3d> dPs;
	std::vector<Eigen::Quaterniond> dQs;
	std::vector<Eigen::Vector3d> dVs;
	std::vector<Eigen::Matrix<double, 9,9>> CovsPVPhis;
	std::vector<Eigen::Matrix<double, 15,15>> CovsPVPhisBag;
	Eigen::Vector3d acc_linearized;
	Eigen::Vector3d gyr_linearized;	
		
// 	std::vector<Eigen::Matrix<double, 7,1> > Poses;
//     std::vector<Eigen::Matrix<double, 9,1> > SpeedBiases;
// 	// for projection factor
//     std::vector<Eigen::Matrix<double, 3,1> > Features;
// 	std::vector<Eigen::Matrix<double, 3,1> > MapPoints;
	
	// jacobian
// 	std::vector<Eigen::Vector3d> linearized_bas;
// 	std::vector<Eigen::Vector3d> linearized_bgs;
	std::vector<Eigen::Vector3d> dBas;
	std::vector<Eigen::Vector3d> dBgs;
	// delta measurements for imu factor
	std::vector<Eigen::Matrix<double, 3,3> > dp_dbas;
	std::vector<Eigen::Matrix<double, 3,3> > dp_dbgs;
	std::vector<Eigen::Matrix<double, 3,3> > dq_dbgs;
	std::vector<Eigen::Matrix<double, 3,3> > dv_dbas;
	std::vector<Eigen::Matrix<double, 3,3> > dv_dbgs;
#endif

	double para_Pose[WINDOW_SIZE + 1][SIZE_POSE];
    double para_SpeedBias[WINDOW_SIZE + 1][SIZE_SPEEDBIAS];
//     double para_Feature[NUM_OF_F][SIZE_FEATURE];
//     double para_Ex_Pose[NUM_OF_CAM][SIZE_POSE];
#if 0
	double deltaT(0.0);
    double dPs;
	double dQs;
	double dVs;
	double acc_linearized;
	double gyr_linearized;	
	
	double para_Pose[WINDOW_SIZE + 1][SIZE_POSE];
    double para_SpeedBias[WINDOW_SIZE + 1][SIZE_SPEEDBIAS];
    double para_Feature[NUM_OF_F][SIZE_FEATURE];
    double para_Ex_Pose[NUM_OF_CAM][SIZE_POSE];
    double para_Retrive_Pose[SIZE_POSE];
	
	// jacobian
	
	// delta measurements
	double dp_dbas[WINDOW_SIZE][9];
	double dp_dbgs[WINDOW_SIZE][9];
	double dq_dbgs[WINDOW_SIZE][9];
	double dv_dbas[WINDOW_SIZE][9];
	double dv_dbgs[WINDOW_SIZE][9];
#endif

};

#endif