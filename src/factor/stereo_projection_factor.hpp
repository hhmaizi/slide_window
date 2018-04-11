#ifndef STEREO_PROJECTION_FACTOR_H
#define STEREO_PROJECTION_FACTOR_H

// #include <ros/assert.h>
#include <ceres/ceres.h>
#include <eigen3/Eigen/Dense> // using find eigen in cmakelist to settle this.
#include "./utility.h"
#include "./tic_toc.h"
#include "./parameters.h"

// 2: residual dimensions, 7: pose_i, 7:pose_j
// 7: extrinsic params, cam2imu
// 1: feature_inverse_depth.
// for stereo preojection: 3-> residual dim, 7-> pose_i, 3->mappoint(w/o?)
class ProjectionFactor : public ceres::SizedCostFunction<3, 7>
{
  public:
    ProjectionFactor(const Eigen::Vector3d &_pts_i, const Eigen::Vector3d & _observation_j);
	ProjectionFactor(const Eigen::Vector3d &_pts_i, const Eigen::Vector3d & _observation_j, \
					 const Eigen::Matrix3d &_info
	);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters);

    Eigen::Vector3d pts_i;
	Eigen::Vector3d observation;
	Eigen::Matrix3d sqrt_info;
//     Eigen::Matrix<double, 2, 3> tangent_base;
//     static Eigen::Matrix3d sqrt_info;
    static double sum_t;
// 	static Eigen::Matrix<double, 3,3>
};

#endif
