#include "stereo_projection_factor.hpp"

Eigen::Matrix3d ProjectionFactor::sqrt_info;
double ProjectionFactor::sum_t;

ProjectionFactor::ProjectionFactor(const Eigen::Vector3d &_pts_i, const Eigen::Vector3d &_observation_j) : pts_i(_pts_i), observation(_observation_j)
{
};

bool ProjectionFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    TicToc tic_toc;
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
//     Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
//     Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
	// project pts_i from world frame to cam frame.
// 	Eigen::Vector3d Pbc = 	
	Eigen::Vector3d pts_cam_i = qic.inverse().toRotationMatrix() * \
								Qi.inverse().toRotationMatrix() *  \
								( pts_i - Pi) - qic.inverse().toRotationMatrix() * tic;
/*
        Eigen::Matrix3d Ri = Qi.toRotationMatrix();
        Eigen::Matrix3d Rj = Qj.toRotationMatrix();*/
//         Eigen::Matrix3d ric = qic.toRotationMatrix();
// infomation matrix?
//     Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
//     Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
//     Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
//     Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
//     Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);

    Eigen::Map<Eigen::Vector3d> residual(residuals);

// 	double bf(30.0);
// 	double focal(300.0);
// 	double fx(300.0); double fy(300.0);
// 	double cx(150.0), cy(159.0);
	double invz = 1.0/pts_cam_i[2];
	double invz_2 = invz * invz;
	double x = pts_cam_i[0];
	double y = pts_cam_i[1];

	// remember cx, cy. project map point to stereo cam
	double predicted_x, predicted_y, predicted_xr;
	predicted_x = fx * x * invz + cx;
	predicted_y = fy * y * invz + cy;
	predicted_xr =  predicted_x - bf * invz;
	
	residual[0] = predicted_x - observation[0];
	residual[1] = predicted_y - observation[1];
	residual[2] = predicted_xr - observation[2];

    residual = sqrt_info * residual;

    if (jacobians)
    {
		// more info about jacobians , see ceres tutorial "Analytic derivatives."
		// prepare your data for for jacobians.

		// stereo pose only projection from orbslam2 EdgeStereoSe3Poseonly,
		// why its' derivative's sign is opposite to my derivation?
		//, we can try stereo projection and paul's stereo projection model later(3.100 phd thesis pp40).
// concerning position(local frame) tx, ty, tz
// tx
//   _jacobianOplusXi(0,3) = -invz *fx;
//   _jacobianOplusXi(1,3) = 0;
//   _jacobianOplusXi(2,3) = _jacobianOplusXi(0,3);
// ty
//   _jacobianOplusXi(0,4) = 0;
//   _jacobianOplusXi(1,4) = -invz *fy;
//   _jacobianOplusXi(2,4) = 0;
// tz
//   _jacobianOplusXi(0,5) = x*invz_2 *fx;
//   _jacobianOplusXi(1,5) = y*invz_2 *fy;
//   _jacobianOplusXi(2,5) = _jacobianOplusXi(0,5)-bf*invz_2;
//
// concerning wx
//   _jacobianOplusXi(0,0) =  x*y*invz_2 *fx;
//   _jacobianOplusXi(1,0) = (1+y*y*invz_2) *fy;
//   _jacobianOplusXi(2,0) = _jacobianOplusXi(0,0)-bf*y*invz_2;
// concerning wy
//   _jacobianOplusXi(0,1) = -(1+(x*x*invz_2)) *fx;
//   _jacobianOplusXi(1,1) = -x*y*invz_2 *fy;
//   _jacobianOplusXi(2,1) = _jacobianOplusXi(0,1)+bf*x*invz_2;
// concerning wz
//   _jacobianOplusXi(0,2) = y*invz *fx;
//   _jacobianOplusXi(1,2) = -x*invz *fy;
//   _jacobianOplusXi(2,2) = _jacobianOplusXi(0,2);
		// derivative concerning tx, ty, tz, slow assignment use eigen map instead.
        if (jacobians[0])
        {
			jacobians[0][0] = -invz *fx;
			jacobians[0][1] = 0.0;
			jacobians[0][2] = -invz *fx;
        }
        if (jacobians[1])
        {
			jacobians[1][0] = 0.0;
			jacobians[1][1] = -invz *fy;
			jacobians[1][2] = 0.0;
        }
        if (jacobians[2])
        {
			double tmp(x*invz_2 *fx);
			jacobians[2][0] = tmp;
			jacobians[2][1] = y*invz_2 *fy;
			jacobians[2][2] = tmp - bf * invz_2;
        }
        // derivative concerning omegax, omegay, omegaz
        if (jacobians[3])
        {
			jacobians[3][0] = x*y*invz_2 *fx;
			jacobians[3][1] = (1+y*y*invz_2) *fy;
			jacobians[3][2] = jacobians[1][0] - bf*y*invz_2;
		}
		if (jacobians[4])
		{
			jacobians[4][0] = -(1+(x*x*invz_2)) *fx;
			jacobians[4][1] = -x*y*invz_2 *fy;
			jacobians[4][2] = jacobians[4][0] + bf*x*invz_2;
		}
		if (jacobians[5])
		{
			jacobians[5][0] =  y*invz *fx;
			jacobians[5][1] = -x*invz *fy;
			jacobians[5][2] = jacobians[5][0];
		}
	}
    sum_t += tic_toc.toc();

    return true;
}

void ProjectionFactor::check(double **parameters)
{
//     double *res = new double[15];
//     double **jaco = new double *[4];
//     std::cout << num_jacobian << std::endl;
}
