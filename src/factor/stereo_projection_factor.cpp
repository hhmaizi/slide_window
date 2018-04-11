#include "stereo_projection_factor.hpp"
#include <eigen3/Eigen/Dense>

// Eigen::Matrix3d ProjectionFactor::sqrt_info;
double ProjectionFactor::sum_t;

ProjectionFactor::ProjectionFactor(const Eigen::Vector3d &_pts_i, const Eigen::Vector3d &_observation_j) : pts_i(_pts_i), observation(_observation_j)
{
};

ProjectionFactor::ProjectionFactor(const Eigen::Vector3d &_pts_i, const Eigen::Vector3d &_observation_j, const Eigen::Matrix3d &_info) : pts_i(_pts_i), observation(_observation_j)
{
	// calc sqrt_info using info, check this!
	 sqrt_info = Eigen::LLT<Eigen::Matrix<double, 3, 3>>(_info).matrixL().transpose();
};

bool ProjectionFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    TicToc tic_toc;
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

	// project pts_i from world frame to cam frame.
// 	Eigen::Vector3d Pbc = 	
#if 0
	std::cout<<"+++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
	std::cout<<"tic"<<tic.transpose()<<std::endl;
	std::cout<<"qic: "<<qic.vec().transpose()<<", "<<qic.w()<<std::endl;
	std::cout<<"ric: "<<qic.toRotationMatrix()<<std::endl;
	std::cout<<"+++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
#endif
	Eigen::Vector3d pts_cam_i = qic.inverse().toRotationMatrix() * \
								Qi.inverse().toRotationMatrix() *  \
								( pts_i - Pi) - qic.inverse().toRotationMatrix() * tic;

#if 0
	std::cout<<"=========================================================="<<std::endl;
	std::cout<<"Residual check in stereo projection factor: "<<std::endl;
	
	std::cout<<"Pi: ";
	for (int i=0; i<3; i++)
	{
		std::cout<<Pi[i]<<", ";
	}
	std::cout<<std::endl;
	
	std::cout<<"Qi: "<<Qi.w() <<", "<<Qi.x() << ", " << Qi.y()<<", "<<Qi.z()<< std::endl;
	
	std::cout<<"rci: "<<std::endl;
	std::cout<< qic.inverse().toRotationMatrix();
	std::cout<<std::endl;
	
	std::cout<<"tic: "<<tic.transpose()<<std::endl;
	
// 	std::cout<<parameters[0][6]<<std::endl;
// 		std::cout<<"================================================"<<std::endl;
	std::cout<<"Pw " << pts_i.transpose()<<std::endl;
	std::cout<<"Pc " << pts_cam_i.transpose()<<std::endl;
#endif 
    Eigen::Map<Eigen::Vector3d> residual(residuals);

	double invz = 1.0/pts_cam_i[2];
	double invz_2 = invz * invz;
	double x = pts_cam_i[0];
	double y = pts_cam_i[1];
	double z = pts_cam_i[2];

	// remember cx, cy. project map point to stereo cam
	double predicted_x, predicted_y, predicted_xr;
	predicted_x = fx * x * invz + cx;
	predicted_y = fy * y * invz + cy;
	predicted_xr =  predicted_x - bf * invz;

#if 0
	Eigen::Vector3d prediction(predicted_x, predicted_y, predicted_xr);
	std::cout<<"prediction " <<prediction.transpose() <<std::endl;
	std::cout<<"observation " << observation.transpose()<< std::endl;
#endif	
	residual[0] = predicted_x - observation[0];
	residual[1] = predicted_y - observation[1];
	residual[2] = predicted_xr - observation[2];
#if 0
	std::cout<<"residuals: "<<residual[0]<<","<<residual[1]<< "," <<residual[2]<<std::endl;
	std::cout<<"sqrt_info "<<sqrt_info<<std::endl;
#endif
    residual = sqrt_info * residual;
#if 0
	std::cout<<"residual sqrt_info: "<<residual.transpose()<<std::endl;
	std::cout<<"======================================================="<<std::endl;
#endif	
    if (jacobians)
    {		
		// more info about jacobians , see ceres tutorial "Analytic derivatives."
		// prepare your data for for jacobians.

		//jacobians[0] pose i
		//jacobians[1] pose j
		// jacobians[2] cam extrinsic(to implement)
		// jacobian of stereo cam Projection
		if (jacobians[0])
		{
			Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);
			jacobian_pose_i.setZero(); // right order for eigen map? damn eigen map!
			
			Eigen::Matrix<double, 3,3> stereoMaux;
			stereoMaux.setZero();
			stereoMaux(0,0) = fx;
			stereoMaux(0, 0) = fx;
			stereoMaux(0, 1) = 0;
			stereoMaux(0, 2) = -x / z * fx;
			stereoMaux(1, 0) = 0;
			stereoMaux(1, 1) = fy;
			stereoMaux(1, 2) = -y / z * fy;
			stereoMaux(2, 0) = fx;
			stereoMaux(2, 1) = 0;
			stereoMaux(2, 2) = -x / z * fx + bf / z;
			Eigen::Matrix<double, 3, 3>stereoJpi = -stereoMaux / z; //重投影误差对于相机坐标系下的点Pc的雅克比矩阵,注意：这里已经加了负号

		//stereo jacobian  of Pc with respect to PR
			Eigen::Matrix<double, 3, 3> Rcb = qic.inverse().toRotationMatrix();
			Eigen::Matrix<double, 3,3> Rwb = Qi.toRotationMatrix();
			Eigen::Vector3d Pw = pts_i;
	// 		Eigen::Vector3d Pwb = Pi;
			Eigen::Matrix<double, 3, 6> JPcPR ;
			JPcPR.setZero();
			JPcPR.block<3, 3>(0, 0) = -Rcb;

	// 		Rwb.transpose();
			Eigen::Vector3d Paux = Rcb * Rwb.transpose() * (Pw - Pi);
			JPcPR.block<3, 3>(0, 3) = Utility::skewSymmetric(Paux) * Rcb;
			//Map of jacobians
	// 		Eigen::Map<Eigen::Matrix<double, 6, 3, Eigen::RowMajor> > copyJacobian(jacobians);
	// 		copyJacobian.setZero();
			jacobian_pose_i.leftCols<6>() = sqrt_info * stereoMaux * JPcPR;
// 			jacobian_pose_i.leftCols<6>() = -stereoMaux * JPcPR;
// 			jacobian_pose_i.leftCols<6>() = stereoMaux * JPcPR;
// 			jacobian_pose_i.leftCols<6>() = -(sqrt_info * stereoMaux * JPcPR); // residual = predict - obs
			jacobian_pose_i.rightCols<1>().setZero();
// 			Eigen::MatrixXd copyJacobian = (stereoMaux * JPcPR).transpose();
		}
// 		if (jacobians[1])
// 		{
// 			std::cout<< "j1 " <<"should not be here."<<std::endl;
// 		}
// 		std::cout<< "jacobians: " <<jacobians[0][0]<<","<<jacobians[0][20]<<std::endl;
// 		check(parameters);
	}
    sum_t += tic_toc.toc();

    return true;
}
// check! check!!
// params, residual, and jacobians
void ProjectionFactor::check(double **parameters)
{	
    double *res = new double[3];// why 15? in accordance with imupreint factor?
    // jacobians dimension determined by **parameters.
    // we have only pose i, so dimension == 1.
    double **jaco = new double *[1]; 
    // dim_of_residual x dim_of_pose
	jaco[0] = new double[3*7];
	puts("in check.");
	
	for (int i = 0; i<7;i++)
	{
		std::cout<<"parameters: "<<parameters[0][i]<<",";
	}

	Evaluate(parameters, res, jaco);
	puts("check begins.");
	
	std::cout << "evaluate residual: " << Eigen::Map<Eigen::Matrix<double,3,1>>(res).transpose() << std::endl;
	
	std::cout << "evaluate jacobians: " << Eigen::Map<Eigen::Matrix<double,3,7>>(jaco[0]) << std::endl;
	std::cout<<"out of check" << std::endl;
	std::cout<<std::endl<<std::endl;
	
	// check the evaluate module--residual part
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

	// project pts_i from world frame to cam frame.
// 	Eigen::Vector3d Pbc = 	
	Eigen::Vector3d pts_cam_i = qic.inverse().toRotationMatrix() * \
								Qi.inverse().toRotationMatrix() *  \
								( pts_i - Pi) - qic.inverse().toRotationMatrix() * tic;
//     Eigen::Map<Eigen::Vector3d> residual(residuals);
	Eigen::Vector3d residual;

	double invz = 1.0/pts_cam_i[2];
	double invz_2 = invz * invz;
	double x = pts_cam_i[0];
	double y = pts_cam_i[1];
	double z = pts_cam_i[2];

	// remember cx, cy. project map point to stereo cam
	double predicted_x, predicted_y, predicted_xr;
	predicted_x = fx * x * invz + cx;
	predicted_y = fy * y * invz + cy;
	predicted_xr =  predicted_x - bf * invz;
	
	residual[0] = predicted_x - observation[0];
	residual[1] = predicted_y - observation[1];
	residual[2] = predicted_xr - observation[2];	
	puts("residual numeric: ");
	std::cout<< residual.transpose()<<std::endl;
	
	// check update using jacobian
	Eigen::Matrix<double, 3, 6> num_jacobian;
	double eps = 1e-6;
	for (int k = 0; k < 6; k++)
	{
		Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
		Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
				
		int a = k / 3, b = k % 3;
		Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;
		
		if (a == 0)
			Pi += delta;
		else if (a == 1)
			Qi = Qi * Utility::deltaQ(delta);
		
		// using updated Pi, Qi to calculate residual
		Eigen::Vector3d pts_cam_i = qic.inverse().toRotationMatrix() * \
								Qi.inverse().toRotationMatrix() *  \
								( pts_i - Pi) - qic.inverse().toRotationMatrix() * tic;
//     Eigen::Map<Eigen::Vector3d> residual(residuals);
		Eigen::Vector3d residual;

		double invz = 1.0/pts_cam_i[2];
		double invz_2 = invz * invz;
		double x = pts_cam_i[0];
		double y = pts_cam_i[1];
		double z = pts_cam_i[2];

		// remember cx, cy. project map point to stereo cam
		double predicted_x, predicted_y, predicted_xr;
		predicted_x = fx * x * invz + cx;
		predicted_y = fy * y * invz + cy;
		predicted_xr =  predicted_x - bf * invz;
		
		Eigen::Vector3d tmp_residual;
		tmp_residual << (predicted_x - observation[0]), (predicted_y - observation[1]), (predicted_xr - observation[2]);
		
// 		tmp_residual = sqrt_info * tmp_residual;
		num_jacobian.col(k) = (tmp_residual - residual) / eps;
	}
	// numerical jacobian
    std::cout << "numerical jacobian: " << num_jacobian << std::endl;
}
