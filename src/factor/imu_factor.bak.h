#ifndef IMU_FACTOR_H
#define IMU_FACTOR_H

#include <iostream>
#include <iomanip> 
#include <eigen3/Eigen/Dense>

#include "./utility.h"
#include "./parameters.h"

#include <ceres/ceres.h>
#include "../IMU/IMUPreintegrator.h"
// 15: dimension for residual pvr ba bg
// 7: pose i dimension;
// 9: pose i velocity, bias of acc and bias of gyro
// 7: pose j; 9 pose j ...
// using wangjing's states definition and state propogation equation
#define NEGJACOB 0
#define CRAZY 0

using namespace ORB_SLAM2;

class IMUFactor : public ceres::SizedCostFunction<15, 7, 9, 7, 9>
{
  public:
    IMUFactor() = delete;

	IMUFactor(const IMUPreintegrator* imupreint_, Eigen::Vector3d linearized_ba_, Eigen::Vector3d linearized_bg_): imupreint(imupreint_), linearized_ba(linearized_ba_), linearized_bg(linearized_bg_)
	{
// 		imupreint(imupreint_);
// 		linearized_ba(linearized_ba_);
// 		linearized_bg(linearized_bg_);
		sum_dt = imupreint->getDeltaTime();
// vins states order: P R V Ba Bg
// viorb states order: P V R Ba Bg
		Eigen::Matrix<double, 9,9> cov_ = imupreint->getCovPVPhi();
		covariance.setZero();
		covariance.block<3,3>(0,0) = cov_.block<3,3>(0,0);
		covariance.block<3,3>(0,3) = cov_.block<3,3>(0,6);
		covariance.block<3,3>(0,6) = cov_.block<3,3>(0,3);
		covariance.block<3,3>(3,0) = cov_.block<3,3>(6,0);
		covariance.block<3,3>(6,0) = cov_.block<3,3>(3,0);
		
		covariance.block<3,3>(3,3) = cov_.block<3,3>(6,6);
		covariance.block<3,3>(3,6) = cov_.block<3,3>(6,3);
		covariance.block<3,3>(6,3) = cov_.block<3,3>(3,6);
		
		covariance.block<3,3>(6,6) = cov_.block<3,3>(3,3);
		
// 		covariance.block<3,3>(9, 9) = IMUData::getAccBiasRWCov();
// 		covariance.block<3,3>(12,12) = IMUData::getGyrBiasRWCov();
		covariance.block<3,3>(9, 9) = IMUData::getAccBiasRWCov() * sum_dt;
		covariance.block<3,3>(12,12) = IMUData::getGyrBiasRWCov() * sum_dt;
// 		covariance.setIdentity();
// 		int idx = 0;
// double IMUData::_gyrBiasRw2 = 2.0e-5*2.0e-5;  // sigma_gw * sigma_gw * dt
// double IMUData::_accBiasRw2 = 5.0e-3*5.0e-3;  // sigma_aw * sigma_aw * dt
		
		Eigen::IOFormat CleanFmt(2, 0, ", ", "\n", "[", "]");
#if 0
// 		std::cout<<"cov in imu factor: "<<std::endl;
// 		std::cout << covariance.block<9,9>(0,0).format(CleanFmt) << std::endl;
// 		std::cout<<"cov_: "<<std::endl;
// 		std::cout << cov_.format(CleanFmt) << std::endl;
		
		std::cout << "covariance diagonal: "<<std::endl;
		std::cout << std::setprecision(4)<< covariance.diagonal().format(CommaFmt)<<std::endl;
		std::cout << "cov_ diagonal: "<<std::endl;
		std::cout << std::setprecision(4) << cov_.diagonal().format(CommaFmt)<<std::endl;
		
		Eigen::Matrix<double, 15, 15> sqrt_infomat = Eigen::LLT<Eigen::Matrix<double, 15, 15> >(covariance.inverse()).matrixL().transpose();
// 		std::cout << "sqrt_info: "<<std::endl;
// 		std::cout << std::setprecision(4) << sqrt_infomat.block<9,9>(0,0).format(CleanFmt) << std::endl;
		std::cout << std::setprecision(4) << sqrt_infomat.block<6,6>(9,9).format(CleanFmt) << std::endl;
#endif
}

    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
    {

        Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
        Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

        Eigen::Vector3d Vi(parameters[1][0], parameters[1][1], parameters[1][2]);
        Eigen::Vector3d Bai(parameters[1][3], parameters[1][4], parameters[1][5]);
        Eigen::Vector3d Bgi(parameters[1][6], parameters[1][7], parameters[1][8]);

        Eigen::Vector3d Pj(parameters[2][0], parameters[2][1], parameters[2][2]);
        Eigen::Quaterniond Qj(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

        Eigen::Vector3d Vj(parameters[3][0], parameters[3][1], parameters[3][2]);
        Eigen::Vector3d Baj(parameters[3][3], parameters[3][4], parameters[3][5]);
        Eigen::Vector3d Bgj(parameters[3][6], parameters[3][7], parameters[3][8]);

//Eigen::Matrix<double, 15, 15> Fd; Eigen::Matrix<double, 15, 12> Gd;
//Eigen::Vector3d pPj = Pi + Vi * sum_t - 0.5 * g * sum_t * sum_t + corrected_delta_p;
//Eigen::Quaterniond pQj = Qi * delta_q;
//Eigen::Vector3d pVj = Vi - g * sum_t + corrected_delta_v;
//Eigen::Vector3d pBaj = Bai; Eigen::Vector3d pBgj = Bgi;
//Vi + Qi * delta_v - g * sum_dt = Vj;
//Qi * delta_q = Qj;

//delta_p = Qi.inverse() * (0.5 * g * sum_dt * sum_dt + Pj - Pi);
//delta_v = Qi.inverse() * (g * sum_dt + Vj - Vi);
//delta_q = Qi.inverse() * Qj;

        Eigen::Map<Eigen::Matrix<double, 15, 1> > residual(residuals);
        residual = Evaluate_Direct(Pi, Qi, Vi, Bai, Bgi, Pj, Qj, Vj, Baj, Bgj);

// 		std::cout<<"residual:"<<residual.transpose()<<std::endl;
// //         Eigen::Matrix<double, 15, 15> sqrt_info = Eigen::LLT<Eigen::Matrix<double, 15, 15> >(covariance.inverse()).matrixL();
        Eigen::Matrix<double, 15, 15> sqrt_info = Eigen::LLT<Eigen::Matrix<double, 15, 15> >(covariance.inverse()).matrixL().transpose();
// 		Eigen::Matrix<double, 15, 15> invCovariance = covariance.inverse();
//         Eigen::Matrix<double, 15, 15> sqrt_info = invCovariance.llt().matrixL();
// 		std::cout << "sqrt_info: "<<std::endl;
// 		std::cout << std::setprecision(4) << sqrt_info.block<9,9>(0,0).format(CleanFmt) << std::endl;
// 		std::cout << std::setprecision(4) << sqrt_info.block<6,>(9,9).format(CleanFmt) << std::endl;
// 		sqrt_info.setIdentity();
        residual = sqrt_info * residual;

        if (jacobians)
        {
			Eigen::Vector3d delta_p = imupreint->getDeltaP();
			Eigen::Quaterniond delta_q = imupreint->getDeltaQ();
			Eigen::Vector3d delta_v = imupreint->getDeltaV();
			
			Eigen::Matrix3d dp_dba = imupreint->getJPBiasa();
			Eigen::Matrix3d dp_dbg = imupreint->getJPBiasg();
			Eigen::Matrix3d dq_dbg = imupreint->getJRBiasg();
			Eigen::Matrix3d dv_dba = imupreint->getJVBiasa();
			Eigen::Matrix3d dv_dbg = imupreint->getJVBiasg();

            if (jacobians[0])
            {
                Eigen::Map<Eigen::Matrix<double, 15, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);
                jacobian_pose_i.setZero();

                jacobian_pose_i.block<3, 3>(O_P, O_P) = -Qi.inverse().toRotationMatrix();
                jacobian_pose_i.block<3, 3>(O_P, O_R) = Utility::skewSymmetric(Qi.inverse() * (0.5 * G * sum_dt * sum_dt + Pj - Pi - Vi * sum_dt));

#if 0
            jacobian_pose_i.block<3, 3>(O_R, O_R) = -(Qj.inverse() * Qi).toRotationMatrix();
#else
//                 Eigen::Quaterniond corrected_delta_q = delta_q * Utility::deltaQ(dq_dbg * (Bgi - linearized_bg));
				Eigen::Quaterniond corrected_delta_q = delta_q;
                jacobian_pose_i.block<3, 3>(O_R, O_R) = -(Utility::Qleft(Qj.inverse() * Qi) * Utility::Qright(corrected_delta_q)).bottomRightCorner<3, 3>();
#endif
                jacobian_pose_i.block<3, 3>(O_V, O_R) = Utility::skewSymmetric(Qi.inverse() * (G * sum_dt + Vj - Vi));

                jacobian_pose_i = sqrt_info * jacobian_pose_i;
				
#if NEGJACOB
				jacobian_pose_i = -jacobian_pose_i;
#endif

                if (jacobian_pose_i.maxCoeff() > 1e8 || jacobian_pose_i.minCoeff() < -1e8)
                {
					std::cout<<"unstable jacobian."<<std::endl;
                }
            }
            if (jacobians[1])
            {
                Eigen::Map<Eigen::Matrix<double, 15, 9, Eigen::RowMajor>> jacobian_speedbias_i(jacobians[1]);
                jacobian_speedbias_i.setZero();
                jacobian_speedbias_i.block<3, 3>(O_P, O_V - O_V) = -Qi.inverse().toRotationMatrix() * sum_dt;
                jacobian_speedbias_i.block<3, 3>(O_P, O_BA - O_V) = -dp_dba;
                jacobian_speedbias_i.block<3, 3>(O_P, O_BG - O_V) = -dp_dbg;

#if 0
            jacobian_speedbias_i.block<3, 3>(O_R, O_BG - O_V) = -dq_dbg;
#else
//                 Eigen::Quaterniond corrected_delta_q = delta_q * Utility::deltaQ(dq_dbg * (Bgi - linearized_bg));
				Eigen::Quaterniond corrected_delta_q = delta_q;
                jacobian_speedbias_i.block<3, 3>(O_R, O_BG - O_V) = -Utility::Qleft(Qj.inverse() * Qi * corrected_delta_q).bottomRightCorner<3, 3>() * dq_dbg;
#endif

                jacobian_speedbias_i.block<3, 3>(O_V, O_V - O_V) = -Qi.inverse().toRotationMatrix();
                jacobian_speedbias_i.block<3, 3>(O_V, O_BA - O_V) = -dv_dba;
                jacobian_speedbias_i.block<3, 3>(O_V, O_BG - O_V) = -dv_dbg;

                jacobian_speedbias_i.block<3, 3>(O_BA, O_BA - O_V) = -Eigen::Matrix3d::Identity();

                jacobian_speedbias_i.block<3, 3>(O_BG, O_BG - O_V) = -Eigen::Matrix3d::Identity();

                jacobian_speedbias_i = sqrt_info * jacobian_speedbias_i;
				
#if NEGJACOB
				jacobian_speedbias_i = -jacobian_speedbias_i;
#endif
            }
            if (jacobians[2])
            {
                Eigen::Map<Eigen::Matrix<double, 15, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[2]);
                jacobian_pose_j.setZero();

                jacobian_pose_j.block<3, 3>(O_P, O_P) = Qi.inverse().toRotationMatrix();

#if 0
            jacobian_pose_j.block<3, 3>(O_R, O_R) = Eigen::Matrix3d::Identity();
#else
//                 Eigen::Quaterniond corrected_delta_q = delta_q * Utility::deltaQ(dq_dbg * (Bgi - linearized_bg));
				Eigen::Quaterniond corrected_delta_q = delta_q;
                jacobian_pose_j.block<3, 3>(O_R, O_R) = Utility::Qleft(corrected_delta_q.inverse() * Qi.inverse() * Qj).bottomRightCorner<3, 3>();
#endif

                jacobian_pose_j = sqrt_info * jacobian_pose_j;
#if NEGJACOB
				jacobian_pose_j = -jacobian_pose_j;
#endif
            }
            if (jacobians[3])
            {
                Eigen::Map<Eigen::Matrix<double, 15, 9, Eigen::RowMajor>> jacobian_speedbias_j(jacobians[3]);
                jacobian_speedbias_j.setZero();

                jacobian_speedbias_j.block<3, 3>(O_V, O_V - O_V) = Qi.inverse().toRotationMatrix();

                jacobian_speedbias_j.block<3, 3>(O_BA, O_BA - O_V) = Eigen::Matrix3d::Identity();

                jacobian_speedbias_j.block<3, 3>(O_BG, O_BG - O_V) = Eigen::Matrix3d::Identity();

                jacobian_speedbias_j = sqrt_info * jacobian_speedbias_j;
#if NEGJACOB
				jacobian_speedbias_j = -jacobian_speedbias_j;
#endif

            }
        }

        return true;
    }

    // see (44, 45) from "on manifold pre_integration ..." for more info
    // obs - est ?
    Eigen::Matrix<double, 15, 1> Evaluate_Direct(const Eigen::Vector3d &Pi, const Eigen::Quaterniond &Qi, const Eigen::Vector3d &Vi, const 		  Eigen::Vector3d &Bai, const Eigen::Vector3d &Bgi, const Eigen::Vector3d &Pj, const Eigen::Quaterniond &Qj, const Eigen::Vector3d &Vj, const Eigen::Vector3d &Baj, const Eigen::Vector3d &Bgj) const
	{
		Eigen::Matrix<double, 15, 1> residuals_;
#if CRAZY
		Eigen:: Vector3d dba(0.0, 0.0, 0.0);
		Eigen:: Vector3d dbg(0.0, 0.0, 0.0);
#else
		Eigen:: Vector3d dba = Bai - linearized_ba;
		Eigen:: Vector3d dbg = Bgi - linearized_bg;
#endif
		Eigen::Vector3d delta_p = imupreint->getDeltaP();
		Eigen::Quaterniond delta_q = imupreint->getDeltaQ();
		Eigen:: Vector3d delta_v = imupreint->getDeltaV();
		
		Eigen::Matrix3d dp_dba = imupreint->getJPBiasa();
		Eigen::Matrix3d dp_dbg = imupreint->getJPBiasg();
		Eigen::Matrix3d dq_dbg = imupreint->getJRBiasg();
		Eigen::Matrix3d dv_dba = imupreint->getJVBiasa();
		Eigen::Matrix3d dv_dbg = imupreint->getJVBiasg();
		
		Eigen::Quaterniond corrected_delta_q = imupreint->getDeltaQ() * Utility::deltaQ(dq_dbg * dbg);
		Eigen::Vector3d corrected_delta_v = imupreint->getDeltaV() + dv_dba * dba + dv_dbg * dbg;
		Eigen::Vector3d corrected_delta_p = imupreint->getDeltaP() + dp_dba * dba + dp_dbg * dbg;
		
		residuals_.block<3,1>(O_P, 0) = Qi.inverse() * (0.5 * G * sum_dt * sum_dt + Pj - Pi - Vi * sum_dt) - corrected_delta_p;
		residuals_.block<3,1>(O_R, 0) = 2 * (corrected_delta_q.inverse() * (Qi.inverse() * Qj)).vec();// returns the imginary part of a quaternion. this approx w the angular velocity
		residuals_.block<3,1>(O_V, 0) = Qi.inverse() * (G * sum_dt + Vj - Vi) - corrected_delta_v;
#if CRAZY
		residuals_.block<3,1>(O_BA, 0).setZero(); // check this
		residuals_.block<3,1>(O_BG, 0).setZero();// check this
#else
		residuals_.block<3,1>(O_BA, 0) = Baj - Bai; // check this
		residuals_.block<3,1>(O_BG, 0) = Bgj - Bgi; // check this
#endif

		return residuals_;
	}

	const IMUPreintegrator* imupreint;
	double sum_dt;
    Eigen::Vector3d linearized_ba, linearized_bg; // form navstate.Get_biasAcc  Get_biasGyr()
    
	int O_P = 0;	int O_R = 3;	int O_V = 6;
	int O_BA = 9;	int O_BG = 12;

	Eigen::Matrix<double, 15, 15> covariance;
};

#endif
