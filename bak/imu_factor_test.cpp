#include <stdio.h>
#include <queue>
#include <map>

#include <opencv2/opencv.hpp>

#include "estimator.h"
#include "parameters.h"

// Estimator estimator;

double latest_time;
Eigen::Vector3d tmp_P;
Eigen::Quaterniond tmp_Q;
Eigen::Vector3d tmp_V;
Eigen::Vector3d tmp_Ba;
Eigen::Vector3d tmp_Bg;
Eigen::Vector3d acc_0;
Eigen::Vector3d gyr_0;

using namespace std;
Estimator estimator;

int global_frame_cnt = 0;

int main(int argc, char **argv)
{
	// check your input.
	std::cout<<"let's start."<<std::endl;
	estimator.setParameter();
	
	// read parameters.
	estimator.readIMUPreint();
	std::cout<<"readIMUPreint done."<<std::endl;
	estimator.readFeatures();
	std::cout<<"read features done."<<std::endl;
	int WIN_SIZE = estimator.Ps.size();
// 		WIN_SIZE = 2;
	// initialize para pose
	std::cout<< "window size: " << WIN_SIZE <<std::endl;
	double para_pose[6][7];
	double para_SpeedBias[6][9];
	for (int i = 0; i < WIN_SIZE; i++)
	{
		Eigen::Vector3d Pi = estimator.Ps[i];
		Eigen::Quaterniond Qi = estimator.Qs[i];
		std::cout<<"pose "<< i << " :"<<Pi.transpose();
		std::cout<<Qi.z()<<", "<<Qi.vec()<<endl;
		Eigen::Vector3d Vi = estimator.Vs[i];
		Eigen::Vector3d Bai = estimator.Bas[i];
		Eigen::Vector3d Bgi = estimator.Bgs[i];
		para_pose[i][0] = Pi[0];
		para_pose[i][1] = Pi[1];
		para_pose[i][2] = Pi[2];
		
		para_pose[i][3] = Qi.x();
		para_pose[i][4] = Qi.y();
		para_pose[i][5] = Qi.z();
		para_pose[i][6] = Qi.w();
		
		para_SpeedBias[i][0] = Vi[0];
		para_SpeedBias[i][1] = Vi[1];
		para_SpeedBias[i][2] = Vi[2];
		
		para_SpeedBias[i][3] = Bai[0];
		para_SpeedBias[i][4] = Bai[1];
		para_SpeedBias[i][5] = Bai[2];
		
		para_SpeedBias[i][6] = Bgi[0];
		para_SpeedBias[i][7] = Bgi[1];
		para_SpeedBias[i][8] = Bgi[2];
	}
	
	for (int i = 0; i<WIN_SIZE; i++)
	{
		cout<<"para_pose init "<< i<<": ";
		for (int j = 0; j<7;j++)
			cout<<para_pose[i][j]<<" ";
		cout<<endl;
	}
	
	const int SIZE_POSE(7);
	const int WINDOW_SIZE(2);
	// 0 current frame, 1 last key frame
	// x,y,z, x, y, z, w
// 	double para_pose[2][SIZE_POSE] = \
// 	{{1.09683394e-01, 3.62280667e-01, 2.86611021e-01, \
// 		1.15474097e-01, -7.28111416e-02,
//              -7.04158127e-01, 6.96796715e-01}, 
// 	{1.09700903e-01, 3.62365097e-01, 2.86691636e-01 , \
// 		1.15474097e-01, -7.28111416e-02,
//              -7.04158127e-01, 6.96796715e-01}};
// 	double para_SpeedBias[2][9] = \
// 	{{4.43052442e-04, 2.23113280e-02, 6.68430654e-03, \
// 		-7.58147165e-02, 2.02566341e-01, -1.52722234e-03, \
// 		-3.15532088e-03, 2.11417694e-02, 7.87554458e-02}, \
// 	 {-2.59988033e-03, 1.80969648e-02, 3.04262550e-03, \
// 	  -7.58147165e-02, 2.02566341e-01, -1.52722234e-03, \
// 	  -3.15532088e-03, 2.11417694e-02, 7.87554458e-02 }};	
	
    ceres::Problem problem;
    ceres::LossFunction *loss_function;
    //loss_function = new ceres::HuberLoss(1.0);
    loss_function = new ceres::CauchyLoss(1.0);
	// fixed the oldest key frame in sliding window.
    for (int i = 0; i < WIN_SIZE - 1; i++)
    {
        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_pose[i], SIZE_POSE, local_parameterization);
    }

    TicToc t_whole, t_prepare;
    // initial estimation.
	// para_Pose
	// constructing imu factor and adding residuals block.
	// adding residuals for imu preint measurements.
// window keyframe imupreint measurements.
	int i = 0;
	for (; i < WIN_SIZE - 1; i++)
	{
		double tsum_dt = estimator.dTs[i];
		Eigen::Vector3d tdeltaP = estimator.dPs[i];		
		Eigen::Quaterniond tdeltaQ = estimator.dQs[i];
// 		std::cout<<"Delta Q: "<< tdeltaQ << std::endl;
		Eigen::Vector3d tdeltaV = estimator.dVs[i];
		Eigen::Matrix<double, 9,9> tCov = estimator.CovsPVPhis[i];

		Eigen::Vector3d tlba(-7.58147165e-02, 2.02566341e-01, -1.52722234e-03) ; // linearized_bas
		Eigen::Vector3d tlbg(-3.15532088e-03, 2.11417694e-02, 7.87554458e-02) ; // linearized_bgs
		
		Eigen::Matrix3d tdp_dba = estimator.dp_dbas[i];
		Eigen::Matrix3d tdp_dbg = estimator.dp_dbgs[i];
		Eigen::Matrix3d tdq_dbg = estimator.dq_dbgs[i];
		Eigen::Matrix3d tdv_dba = estimator.dv_dbas[i];
		Eigen::Matrix3d tdv_dbg = estimator.dv_dbgs[i];
#if 0
		cout<<"sum_dt" << tsum_dt <<endl;
		cout<<"deltaP "<< tdeltaP.transpose() <<endl;
		cout<<"deltaQ " << "w: "<<tdeltaQ.w() << ", x: "<< tdeltaQ.x();
		cout<< ", y: "<< tdeltaQ.y() << ", z: " <<tdeltaQ.z() <<endl;
		cout<<"deltaV " << tdeltaV.transpose() << endl;
		cout << "Cov9x9 " << tCov << endl;
		cout << "lba " << tlba.transpose() << endl;
		cout << "lbg " << tlbg.transpose() << endl;
		cout << "now for jacobians"<<endl<<endl;
		cout << "dpdba " << tdp_dba <<endl;
		cout << "dpdbg " << tdp_dbg <<endl;
		cout << "dqdbg " << tdq_dbg << endl;
		cout << "dvdba" << tdv_dba << endl;
		cout << "dvdbg" << tdv_dbg << endl;
#endif
        IMUFactor* imu_factor = new IMUFactor(tsum_dt, tdeltaP, tdeltaQ, tdeltaV, tCov, tlba, tlbg, tdp_dba, tdp_dbg, tdq_dbg, tdv_dba, tdv_dbg);
		
        if (imu_factor->sum_dt > 10.0)
		{
			std::cout<<"check your preint data, delta time too large"<<std::endl;
		}
        problem.AddResidualBlock(imu_factor, NULL, para_pose[i+1], para_SpeedBias[i+1], para_pose[i], para_SpeedBias[i]);				
	}

// for curent frame and last keyframe.	
// 	std::cout<< "constructing imu factor."<<std::endl;
// 	{
// 		// params for imu factor constructor
// 		double tsum_dt = 4.9999952316284180e-02;
// 		Eigen::Vector3d tdeltaP(1.30694984e-02, -9.76127107e-04, -3.55784921e-03);		
// 		Eigen::Matrix3d deltaR;
// 		deltaR << 1., -3.28537571e-05, -2.90895514e-05, 3.28543465e-05, 1.,
//           2.02602605e-05, 2.90888856e-05, -2.02612173e-05, 1.;
// 		Eigen::Quaterniond tdeltaQ(deltaR);
// // 		std::cout<<"Delta Q: "<< tdeltaQ << std::endl;
// 		Eigen::Vector3d tdeltaV(4.74552274e-01, -3.53821777e-02, -1.29062742e-01);
// 		Eigen::Matrix<double, 9,9> tCov;
// 		tCov << 1.25726046e-07, 3.84634300e-13, 1.42010210e-12,
//           3.78125196e-06, 2.02881080e-11, 7.49551116e-11,
//           -8.37472740e-15, -2.22350013e-10, 6.02108363e-11,
//           3.84634300e-13, 1.25731233e-07, -1.04635620e-13,
//           2.04352750e-11, 3.78152617e-06, -5.55754531e-12,
//           2.22371357e-10, -2.37357165e-14, 8.17918400e-10,
//           1.42010210e-12, -1.04635620e-13, 1.25730878e-07,
//           7.48710260e-11, -5.51135309e-12, 3.78150730e-06,
//           -6.02514288e-11, -8.17921786e-10, -1.20669245e-14,
//           3.78125196e-06, 2.04352750e-11, 7.48710260e-11, 1.51250948e-04,
//           1.15098986e-09, 4.21208579e-09, -5.26764288e-13,
//           -1.40188696e-08, 3.84196275e-09, 2.02881080e-11,
//           3.78152617e-06, -5.51135309e-12, 1.15098986e-09,
//           1.51266358e-04, -3.12591869e-10, 1.40202179e-08,
//           -1.49760456e-12, 5.16342382e-08, 7.49551116e-11,
//           -5.55754531e-12, 3.78150730e-06, 4.21208579e-09,
//           -3.12591869e-10, 1.51265311e-04, -3.84452514e-09,
//           -5.16344478e-08, -7.60651904e-13, -8.37472740e-15,
//           2.22371357e-10, -6.02514288e-11, -5.26764288e-13,
//           1.40202179e-08, -3.84452514e-09, 2.41999544e-07,
//           1.54646343e-20, -5.44241571e-19, -2.22350013e-10,
//           -2.37357165e-14, -8.17921786e-10, -1.40188696e-08,
//           -1.49760456e-12, -5.16344478e-08, 1.54646359e-20,
//           2.41999544e-07, 3.19610672e-19, 6.02108363e-11, 8.17918400e-10,
//           -1.20669245e-14, 3.84196275e-09, 5.16342382e-08,
//           -7.60651904e-13, -5.44241571e-19, 3.19610672e-19,
//           2.41999544e-07 ; 
// 
// 		Eigen::Vector3d tlba(-7.58147165e-02, 2.02566341e-01, -1.52722234e-03) ; // linearized_bas
// 		Eigen::Vector3d tlbg(-3.15532088e-03, 2.11417694e-02, 7.87554458e-02) ; // linearized_bgs
// 		Eigen::Matrix3d tdp_dba ;
// 		tdp_dba << -1.37499743e-03, 5.83257425e-08, -3.42100126e-09,
//           -5.83257425e-08, -1.37499743e-03, -1.44103982e-08,
//           3.42183837e-09, 1.44101016e-08, -1.37499743e-03;
// 		Eigen::Matrix3d tdp_dbg;
// 		tdp_dbg << -3.43534407e-06, 5.37242231e-05, -1.45732292e-05,
//           -5.37234482e-05, -3.43436318e-06, -1.97665344e-04,
//           1.45812546e-05, 1.97664587e-04, -3.43641818e-06;
// 		Eigen::Matrix3d tdq_dbg;
// 		tdq_dbg << -4.99999523e-02, 4.84702241e-07, -1.35268078e-06,
//           -4.84710938e-07, -4.99999523e-02, 1.86567988e-07,
//           1.35268431e-06, -1.86587670e-07, -4.99999523e-02;
// 		Eigen::Matrix3d tdv_dba;
// 		tdv_dba << -4.99999523e-02, 2.12739496e-06, 1.01791080e-07,
//           -2.12740042e-06, -4.99999523e-02, -8.26487337e-07,
//           -1.01750437e-07, 8.26485120e-07, -4.99999523e-02;
// 		Eigen::Matrix3d tdv_dbg;
// 		tdv_dbg << 1.39378798e-07, 3.21797910e-03, -8.81392218e-04,
//           -3.21793696e-03, 2.55085780e-07, -1.18444739e-02,
//           8.81928776e-04, 1.18444273e-02, 1.20437676e-07 ;
// 
// 		cout<<"sum_dt" << tsum_dt <<endl;
// 		cout<<"deltaP "<< tdeltaP.transpose() <<endl;
// 		cout<<"deltaQ " << "w: "<<tdeltaQ.w() << ", x: "<< tdeltaQ.x();
// 		cout<< ", y: "<< tdeltaQ.y() << ", z: " <<tdeltaQ.z() <<endl;
// 		cout<<"deltaV " << tdeltaV.transpose() << endl;
// 		cout << "Cov9x9 " << tCov << endl;
// 		cout << "lba " << tlba.transpose() << endl;
// 		cout << "lbg " << tlbg.transpose() << endl;
// 		cout << "now for jacobians"<<endl<<endl;
// 		cout << "dpdba " << tdp_dba <<endl;
// 		cout << "dpdbg " << tdp_dbg <<endl;
// 		cout << "dqdbg " << tdq_dbg << endl;
// 		cout << "dvdba" << tdv_dba << endl;
// 		cout << "dvdbg" << tdv_dbg << endl;
//         IMUFactor* imu_factor = new IMUFactor(tsum_dt, tdeltaP, tdeltaQ, tdeltaV, tCov, tlba, tlbg, tdp_dba, tdp_dbg, tdq_dbg, tdv_dba, tdv_dbg);
// 		int i = 0; int j = 1; // the order?
//         if (imu_factor->sum_dt > 10.0)
// 		{
// 			std::cout<<"check your preint data, delta time too large"<<std::endl;
// 		}
//         problem.AddResidualBlock(imu_factor, NULL, para_Pose[i], para_SpeedBias[i], para_Pose[j], para_SpeedBias[j]);		
// 	}

    ceres::Solver::Options options;

    options.linear_solver_type = ceres::DENSE_SCHUR;
    //options.num_threads = 2;
    options.trust_region_strategy_type = ceres::DOGLEG;
    options.max_num_iterations = NUM_ITERATIONS;
    //options.use_explicit_schur_complement = true;
    //options.minimizer_progress_to_stdout = true;
    //options.use_nonmonotonic_steps = true;
    options.max_solver_time_in_seconds = SOLVER_TIME;
    TicToc t_solver;
    ceres::Solver::Summary summary;
	std::cout<<"gonna solve."<< std::endl;
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << std::endl;
//     ROS_DEBUG("Iterations : %d", static_cast<int>(summary.iterations.size()));
//     ROS_DEBUG("solver costs: %f", t_solver.toc());

    // relative info between two loop frame, removed by ted.
//     double2vector();

	// marginalization! remove this for now!
	std::cout<<"optimization done."<<std::endl;

	// prepare imu factor data and projection factor data
//     std::thread measurement_process{process};
// 	estimator.optimization_simple();
	for (int i = 0; i<WIN_SIZE; i++)
	{
		cout<<"para_pose "<< i<<": ";
		for (int j = 0; j<7;j++)
			cout<<para_pose[i][j]<<" ";
		cout<<endl;
	}

    return 0;
}
