#include "estimator.h"
// #include "utility/Converter.h"
#include <boost/concept_check.hpp>

#define LOG_EST

using namespace std;
using namespace Eigen;

Estimator::Estimator()
{
    std::cout<<"init begins"<<std::endl;
//     clearState();
// 	setParameter();
//     failure_occur = 0;
}

void Estimator::setParameter()
{
	// Gravity vector
// 	G << 1.00370258e-01, 9.30160713e+00, 3.11546779e+00;
// 	G <<  5.45640774e-02, 9.25659943e+00, 3.24784493e+00;
// 	G << -6.07326210e-01, 9.08026218e+00, 3.66279817e+00;
	G << -5.75561345e-01, 9.08304977e+00, 3.66101575e+00;
// 	G(0) = -6.07326210e-01; // G[1] = 9.08026218e+00; G[2] = 3.66279817e+00;
	G = -G;
// 	G[0] = 0.0; G[1] = 0.0; G[2] = 0.0; 
	std::cout<<"gravity " << G.transpose()<<std::endl;
	// cam imu extrinsic
	ric << 1.48655428e-02, -9.99880910e-01, 4.14029695e-03,
       9.99557257e-01, 1.49672134e-02, 2.57155299e-02, -2.57744361e-02,
       3.75618832e-03, 9.99660730e-01;
	tic << -2.16401462e-02, -6.46769851e-02, 9.81073081e-03;
	std::cout<<"tic in setParameter: "<<tic.transpose()<<std::endl;

	qic = ric; // weirdoooo.
// 	qic(ric);// rotation matrix to quaternion , check ric or rci
	// cam intrinsics
	fx = 4.35204681e+02;
	fy = 4.35204681e+02;
	cx = 3.67451721e+02;
	cy = 2.52200851e+02;
	bf = 4.7906394e+01;

}

void Estimator::clearState()
{
    for (int i = 0; i < WINDOW_SIZE + 1; i++)
    {
        Rs[i].setIdentity();
        Ps[i].setZero();
        Vs[i].setZero();
        Bas[i].setZero();
        Bgs[i].setZero();
        dt_buf[i].clear();
        linear_acceleration_buf[i].clear();
        angular_velocity_buf[i].clear();
    }

    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        tic = Eigen::Vector3d::Zero();
        ric = Eigen::Matrix3d::Identity();
		qic = ric;
    }

    solver_flag = INITIAL;
//     first_imu = false;

}

bool Estimator::readIMUPreint()
{
// 	std::string  imuFileName = "IMUFactor.yaml";
	std::string  imuFileName = "imupreint.yml";
// 	std::string  imuFileName = "imupreint1220.yml";
// 	std::ifstream imufile;
	cv::FileStorage imuFS(imuFileName, cv::FileStorage::READ);
	cv::FileNode preint_measures = imuFS["PreintData"];
	cv::FileNodeIterator it = preint_measures.begin();
	cv::FileNodeIterator it_end = preint_measures.end();
	
	// 
	std::cout<<"#FileNodes: "<<preint_measures.size()<<endl;

	// for current frame
	// reading state P Q V Ba Bg
	// reading imupreint delta measure and CovsPVPhisBag
	// reading jacob of delta measurement with respect to ba bg, that is
	// dp/dba dp/dbg dq/dba dq/dbg dv/dba dv/dbg
	int idx = 0;
	for(; it != it_end; ++it, idx++)
	{
		cv::Mat acc; cv::Mat gyr;

		cv::Mat Pi; cv::Mat Qi; cv::Mat Vi;
		cv::Mat Bai; cv::Mat Bgi;
		//dt_buf, linear_acceleration_buf, angular_velocity_buf
		// Pi, Qi, Vi, Bai, Bgi, dPi, dQi, dVi
		(*it)["Pi"] >> Pi;
// 		std::cout<<"Pi in read imupreint: "<<Pi<<std::endl;
		Ps.push_back(Converter::toVector3d(Pi));
		(*it)["Qi"] >> Qi;
		Eigen::Quaterniond eiQi = Converter::toEiQuaternion(Qi);
		Qs.push_back(eiQi);
// 		Rs.push_back(eiQi.toRotationMatrix());
		(*it)["Vi"] >> Vi;
		Vs.push_back(Converter::toVector3d(Vi));		
		(*it)["Bai"] >> Bai;
		Bas.push_back(Converter::toVector3d(Bai));
		(*it)["Bgi"] >> Bgi;
		Bgs.push_back(Converter::toVector3d(Bgi));
		
		double dT(0.0);
		cv::Mat dPi; /*cv::Mat dQi;*/ cv::Mat dRi; cv::Mat dVi; cv::Mat  dBai; cv::Mat  dBgi; 
		cv::Mat Cov;		
		(*it)["DeltaT"] >> dT;
		dTs.push_back(dT);
		(*it)["dPi"] >> dPi;
		dPs.push_back(Converter::toVector3d(dPi));
		(*it)["dRi"] >> dRi;
// 		Eigen::Matrix3d dRi_ = Converter::toMatrix3d(dRi);
// 		Eigen::Quaterniond dQi = dRi;
		Eigen::Quaterniond dQi(Converter::toMatrix3d(dRi));
		dQs.push_back(dQi);
// 		(*it)["dQi"] >> dQi;
// 		dQs.push_back(Converter::toEiQuaternion(dQi));
		(*it)["dVi"] >> dVi;
		dVs.push_back(Converter::toVector3d(dVi));
// 		(*it)["CovPVPhi"] >> Cov;
// 		CovsPVPhis.push_back(Converter::toMatrix9d(Cov));
		(*it)["CovPVPhi"] >> Cov;
		CovsPVPhisBag.push_back(Converter::toMatrix15d(Cov));
		(*it)["dBai"] >> dBai;
		dBas.push_back(Converter::toVector3d(dBai));		
		(*it)["dBgi"] >> dBgi;
		dBgs.push_back(Converter::toVector3d(dBgi));
		
		cv::Mat dp_dba; cv::Mat dp_dbg;
		cv::Mat dq_dbg;
		cv::Mat dv_dba; cv::Mat dv_dbg;
		(*it)["dpdba"] >> dp_dba;
		dp_dbas.push_back(Converter::toMatrix3d(dp_dba));
		(*it)["dpdbg"] >> dp_dbg;
		dp_dbgs.push_back(Converter::toMatrix3d(dp_dbg));
		(*it)["dqdbg"] >> dq_dbg;
		dq_dbgs.push_back(Converter::toMatrix3d(dq_dbg));
		(*it)["dvdba"] >> dv_dba;
		dv_dbas.push_back(Converter::toMatrix3d(dv_dba));
		(*it)["dvdbg"] >> dv_dbg;
		dv_dbgs.push_back(Converter::toMatrix3d(dv_dbg));
	}	
	
	// check data order
// 	for (int i = 0; i<Ps.size(); i++)
// 	{
// 		std::cout<<"Pi: "<<Ps[i].transpose()<<std::endl;
// 	}
	
	return true;
}

bool readIMUPreint(char* imuFileName)
{
	return false;
}

bool Estimator::readFeatures()
{
// 	std::string  prjFileName = "StereoPrjFactor.yaml";
// 	std::string  prjFileName = "stereo_prj.yml";
// 	std::cout<<"error opening yaml."<<std::endl;
	std::string  prjFileName = "str_prj_testdata.yml";
// 	std::cout<<"succed?";
	cv::FileStorage prjFS(prjFileName, cv::FileStorage::READ);
	
	// read nav pose
	// read projection measurements
// 	cv::FileNode prj_measures = prjFS["ProjectData"];
	cv::FileNode prj_measures = prjFS["Projections"];
	cv::FileNodeIterator it = prj_measures.begin();
	cv::FileNodeIterator it_end = prj_measures.end();	
	
	int idx = 0;
// 	std::cout<<"gonna read filenode"<<std::endl;
	for (; it != it_end; ++it, idx++)
	{
		cv::Mat Pw; cv::Mat obs; cv::Mat info;
		int pid; int fid;
		(*it)["FrameId"] >> fid; 
		frame_ids.push_back(fid);
		(*it)["PtId"] >> pid;
		(*it)["Pw"] >> Pw;
// 		point_cloud.insert(pair<int, Eigen::Vector3d>(fid, Converter::toVector3d(Pw)));
		point_cloud.push_back(Converter::toVector3d(Pw));
		(*it)["Obs"] >> obs;
// 		features.insert(pair<int, Eigen::Vector3d>(fid, Converter::toVector3d(obs)));
		features.push_back(Converter::toVector3d(obs));
		(*it)["Info"] >> info;
// 		stereo_infos.insert(pair<int, Eigen::Matrix3d>(fid, Converter::toMatrix3d(info)));
		stereo_infos.push_back(Converter::toMatrix3d(info)); 
	}
	std::cout<<"Pw size: "<<point_cloud.size()<<std::endl;
// 	std::cout<<"#FileNodes: "<<idx<<std::endl;	
	return true;
}

bool readFeatures(char* prjFileName)
{
	return false;
}

void Estimator::vector2double()
{
    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
		// using map instead of this cubersome way.
// 		std::cout<<i<<"th element in vector2double."<<endl;		
        para_Pose[i][0] = Ps[i][0];
        para_Pose[i][1] = Ps[i][1];
        para_Pose[i][2] = Ps[i][2];
// 		std::cout<<"vector2double pose"<<i<<": ";
// 		std::cout<<para_Pose[i][0]<<", "<<para_Pose[i][1]<<", "<<para_Pose[i][2]<<std::endl;
//         Eigen::Quaterniond q{Rs[i]};
//         para_Pose[i][3] = q.x();
//         para_Pose[i][4] = q.y();
//         para_Pose[i][5] = q.z();
//         para_Pose[i][6] = q.w();        
        para_Pose[i][3] = Qs[i].x();
        para_Pose[i][4] = Qs[i].y();
        para_Pose[i][5] = Qs[i].z();
        para_Pose[i][6] = Qs[i].w();

        para_SpeedBias[i][0] = Vs[i][0];
        para_SpeedBias[i][1] = Vs[i][1];
        para_SpeedBias[i][2] = Vs[i][2];

        para_SpeedBias[i][3] = Bas[i][0];
        para_SpeedBias[i][4] = Bas[i][1];
        para_SpeedBias[i][5] = Bas[i][2];

        para_SpeedBias[i][6] = Bgs[i][0];
        para_SpeedBias[i][7] = Bgs[i][1];
        para_SpeedBias[i][8] = Bgs[i][2];
#ifdef LOG_EST
		std::cout<<"init para_pose "<<i<<": ";
		for (int j = 0; j<7;j++)
		{
			std::cout<<para_Pose[i][j]<<", ";
		}
		std::cout<<std::endl;
		
// 		std::cout<<"para_SpeedBias: ";
// 		for (int j = 0; j<9;j++)
// 		{
// 			std::cout<<para_SpeedBias[i][j]<<", ";
// 		}
// 		std::cout<<std::endl;
#endif
    }

}

void Estimator::double2vector()
{
    Eigen::Vector3d origin_R0 = Utility::R2ypr(Rs[0]);
    Eigen::Vector3d origin_P0 = Ps[0];

    Eigen::Vector3d origin_R00 = Utility::R2ypr(Eigen::Quaterniond(para_Pose[0][6],
                                                      para_Pose[0][3],
                                                      para_Pose[0][4],
                                                      para_Pose[0][5]).toRotationMatrix());
    double y_diff = origin_R0(0) - origin_R00(0);
    //TODO
    Matrix3d rot_diff = Utility::ypr2R(Vector3d(y_diff, 0, 0));

    for (int i = 0; i <= WINDOW_SIZE; i++)
    {

        Rs[i] = rot_diff * Quaterniond(para_Pose[i][6], para_Pose[i][3], para_Pose[i][4], para_Pose[i][5]).normalized().toRotationMatrix();
        Ps[i] = rot_diff * Vector3d(para_Pose[i][0] - para_Pose[0][0],
                                para_Pose[i][1] - para_Pose[0][1],
                                para_Pose[i][2] - para_Pose[0][2]) + origin_P0;
        Vs[i] = rot_diff * Vector3d(para_SpeedBias[i][0],
                                    para_SpeedBias[i][1],
                                    para_SpeedBias[i][2]);

        Bas[i] = Vector3d(para_SpeedBias[i][3],
                          para_SpeedBias[i][4],
                          para_SpeedBias[i][5]);

        Bgs[i] = Vector3d(para_SpeedBias[i][6],
                          para_SpeedBias[i][7],
                          para_SpeedBias[i][8]);
    }
	// no online cam imu calibration for now.
}

void Estimator::optimization_simple()
{
    ceres::Problem problem;
    ceres::LossFunction *loss_function;
    //loss_function = new ceres::HuberLoss(1.0);
    loss_function = new ceres::CauchyLoss(1.0);
	// fixed the oldest key frame in sliding window.
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_Pose[i], SIZE_POSE, local_parameterization);
        problem.AddParameterBlock(para_SpeedBias[i], SIZE_SPEEDBIAS);
    }

    TicToc t_whole, t_prepare;
    // initial estimation.
// 	std::cout<<"para_Pose size "<<Ps.size()<<std::endl;
	std::cout<<"vector2double Ps size "<<Ps.size()<<std::endl;
	// initial estimation
	vector2double();
	// adding marginalization factor, no marginalization for now.
    // adding imu factor
	// i = 0 for current frame, i = window_size - 1 for the last key frame
#if 1
	std::cout<<"imu factor construction."<<std::endl;
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
		// params for imu factor constructor
		double tsum_dt = dTs[i];
		Eigen::Vector3d tdeltaP = dPs[i];
		Eigen::Quaterniond tdeltaQ = dQs[i];
		std::cout<<"Measure Delta Q: "<< tdeltaQ.w()<<", "<<tdeltaQ.vec().transpose() << std::endl;
		Eigen::Vector3d tdeltaV = dVs[i];
// 		Eigen::Matrix<double, 9,9> tCov = CovsPVPhis[i];
		Eigen::Matrix<double, 15,15> tCov = CovsPVPhisBag[i];

		Eigen::Vector3d tlba = dBas[i]; // linearized_bas
		Eigen::Vector3d tlbg = dBgs[i]; // linearized_bgs
		Eigen::Matrix3d tdp_dba = dp_dbas[i];
		Eigen::Matrix3d tdp_dbg = dp_dbgs[i];
		Eigen::Matrix3d tdq_dbg = dq_dbgs[i];
		Eigen::Matrix3d tdv_dba = dv_dbas[i];
		Eigen::Matrix3d tdv_dbg = dv_dbgs[i];
				
		// using these params to construct imu factor
//         int j = i + 1;
        IMUFactor* imu_factor = new IMUFactor(tsum_dt, tdeltaP, tdeltaQ, tdeltaV, tCov, tlba, tlbg, tdp_dba, tdp_dbg, tdq_dbg, tdv_dba, tdv_dbg);
        if (imu_factor->sum_dt > 10.0)
            continue;
        problem.AddResidualBlock(imu_factor, NULL, para_Pose[i+1], para_SpeedBias[i+1], para_Pose[i], para_SpeedBias[i]);
    }    
#endif
	// adding current frame residual block.
	// adding slide window residual block.
	// frameid = 0 for current frame, the largest frame id for the oldest keyframe
	std::cout<<"projection factor construction."<<std::endl;
// 	std::map<int, int> frameIdtoPoseId {{498,0},{101,1},{100,2},{99,3},{98,4}, {97,5}};
// 	std::map<int, int> frameIdtoPoseId {{497,0},{101,1},{100,2},{99,3},{98,4}, {97,5}};
	std::map<int, int> frameIdtoPoseId {{419,0},{98,1},{97,2},{96,3},{95,4}, {94,5}};
	for (int i = 0; i < point_cloud.size(); i++)
	{
// make sure frame id and para_Pose index match
// initialize StereoPrjFactor using point and its corresponding feature
		int fid = frame_ids[i];
		Eigen::Vector3d pw = point_cloud[i];
		Eigen::Vector3d ft = features[i];
		Eigen::Matrix3d info = stereo_infos[i];
		int idx = frameIdtoPoseId[fid];
		ProjectionFactor* prj_factor = new ProjectionFactor(pw, ft, info);
// 		double **ppara_pose = para_Pose[idx];
		double *ppara_pose = &para_Pose[idx][0]; 
// 		prj_factor->check(&(ppara_pose));
		// I deleted para_Feature for now
		problem.AddResidualBlock(prj_factor, NULL, para_Pose[idx]);// could do ba here. I mean cam extrin calib, point refining.
	}

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
	std::cout<<"gonna solve."<<endl;
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << endl;
//     ROS_DEBUG("Iterations : %d", static_cast<int>(summary.iterations.size()));
//     ROS_DEBUG("solver costs: %f", t_solver.toc());

    // relative info between two loop frame, removed by ted.
//     double2vector();

	// marginalization! remove this for now!

}
