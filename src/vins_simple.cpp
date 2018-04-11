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

Estimator estimator;

int global_frame_cnt = 0;
using namespace std;

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
		
	estimator.optimization_simple();
	std::cout<<"optimization done."<<std::endl;

	cout<<"optimization results:"<<endl;
	cout<<"P: x y z; Q: x y z w"<<endl;
	for (int i =0; i<WINDOW_SIZE; i++)
	{
		cout<<"pose"<<i<<": ";
		for (int j = 0; j < 7; j++)
		{
			cout<<estimator.para_Pose[i][j]<<", ";
		}
		cout<<endl;
	}
	// prepare imu factor data and projection factor data
//     std::thread measurement_process{process};
// 	estimator.optimization_simple();

    return 0;
}
