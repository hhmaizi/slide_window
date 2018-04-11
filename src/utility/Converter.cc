/**
* This file is part of ORB-SLAM2.
*
* Copyright (C) 2014-2016 Raúl Mur-Artal <raulmur at unizar dot es> (University of Zaragoza)
* For more information see <https://github.com/raulmur/ORB_SLAM2>
*
* ORB-SLAM2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with ORB-SLAM2. If not, see <http://www.gnu.org/licenses/>.
*/

#include <iomanip>
#include "Converter.h"
using namespace std;

// cv::Mat Converter::toCvMatInverse(const cv::Mat &Tcw)
// {
//     cv::Mat Rcw = Tcw.rowRange(0,3).colRange(0,3);
//     cv::Mat tcw = Tcw.rowRange(0,3).col(3);
//     cv::Mat Rwc = Rcw.t();
//     cv::Mat twc = -Rwc*tcw;
// 
//     cv::Mat Twc = cv::Mat::eye(4,4,Tcw.type());
//     Rwc.copyTo(Twc.rowRange(0,3).colRange(0,3));
//     twc.copyTo(Twc.rowRange(0,3).col(3));
// 
//     return Twc.clone();
// }

cv::Mat Converter::toCvMat(const Eigen::Matrix<double,4,4> &m)
{
    cv::Mat cvMat(4,4,CV_32F);
    for(int i=0;i<4;i++)
        for(int j=0; j<4; j++)
            cvMat.at<float>(i,j)=m(i,j);

    return cvMat.clone();
}

cv::Mat Converter::toCvMat(const Eigen::Matrix3d &m)
{
    cv::Mat cvMat(3,3,CV_32F);
    for(int i=0;i<3;i++)
        for(int j=0; j<3; j++)
            cvMat.at<float>(i,j)=m(i,j);

    return cvMat.clone();
}

cv::Mat Converter::toCvMat(const Eigen::Matrix<double,3,1> &m)
{
    cv::Mat cvMat(3,1,CV_32F);
    for(int i=0;i<3;i++)
            cvMat.at<float>(i)=m(i);

    return cvMat.clone();
}

cv::Mat Converter::toCvSE3(const Eigen::Matrix<double,3,3> &R, const Eigen::Matrix<double,3,1> &t)
{
    cv::Mat cvMat = cv::Mat::eye(4,4,CV_32F);
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            cvMat.at<float>(i,j)=R(i,j);
        }
    }
    for(int i=0;i<3;i++)
    {
        cvMat.at<float>(i,3)=t(i);
    }

    return cvMat.clone();
}

Eigen::Matrix<double,3,1> Converter::toVector3d(const cv::Mat &cvVector)
{
    Eigen::Matrix<double,3,1> v;
    v << cvVector.at<float>(0), cvVector.at<float>(1), cvVector.at<float>(2);

    return v;
}

Eigen::Matrix<double,3,1> Converter::toVector3d(const cv::Point3f &cvPoint)
{
    Eigen::Matrix<double,3,1> v;
    v << cvPoint.x, cvPoint.y, cvPoint.z;

    return v;
}

Eigen::Vector4d Converter::toVector4d(const cv::Mat& M)
{
	Eigen::Vector4d v;
	v << M.at<float>(0,0), M.at<float>(0,1), M.at<float>(0,2), M.at<float>(0,3);
	return v;
}

Eigen::Quaterniond Converter::toEiQuaternion(const cv::Mat& M)
{
	Eigen::Quaterniond Q(M.at<float>(0,0), M.at<float>(0,1), M.at<float>(0,2), M.at<float>(0,3));
	Q.normalize();
	return Q;

}


Eigen::Matrix<double,3,3> Converter::toMatrix3d(const cv::Mat &cvMat3)
{
    Eigen::Matrix<double,3,3> M;

    M << cvMat3.at<float>(0,0), cvMat3.at<float>(0,1), cvMat3.at<float>(0,2),
         cvMat3.at<float>(1,0), cvMat3.at<float>(1,1), cvMat3.at<float>(1,2),
         cvMat3.at<float>(2,0), cvMat3.at<float>(2,1), cvMat3.at<float>(2,2);

    return M;
}

Eigen::Matrix<double,9,9> Converter::toMatrix9d(const cv::Mat &M9x9)
{
// 	Eigen::Matrix<double,9,9> M;
	double m[81];
	for(int i = 0; i<9; i++)
	{
		for (int j=0; j<9; j++)
		{
			m[i*9 +j] = M9x9.at<float>(i,j);
		}		
	}
	
	Eigen::Map<Eigen::Matrix<double,9,9, Eigen::RowMajor> > M(m);
// 	Eigen::Matrix<double,9,9> M = Eigen::Map<Eigen::Matrix<double,9,9, Eigen::RowMajor>>(m);
	
	return M;// clone or just a map? what happens if m was changed?
}

Eigen::Matrix<double, 15, 15> Converter::toMatrix15d(const cv::Mat & M15x15)
{
	double m[225];
	for(int i = 0; i<15; i++)
	{
		for (int j=0; j<15; j++)
		{
			m[i*15 +j] = M15x15.at<float>(i,j);
		}		
	}
	
	Eigen::Map<Eigen::Matrix<double,15,15, Eigen::RowMajor> > M(m);
// 	Eigen::Matrix<double,9,9> M = Eigen::Map<Eigen::Matrix<double,9,9, Eigen::RowMajor>>(m);
	
	return M;// clone or just a map? what happens if m was changed?

}

std::vector<float> Converter::toQuaternion(const cv::Mat &M)
{
    Eigen::Matrix<double,3,3> eigMat = toMatrix3d(M);
    Eigen::Quaterniond q(eigMat);

    std::vector<float> v(4);
    v[0] = q.x();
    v[1] = q.y();
    v[2] = q.z();
    v[3] = q.w();

    return v;
}
