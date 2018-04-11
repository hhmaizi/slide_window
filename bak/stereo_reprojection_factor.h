#ifndef CERES_EXAMPLES_STEREO_REPROJECTION_ERROR_H_
#define CERES_EXAMPLES_STEREO_REPROJECTION_FACTOR_H_

#include "ceres/rotation.h"
#include <ceres/ceres.h>
#include <Eigen/Dense>
#include "../utility/utility.h"
#include "../utility/tic_toc.h"
#include "../parameters.h"

namespace ceres {
namespace examples {

// Templated pinhole camera model for used with Ceres.  The camera is
// parameterized using 9 parameters: 3 for rotation, 3 for translation, 1 for
// focal length and 2 for radial distortion. The principal point is not modeled
// (i.e. it is assumed be located at the image center).
struct SnavelyReprojectionError {
  SnavelyReprojectionError(double observed_ul, double observed_vl, double observed_ur)
      : observed_ul(observed_ul), observed_vl(observed_vl), observed_ur(observed_ur) {}

  template <typename T>
  bool operator()(const T* const camera,
                  const T* const point,
                  T* residuals) const {
    // camera[0,1,2] are the angle-axis rotation.
    T p[3];
    AngleAxisRotatePoint(camera, point, p);

    // camera[3,4,5] are the translation.
    p[0] += camera[3];
    p[1] += camera[4];
    p[2] += camera[5];

    // Compute the center of distortion. The sign change comes from
    // the camera model that Noah Snavely's Bundler assumes, whereby
    // the camera coordinate system has a negative z axis.
    const T xp = - p[0] / p[2];
    const T yp = - p[1] / p[2];
	const T bf(30.0);
 
    // Compute final projected point position.
    const T& focal = camera[6];
    const T predicted_x = focal * xp;
    const T predicted_y = focal * yp;
	const T predicted_xr = predicted_x - bf/p[2];

    // The error is the difference between the predicted and observed position.
    residuals[0] = predicted_x - observed_ul;
    residuals[1] = predicted_y - observed_vl;
	residuals[2] = predicted_xr - observed_ur;

    return true;
  }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create(const double observed_ul,
                                     const double observed_vl,
									 const double observed_ur) {
    return (new ceres::AutoDiffCostFunction<SnavelyReprojectionError, 3, 6, 3>(
                new SnavelyReprojectionError(observed_ul, observed_vl, observed_ur)));
  }

  double observed_ul;
  double observed_vl;
  double observed_ur;
};

// Templated pinhole camera model for used with Ceres.  The camera is
// parameterized using 10 parameters. 4 for rotation, 3 for
// translation, 1 for focal length and 2 for radial distortion. The
// principal point is not modeled (i.e. it is assumed be located at
// the image center).
struct SnavelyReprojectionErrorWithQuaternions {
  // (u, v): the position of the observation with respect to the image
  // center point.
  SnavelyReprojectionErrorWithQuaternions(double observed_ul, double observed_vl, double observed_ur)
      : observed_ul(observed_ul), observed_vl(observed_vl), observed_ur(observed_ur) {}

  template <typename T>
  bool operator()(const T* const camera,
                  const T* const point,
                  T* residuals) const {
    // camera[0,1,2,3] is are the rotation of the camera as a quaternion.
    //
    // We use QuaternionRotatePoint as it does not assume that the
    // quaternion is normalized, since one of the ways to run the
    // bundle adjuster is to let Ceres optimize all 4 quaternion
    // parameters without a local parameterization.
    T p[3];
    QuaternionRotatePoint(camera, point, p);

    p[0] += camera[4];
    p[1] += camera[5];
    p[2] += camera[6];
	
	// Compute the center of distortion. The sign change comes from
    // the camera model that Noah Snavely's Bundler assumes, whereby
    // the camera coordinate system has a negative z axis.
    const T xp = - p[0] / p[2];
    const T yp = - p[1] / p[2];
	const T bf(30.0);
 
    // Compute final projected point position.
    const T& focal = camera[7];
    const T predicted_x = focal * xp;
    const T predicted_y = focal * yp;
	const T predicted_xr = predicted_x - bf/p[2];

    // The error is the difference between the predicted and observed position.
    residuals[0] = predicted_x - observed_ul;
    residuals[1] = predicted_y - observed_vl;
	residuals[2] = predicted_xr - observed_ur;

    return true;
  }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create(const double observed_ul,
                                     const double observed_vl) {
    return (new ceres::AutoDiffCostFunction<
            SnavelyReprojectionErrorWithQuaternions, 3, 10, 3>(
                new SnavelyReprojectionErrorWithQuaternions(observed_ul,
                                                            observed_vl)));
  }

  double observed_ul;
  double observed_vl;
  double observed_ur;
};

}  // namespace examples
}  // namespace ceres

#endif  // CERES_EXAMPLES_SNAVELY_REPROJECTION_ERROR_H_
