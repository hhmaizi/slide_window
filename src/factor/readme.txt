## stereo pose only projection jacobian
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

