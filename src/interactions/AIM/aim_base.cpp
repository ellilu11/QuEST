#include "aim_base.h"

const Eigen::Vector3cd AIM::AimBase::FDTD_Del_Del( 
    // const std::vector<Eigen::Vector3cd> stencil 
    boost::multi_array<Eigen::Vector3cd, 1> stencil
    )
    
{
  constexpr int O = 0;
  
  constexpr int XP = 9;
  constexpr int XM = 18; 
  constexpr int YP = 3;
  constexpr int YM = 6;
  constexpr int ZP = 1;
  constexpr int ZM = 2;

  constexpr int XYPP = 12;
  constexpr int XYMP = 21;
  constexpr int XYPM = 15;
  constexpr int XYMM = 24;
  
  constexpr int XZPP = 10;
  constexpr int XZMP = 19;
  constexpr int XZPM = 11;
  constexpr int XZMM = 20;

  constexpr int YZPP = 4;
  constexpr int YZMP = 7;
  constexpr int YZPM = 5;
  constexpr int YZMM = 8;

  Eigen::Vector3cd D_XX, D_YY, D_ZZ, D_XY, D_XZ, D_YZ;
  D_XX = ( stencil[XP] - 2.0*stencil[O] + stencil[XM] ) / pow(h_, 2);
  D_YY = ( stencil[YP] - 2.0*stencil[O] + stencil[YM] ) / pow(h_, 2);
  D_ZZ = ( stencil[ZP] - 2.0*stencil[O] + stencil[ZM] ) / pow(h_, 2);

  D_XY = ( stencil[XYPP] - stencil[XYMP] - stencil[XYPM] + stencil[XYMM] ) / ( 4.0*pow(h_, 2) );
  D_XZ = ( stencil[XZPP] - stencil[XZMP] - stencil[XZPM] + stencil[XZMM] ) / ( 4.0*pow(h_, 2) );
  D_YZ = ( stencil[YZPP] - stencil[YZMP] - stencil[YZPM] + stencil[YZMM] ) / ( 4.0*pow(h_, 2) );

  cmplx E_X = ( D_XX[0] + D_XY[1] + D_XZ[2] );
  cmplx E_Y = ( D_XY[0] + D_YY[1] + D_YZ[2] );
  cmplx E_Z = ( D_XZ[0] + D_YZ[1] + D_ZZ[2] );

  Eigen::Vector3cd field(E_X, E_Y, E_Z);

  return field;
}
