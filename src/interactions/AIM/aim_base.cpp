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
  D_XX = stencil[XP] - 2.0*stencil[O] + stencil[XM];
  D_YY = stencil[YP] - 2.0*stencil[O] + stencil[YM];
  D_ZZ = stencil[ZP] - 2.0*stencil[O] + stencil[ZM];

  D_XY = stencil[XYPP] - stencil[XYMP] - stencil[XYPM] + stencil[XYMM];
  D_XZ = stencil[XZPP] - stencil[XZMP] - stencil[XZPM] + stencil[XZMM];
  D_YZ = stencil[YZPP] - stencil[YZMP] - stencil[YZPM] + stencil[YZMM];

   
/*  for(int i=0; i < 3; ++i){
    D_XX[i] = (stencil[XP])[i] - 2.0*(stencil[O])[i] + (stencil[XM])[i];
    D_YY[i] = (stencil[YP])[i] - 2.0*(stencil[O])[i] + (stencil[YM])[i];
    D_ZZ[i] = (stencil[ZP])[i] - 2.0*(stencil[O])[i] + (stencil[ZM])[i];

    D_XY[i] = (stencil[XYPP])[i] - (stencil[XYMP])[i] - (stencil[XYPM])[i] + (stencil[XYMM])[i];
    D_XZ[i] = (stencil[XZPP])[i] - (stencil[XZMP])[i] - (stencil[XZPM])[i] + (stencil[XZMM])[i];
    D_YZ[i] = (stencil[YZPP])[i] - (stencil[YZMP])[i] - (stencil[YZPM])[i] + (stencil[YZMM])[i];
  }
*/

  cmplx E_X = ( D_XX[0] + D_XY[1] + D_XZ[2] ) / pow(h_, 2);
  cmplx E_Y = ( D_XY[0] + D_YY[1] + D_YZ[2] ) / pow(h_, 2);
  cmplx E_Z = ( D_XZ[0] + D_YZ[1] + D_ZZ[2] ) / pow(h_, 2);

  Eigen::Vector3cd field(E_X, E_Y, E_Z);

  return field;
}
