//------------------------------------------------------------------------------
//
// APC_PrecNut.cpp
// 
// Implements precession and nutation routines
//
//------------------------------------------------------------------------------

#include <cmath>

#include <math.h>
#include <vcl.h>
#include <windows.h>    // Header file for windows
#include <stdio.h>      // Header file for standard Input/Output


#include "APC_Const.h"
#include "APC_Math.h"
#include "APC_PrecNut.h"
#include "APC_VecMat3D.h"


//------------------------------------------------------------------------------
//
// PrecMatrix_Ecl: Precession of ecliptic coordinates
//
// Input:
//
//   T1        Epoch given
//   T2        Epoch to precess to
//
// <return>:   Precession transformation matrix
//
// Note: T1 and T2 in Julian centuries since J2000
//
//------------------------------------------------------------------------------
Mat3D PrecMatrix_Ecl (double T1, double T2)
{
  //
  // Constants
  //
  const double dT = T2-T1;
  

  //
  // Variables
  //
  double Pi, pi, p_a;

  
  Pi  = 174.876383889*Rad + 
        ( ((3289.4789+0.60622*T1)*T1) +
              ((-869.8089-0.50491*T1) + 0.03536*dT)*dT )/Arcs;
  pi  = ( (47.0029-(0.06603-0.000598*T1)*T1)+
              ((-0.03302+0.000598*T1)+0.000060*dT)*dT )*dT/Arcs;
  p_a = ( (5029.0966+(2.22226-0.000042*T1)*T1)+
              ((1.11113-0.000042*T1)-0.000006*dT)*dT )*dT/Arcs;
    
  return R_z(-(Pi+p_a)) * R_x(pi) * R_z(Pi);
}    


//------------------------------------------------------------------------------
//
// PrecMatrix_Equ: Precession of equatorial coordinates
//
// Input:
//
//   T1        Epoch given
//   T2        Epoch to precess to
//
// <return>:   Precession transformation matrix
//
// Note: T1 and T2 in Julian centuries since J2000
//
//------------------------------------------------------------------------------
Mat3D PrecMatrix_Equ (double T1, double T2)
{
  //
  // Constants
  //
  const double dT = T2-T1;

  
  //
  // Variables
  //
  double zeta,z,theta;


  zeta  =  ( (2306.2181+(1.39656-0.000139*T1)*T1)+
                ((0.30188-0.000344*T1)+0.017998*dT)*dT )*dT/Arcs;
  z     =  zeta + ( (0.79280+0.000411*T1)+0.000205*dT)*dT*dT/Arcs;
  theta =  ( (2004.3109-(0.85330+0.000217*T1)*T1)-
                ((0.42665+0.000217*T1)+0.041833*dT)*dT )*dT/Arcs;

  return R_z(-z) * R_y(theta) * R_z(-zeta);
}    


//------------------------------------------------------------------------------
//
// NutMatrix: Transformation from mean to true equator and equinox
//
// Input:
//
//   T         Time in Julian centuries since epoch J2000
//
// <return>:   Nutation matrix
//
//------------------------------------------------------------------------------
Mat3D NutMatrix (double T)
{
  //
  // Variables
  //
  double  ls, D, F, N;
  double  eps, dpsi, deps;

  // Mean arguments of lunar and solar motion
  ls = pi2*Frac(0.993133+  99.997306*T);   // mean anomaly Sun          
  D  = pi2*Frac(0.827362+1236.853087*T);   // diff. longitude Moon-Sun  
  F  = pi2*Frac(0.259089+1342.227826*T);   // mean argument of latitude 
  N  = pi2*Frac(0.347346-   5.372447*T);   // longit. ascending node    

  // Nutation angles
  dpsi = ( -17.200*sin(N)   - 1.319*sin(2*(F-D+N)) - 0.227*sin(2*(F+N))
           + 0.206*sin(2*N) + 0.143*sin(ls) ) / Arcs;
  deps = ( + 9.203*cos(N)   + 0.574*cos(2*(F-D+N)) + 0.098*cos(2*(F+N))
           - 0.090*cos(2*N)                 ) / Arcs;

  // Mean obliquity of the ecliptic 
  eps  = 0.4090928-2.2696E-4*T;

  return  R_x(-eps-deps)*R_z(-dpsi)*R_x(+eps);
}
