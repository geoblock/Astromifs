//------------------------------------------------------------------------------
//
// APC_PrecNut.h
// 
// Precession and nutation routines
//
//------------------------------------------------------------------------------

#ifndef INC_APC_PRECNUT_H
#define INC_APC_PRECNUT_H


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
Mat3D PrecMatrix_Ecl (double T1, double T2);


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
Mat3D PrecMatrix_Equ (double T1, double T2);


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
Mat3D NutMatrix (double T);


#endif  // include blocker
