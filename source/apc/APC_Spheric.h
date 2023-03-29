//------------------------------------------------------------------------------
//
// APC_Spheric.h
// 
// Transformations frequently found in spherical astronomy
//
//------------------------------------------------------------------------------

#ifndef INC_APC_SPHERIC_H
#define INC_APC_SPHERIC_H


#include "APC_VecMat3D.h"


//------------------------------------------------------------------------------
//
// Equ2EclMatrix: Transformation of equatorial to ecliptical coordinates
//
// Input:
//
//   T         Time in Julian centuries since J2000
//
// <return>:   Transformation matrix
//
//------------------------------------------------------------------------------
Mat3D Equ2EclMatrix (double T);


//------------------------------------------------------------------------------
//
// Ecl2EquMatrix: Transformation of ecliptical to equatorial coordinates
//
// Input:
//
//   T         Time in Julian centuries since J2000
//
// <return>:   Transformation matrix
//
//------------------------------------------------------------------------------
Mat3D Ecl2EquMatrix (double T);


//------------------------------------------------------------------------------
//
// Equ2Hor: Transformation of equatorial coordinates to the horizon system
//
// Input:
//
//   Dec       Declination
//   tau       Hour angle
//   lat       Geographical latitude of the observer
//
// Output :
//
//   h         Altitude
//   Az        Azimuth
//
// Note: all parameters in [rad]
//
//------------------------------------------------------------------------------
void Equ2Hor ( double Dec, double tau, double lat, 
               double& h, double& Az );


//------------------------------------------------------------------------------
//
// Hor2Equ: Transformation of horizon system coordinates 
//          to equatorial coordinates
//
// Input:
//
//   h         Altitude
//   Az        Azimuth
//   lat       Geographical latitude of the observer
//
// Output:
//
//   Dec       Declination
//   tau       Hour angle
//
// Note: all parameters in [rad]
//
//------------------------------------------------------------------------------
void Hor2Equ ( double h, double Az, double lat,
               double& Dec, double& tau );


//------------------------------------------------------------------------------
//
// Site: Calculates the geocentric position of a site on the Earth's surface
//
// Input:
//
//   lambda    Geographical longitude (east positive) in [rad]
//   phi       Geographical latitude  in [rad]
//
// <return>:   Geocentric position in [km]
//
//------------------------------------------------------------------------------
Vec3D Site (double lambda, double phi);


//------------------------------------------------------------------------------
//
// StdEqu: Calculation of equatorial coordinates from standard coordinates
//         used in astrographic plate measurement
//
// Input:
//
//   RA0       Right ascension of optical axis
//   Dec0      Declination of optical axis
//   X         Standard coordinate X [no dimension]
//   Y         Standard coordinate Y [no dimension]
//
// Output :
//
//   Ra        Right ascension
//   Dec       Declination
//
// Note: all angles in [rad]
//
//------------------------------------------------------------------------------
void StdEqu ( double RA0, double Dec0, double X, double Y,
              double& Ra, double& Dec );


//------------------------------------------------------------------------------
//
// EquStd: Calculation of standard coordinates from equatorial coordinates
//
// Input:
//
//   RA0       Right ascension of optical axis
//   Dec0      Declination of optical axis
//   Ra        Right ascension
//   Dec       Declination
//
// Output:
//
//   X         Standard coordinate X [no dimension]
//   Y         Standard coordinate Y [no dimension]
//
// Note: all angles in [rad]
//
//------------------------------------------------------------------------------
void EquStd ( double RA0, double Dec0, double Ra, double Dec,
              double& X, double& Y );


#endif  // include blocker
