//------------------------------------------------------------------------------
//
// APC_Spheric.cpp
//
// Transformations frequently found in spherical astronomy
//
//------------------------------------------------------------------------------

#include <cmath>

#include <math.h>
#include <stdio.h>      // Header file for standard Input/Output


#include "APC_Const.h"
#include "APC_Spheric.h"
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
Mat3D Equ2EclMatrix (double T)
{
  //
  // Constants
  //
  const double 
    eps = ( 23.43929111-(46.8150+(0.00059-0.001813*T)*T)*T/3600.0 ) * Rad;
   
  return R_x(eps);
}


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
Mat3D Ecl2EquMatrix (double T)
{
  return Transp(Equ2EclMatrix(T));
}


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
// Output:
//
//   h         Altitude
//   Az        Azimuth
//
// Note: all parameters in [rad]
//
//------------------------------------------------------------------------------
void Equ2Hor ( double Dec, double tau, double lat, 
               double& h, double& Az )
{
  //
  // Variables
  //
  Vec3D e_equ, e_hor;

  
  e_equ = Vec3D(Polar(tau,Dec));       // unit vector in horizontal system
  e_hor = R_y(pi/2.0-lat) * e_equ;     // unit vector in equatorial system

  Az = e_hor[phi];                     // polar angles
  h  = e_hor[theta];
}


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
               double& Dec, double& tau )
{
  //
  // Variables
  //
  Vec3D e_equ, e_hor;

  
  e_hor = Vec3D(Polar(Az,h));          // unit vector in horizontal system
  e_equ = R_y(-(pi/2.0-lat)) * e_hor;  // unit vector in equatorial system

  tau = e_equ[phi];                    // polar angles
  Dec = e_equ[theta];
}

    
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
Vec3D Site (double lambda, double phi)
{
  //
  // Constants
  //
  const double f     = 1.0/298.257;   // Flattening 
  const double e_sqr = f*(2.0-f);     // Square of eccentricity
  const double cos_phi = cos(phi);    // (Co)sine of geographical latitude
  const double sin_phi = sin(phi);


  //
  // Variables
  //
  double  N = R_Earth / sqrt (1.0-e_sqr*(sin_phi*sin_phi));
      

  // Cartesian position vector [km]
  return Vec3D ( N*cos_phi*cos(lambda),
                 N*cos_phi*sin(lambda),
                 (1.0-e_sqr)*N*sin_phi  );
}


//------------------------------------------------------------------------------
//
// StdEqu: Calculation of equatorial coordinates from standard coordinates
//         used in astrographic plate measurement
//
// Input:
//
//   RA0       Right ascension of optical axis [rad]
//   Dec0      Declination of optical axis [rad]
//   X         Standard coordinate X
//   Y         Standard coordinate Y
//
// Output:
//
//   Ra        Right ascension [rad]
//   Dec       Declination [rad]
//
// Note: all angles in [rad]
//
//------------------------------------------------------------------------------
void StdEqu ( double RA0, double Dec0, double X, double Y,
              double& Ra, double& Dec )
{
  Ra  = RA0 + atan ( -X / (cos(Dec0)-Y*sin(Dec0)) );
  Dec = asin ( (sin(Dec0)+Y*cos(Dec0))/sqrt(1.0+X*X+Y*Y) );
}


//------------------------------------------------------------------------------
//
// EquStd: Calculation of standard coordinates from equatorial coordinates
//
// Input:
//
//   RA0       Right ascension of optical axis [rad]
//   Dec0      Declination of optical axis [rad]
//   Ra        Right ascension [rad]
//   Dec       Declination [rad]
//
// Output:
//
//   X         Standard coordinate X 
//   Y         Standard coordinate Y 
//
//------------------------------------------------------------------------------
void EquStd ( double RA0, double Dec0, double Ra, double Dec,
              double& X, double& Y ) 
{
  //
  // Constants
  //
  const double c =
    cos(Dec0)*cos(Dec)*cos(Ra-RA0)+sin(Dec0)*sin(Dec);


  X = - ( cos(Dec)*sin(Ra-RA0) ) / c;
  Y = - ( sin(Dec0)*cos(Dec)*cos(Ra-RA0)-cos(Dec0)*sin(Dec) ) / c;
}
