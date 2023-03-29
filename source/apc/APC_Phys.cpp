//------------------------------------------------------------------------------
//
// APC_Phys.cpp
// 
//    Implements physical ephemeris of the Sun and the planets
//
//------------------------------------------------------------------------------

#include <cmath>
#include <math.h>
#include <vcl.h>
#include <windows.h>    // Header file for windows
#include <stdio.h>      // Header file for standard Input/Output


#include "APC_Const.h"
#include "APC_Math.h"
#include "APC_Phys.h"
#include "APC_Planets.h"
#include "APC_VecMat3D.h"


//------------------------------------------------------------------------------
//
// PosAng: Computes the position angle of direction d w.r.t. the position r
//
// Input:
//
//   r         Position of observed body
//   d         Direction for which the position angle is to be calculated
//
// <return>:   Position angle in [rad]; range [0, 2*pi]
//
// Notes: r and d must refer to the same coordinate system and equinox
//
//------------------------------------------------------------------------------
double PosAng (const Vec3D& r, const Vec3D& d)
{
  //
  // Variables
  //
  Vec3D   e_1, e_2, e_3;
  double  c, s, phi;


  // Unit vectors in radial, East, North direction
  e_1 = r / Norm(r);
  e_2 = Cross ( Vec3D(0, 0, 1), r );
  e_2 = e_2 / Norm(e_2);
  e_3 = Cross (e_1, e_2);

  // Position angle
  c = Dot(d, e_3);
  s = Dot(d, e_2);
  phi = atan2(s, c);
  
  return Modulo(phi, 2.0*pi);
}


//------------------------------------------------------------------------------
//
// Shape: Gives equatorial radius and flattening of the planets
//
// Input:
//
//   Planet    Identifies the planet (see APC_Planets.h)
//
// Output:
//
//   R_equ     Equatorial radius in [km]
//   f         Geometrical flattening
//
//------------------------------------------------------------------------------
void Shape (PlanetType Planet, double& R_equ, double& f)
{
  switch (Planet) {
    case Sun:      R_equ =696000.0;  f = 0.0;        break;
    case Mercury:  R_equ =  2439.0;  f = 0.0;        break;
    case Venus:    R_equ =  6051.0;  f = 0.0;        break;
    case Earth:    R_equ =  6378.14; f = 0.00335281; break;
    case Mars:     R_equ =  3393.4;  f = 0.0051865;  break;
    case Jupiter:  R_equ = 71398.0;  f = 0.0648088;  break;
    case Saturn:   R_equ = 60000.0;  f = 0.1076209;  break;
    case Uranus:   R_equ = 25400.0;  f = 0.030;      break;
    case Neptune:  R_equ = 24300.0;  f = 0.0259;     break;
    case Pluto:    R_equ =  1500.0;  f = 0.0;
  }
}


//------------------------------------------------------------------------------
//
// Orient: Calculates orientation of planetocentric reference frame and
//         sense of rotation
//
// Input:
//
//   Planet    Identifies the planet (see APC_Planets.h)
//   System    Specifies rotation system (applies to Jupiter and Saturn only)
//   T         Time in Julian centuries since epoch J2000
//
// Output:
//
//   E         Transformation matrix from Earth mean equator and equinox of
//             J2000 to body fixed equator and prime meridian system
//   Sense     Sense of rotation (direct or retrograde)
//
// Note:
//
//   The current implementation is based on the IAU 1994 rotational elements 
//   (Davies et al. 1996). The implicit light-time correction of the solar
//   prime meridian longitude W0 has been removed for consistency with the 
//   planetary rotation elements.
//
// References:
//
//   Davies M.E., Abalakin V.K., Bursa M., Lieske J.H., Morando B.,
//   Morrison D., Seidelmann P.K., Sinclair A.T., Yallop B., Tjuflin Y.S.;
//   Report of the IAU/IAG/COSPAR Working Group on Cartographic Coordinates
//   and Rotational Elements of the Planets and Satellites: 1994;
//   Celestial Mechanics and Dynamical Astronomy 63, 127-148 (1996)
//
//------------------------------------------------------------------------------
void Orient ( PlanetType Planet, SystemType System, double T,
              Mat3D& E, RotationType& Sense )
{
  //
  // Variables
  //
  double N, Ra, Dec, W;
  double d = 36525.0*T;
    

  // Right ascension and declination of the axis of rotation 
  // with respect to the equator and equinox of J2000;   
  // Prime meridian orientation
  switch (Planet) {
    case Sun:      Ra  = 286.13; 
                   Dec =  63.87; 
                   W   =  84.182 +  14.1844000*d; break;
    case Mercury:  Ra  = 281.01  -   0.033*T; 
                   Dec =  61.45  -   0.005*T;    
                   W   = 329.68  +   6.1385025*d; break;
    case Venus:    Ra  = 272.76;
                   Dec =  67.16;
                   W   = 160.20  -   1.4813688*d; break;
    case Earth:    Ra  =   0.00  -   0.641*T; 
                   Dec =  90.00  -   0.557*T;
                   W   = 190.16  + 360.9856235*d; break;
    case Mars:     Ra  = 317.681 -   0.108*T;
                   Dec =  52.886 -   0.061*T;
                   W   = 176.901 + 350.8919830*d; break;
    case Jupiter:  Ra  = 268.05  -   0.009*T; 
                   Dec =  64.49  +   0.003*T;
                   switch (System) {
                     case Sys_I   : W =  67.10  + 877.900*d; break;
                     case Sys_II  : W =  43.30  + 870.270*d; break;
                     case Sys_III : W = 284.695 + 870.536*d; break;
                   }
                   break;
    case Saturn:   Ra  =  40.589 -   0.036*T;    
                   Dec =  83.537 -   0.004*T;
                   switch (System) {
                     case Sys_I   :
                     case Sys_II  : W = 227.2037 + 844.3000000*d; break;
                     case Sys_III : W =  38.90   + 810.7939024*d; break;
                   }
                   break;
    case Uranus:   Ra  = 257.311;
                   Dec = -15.175;
                   W   = 203.81  - 501.1600928*d; break; // System III
    case Neptune:  N   = Rad*(357.85+52.316*T);    
                   Ra  = 299.36  + 0.70*sin(N); 
                   Dec =  43.46  - 0.51*cos(N); 
                   W   = 253.18  + 536.3128492*d - 0.48*sin(N); break;
    case Pluto:    Ra  = 313.02; 
                   Dec =   9.09;
                   W   = 236.77  -  56.3623195*d;
  }

  Ra*=Rad; Dec*=Rad; W=Rad*Modulo(W,360.0);  
  
  
  // Transformation from Earth mean equator and equinox of J2000 to
  // body fixed equator and prime meridian system                   
  E = R_z(W) * R_x(pi/2.0-Dec) * R_z(pi/2.0+Ra);

  
  // Sense of rotation
  if ( Planet==Venus || Planet==Uranus || Planet==Pluto )
    Sense = Retrograde;
  else
    Sense = Direct;
}  


//------------------------------------------------------------------------------
//
// Rotation: Calculates rotational parameters for the planets
//
// Input:
//
//   r         Planetocentric coordinates of the Earth referred to the Earth 
//             equator and equinox J2000
//   E         Transformation matrix from Earth mean equator and equinox of
//             J2000 to body fixed equator and prime meridian system
//   Sense     Sense of rotation (direct or retrograde)
//   f         Geometrical flattening
//
// Output:
//
//   lon       Planetographic longitude of the Earth in [rad]
//   lat_g     Planetographic latitude of the Earth in [rad]
//   lat_c     Planetocentric latitude of the Earth in [rad]
//
//------------------------------------------------------------------------------
void Rotation ( const Vec3D& r, const Mat3D& E, RotationType Sense, double f,
                double& lon, double& lat_g, double& lat_c )
{
  //
  // Variables
  //
  Vec3D s;

  
  // Planetocentric latitude and longitude 
  s = E*r;

  lat_c = s[theta];
  lon   = s[phi];

  
  // Planetographic latitude and longitude 
  if (Sense==Direct) lon=-lon;
  lon = Modulo(lon,2*pi);
  lat_g = atan ( tan(lat_c) / ((1.0-f)*(1.0-f)) );
}


//------------------------------------------------------------------------------
//
// Illum: Calculates illumination parameters for the planets
//
// Input:
//
//   r         Heliocentric position of body
//   r_Earth   Heliocentric position of the Earth
//
// Output:
//
//   Elong     Elongation in [rad]
//   phi       Phase angle in [rad]
//   k         Phase
//
// Notes: r and r_Earth must refer to the same coordinate system and equinox
//
//------------------------------------------------------------------------------
void Illum ( const Vec3D& r, const Vec3D& r_Earth, 
             double& Elong, double& phi, double& k )
{
  //
  // Variables
  //
  Vec3D   r_geoc;
  double  R,RE,D,c_phi;

  
  // Geocentric position of the planet
  r_geoc = r - r_Earth;

  
  // Distances
  R  = Norm(r);         // Sun-planet distance
  RE = Norm(r_Earth);   // Sun-Earth  distance
  D  = Norm(r_geoc);    // Earth-planet distance


  // Elongation, phase angle and phase
  Elong = acos ( (D*D + RE*RE - R*R) / (2.0*D*RE) );

  c_phi = (D*D + R*R - RE*RE) / (2.0*D*R);
  phi   = acos (c_phi);
  k     = 0.5*(1.0+c_phi);
}


//------------------------------------------------------------------------------
//
// Bright: Calculates the visual magnitude of a planet
//
// Input:
//
//   Planet    Identifies the planet (see APC_Planets.h)
//   r         Heliocentric distance of the planet  in [AU]
//   Delta     Distance of the planet from the Earth  in [AU]
//   phi       Phase angle in [rad]
//   lat       Planetocentric latitude of the Sun in [rad]
//   Dlong     Difference of planetocentric longitudes of the Sun and the 
//             Earth in [rad]
// <return>:   Brightness of the planet in stellar magnitudes [mag]
//
// Notes: lat and Dlong must be given for Saturn only and may be omitted
//        otherwise
//
//------------------------------------------------------------------------------
double Bright ( PlanetType Planet, double r, double Delta, double phi, 
                double lat, double Dlong )
{
  //
  // Variables
  //
  double  p, sd, dl, mag;

  // Normalized phase angle
  p = Deg*phi/100.0;

  // Brightness at unit distance
  switch ( Planet ) {
    case Sun:      return 0.0;
    case Mercury:  mag = -0.42+(3.80-(2.73-2.00*p)*p)*p; break;
    case Venus:    mag = -4.40+(0.09+(2.39-0.65*p)*p)*p; break;
    case Earth:    mag = -3.86;       break;
    case Mars:     mag = -1.52+1.6*p; break;
    case Jupiter:  mag = -9.40+0.5*p; break;
    case Saturn:
      sd = fabs(sin(lat));
      dl = fabs(Modulo(Deg*Dlong+180.0,360.0)-180.0)/100.0;
      mag = -8.88 - 2.60*sd + 1.25*sd*sd + 4.40*dl;
      break;
    case Uranus:   mag = -7.19; break;
    case Neptune:  mag = -6.87; break;
    case Pluto:    mag = -1.0;
  }

  // Apparent brightness
  return ( mag + 5.0*log10(r*Delta) );
}
