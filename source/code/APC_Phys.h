//------------------------------------------------------------------------------
//
// APC_Phys.h
// 
//    Physical ephemeris of the Sun and the planets
//
//------------------------------------------------------------------------------

#ifndef INC_APC_PHYS_H
#define INC_APC_PHYS_H


#include "APC_Planets.h"
#include "APC_VecMat3D.h"


//
// Types
//
enum SystemType   { Sys_I, Sys_II, Sys_III }; // rotation systems (Jup./Sat.)
enum RotationType { Direct, Retrograde };     // sense of rotation


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
double PosAng (const Vec3D& r, const Vec3D& d);


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
void Shape (PlanetType Planet, double& R_equ, double& f);


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
//------------------------------------------------------------------------------
void Orient ( PlanetType Planet, SystemType System, double T,
              Mat3D& E, RotationType& Sense );


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
//   lon       Planetographic longitude of the earth in [rad]
//   lat_g     Planetographic latitude of the earth in [rad]
//   lat_c     Planetocentric latitude of the earth in [rad]
//
//------------------------------------------------------------------------------
void Rotation ( const Vec3D& r, const Mat3D& E, RotationType Sense, double f,
                double& lon, double& lat_g, double& lat_c );


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
             double& Elong, double& phi, double& k );


//------------------------------------------------------------------------------
//
// Bright: Calculates the visual magnitude of a planet
//
// Input:
//
//   Planet    Identifies the planet (see APC_Planets.h)
//   r         Heliocentric distance of the planet in [AU]
//   Delta     Distance of the planet from the Earth in [AU]
//   phi       Phase angle in [rad]
//   lat       Planetocentric latitude of the sun in [rad]
//   Dlong     Difference of planetocentric longitudes of the Sun and the 
//             Earth in [rad]
// <return>:   Brightness of the planet in stellar magnitudes [mag]
//
// Notes: lat and Dlong must be given for Saturn only and may be omitted
//        otherwise
//
//------------------------------------------------------------------------------
double Bright ( PlanetType Planet, double r, double Delta, double phi,
                double lat=0.0, double Dlong=0.0 );


#endif   // include blocker