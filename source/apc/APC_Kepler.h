//------------------------------------------------------------------------------
//
// APC_Kepler.h
// 
// Purpose:
//
//    Computations dealing with Keplerian orbits
//
//------------------------------------------------------------------------------

#ifndef INC_APC_KEPLER_H
#define INC_APC_KEPLER_H


#include "APC_VecMat3D.h"


//------------------------------------------------------------------------------
//
// EccAnom: computes the eccentric anomaly for elliptic orbits
//
// Input:
//
//   M        Mean anomaly in [rad]
//   e        Eccentricity of the orbit [0,1[
//
// <return>   Eccentric anomaly in [rad]
//
//------------------------------------------------------------------------------
double EccAnom (double M, double e);


//------------------------------------------------------------------------------
//
// HypAnom: computes the eccentric anomaly for hyperbolic orbits
//
// Input:
//
//   Mh       Mean anomaly in [rad]
//   e        Eccentricity of the orbit (>1)
//
// <return>   Eccentric anomaly in [rad]
//
//------------------------------------------------------------------------------
double HypAnom (double Mh, double e);


//------------------------------------------------------------------------------
//
// Ellip: computes position and velocity vectors for elliptic orbits
//
// Input:
//
//   GM       Product of gravitational constant and centre mass in [AU^3*d^-2]
//   M        Mean anomaly in [rad]
//   a        Semi-major axis of the orbit in [AU]
//   e        Eccentricity of the orbit (<1)
//
// Output:
//
//   r        Position w.r.t. orbital plane in [AU]
//   v        Velocity w.r.t. orbital plane in [AU/d]
//
//------------------------------------------------------------------------------
void Ellip ( double GM, double M, double a, double e,
			 Vec3D& r, Vec3D& v );


//------------------------------------------------------------------------------
//
// Hyperb: computes position and velocity vectors for hyperbolic orbits
//
// Input:
//
//   GM       Product of gravitational constant and centre mass in [AU^3*d^-2]
//   t0       Time of perihelion passage
//   t        Time for calculation
//   a        Semimajor axis of the orbit in [AU]
//   e        Eccentricity of the orbit (>1)
//
// Output:
//
//   r        Position w.r.t. orbital plane in [AU]
//   v        Velocity w.r.t. orbital plane in [AU/d]
//
// Note: t0 and t in Julian centuries since J2000
//
//------------------------------------------------------------------------------
void Hyperb ( double GM, double t0, double t, double a, double e,
              Vec3D& r, Vec3D& v );


//------------------------------------------------------------------------------
//
// Stumpff: computes values for the Stumpff functions C1, C2 and C3
//
// Input:
//
//   E2       Square of eccentric anomaly (E2=E*E) in [rad^2]
//
// Output:
//
//   c1       Value of C1 = sin(E)/E
//   c2       Value of C2 = (1-cos(E))/(E*E)
//   c3       Value of C3 = (E-sin(E))/(E^3)
//
//------------------------------------------------------------------------------
void Stumpff (double E2, double& c1, double& c2, double& c3);


//------------------------------------------------------------------------------
//
// Parab: computes position and velocity vectors for parabolic and near 
//        parabolic orbits
//
// Input:
//
//   GM       Product of gravitational constant and centre mass in [AU^3*d^-2]
//   t0       Time of perihelion passage
//   t        Time for calculation
//   q        Perihelion distance [AU]
//   e        Eccentricity of the orbit (~1)
//
// Output:
//
//   r        Position w.r.t. orbital plane in [AU]
//   v        Velocity w.r.t. orbital plane in [AU/d]
//
// Note: t0 and t in Julian centuries since J2000
//
//------------------------------------------------------------------------------
void Parab ( double GM, double t0, double t, double q, double e,
             Vec3D& r, Vec3D& v );


//------------------------------------------------------------------------------
//
// Kepler: computes position and velocity vectors for Keplerian orbits w.r.t.
//         the ecliptic
//
// Input:
//
//   GM       Product of gravitational constant and centre mass [AU^3*d^-2]
//   t0       Time of perihelion passage
//   t        Time for calculation
//   q        Perihelion distance [AU]
//   e        Eccentricity of the orbit
//   PQR      Transformation orbital plane -> ecliptic (Gaussian vectors)
//
// Output:
//
//   r        Heliocentric ecliptical position in [AU]
//   v        Heliocentric ecliptical velocity in [AU/d]
//
// Note: t0 and t in Julian centuries since J2000
//
//------------------------------------------------------------------------------
void Kepler ( double GM, double t0, double t, 
              double q, double e, const Mat3D& PQR,
              Vec3D& r, Vec3D& v );


//------------------------------------------------------------------------------
//
// GaussVec: computes the transformation matrix from the orbital plane 
//           coordinate system to the ecliptic
//
// Input:
//
//   Omega    Longitude of the ascending node of the orbit in [rad]
//   i        Inclination of the orbit to the ecliptic in [rad]
//   omega    Argument of perihelion in [rad]
//
// <return>   Transformation matrix containing the Gaussian vectors P, Q and R
//
//------------------------------------------------------------------------------
Mat3D GaussVec (double Omega, double i, double omega);


//------------------------------------------------------------------------------
//
// Elements:  computes the elements of an elliptical orbit from position 
//            and velocity vectors
//
// Input:
//
//   GM       Product of gravitational constant and centre mass [AU^3*d^-2]
//   r        Heliocentric ecliptical position in [AU]
//   v        Heliocentric ecliptical velocity vector in [AU/d]
//
// Output:
//
//   a        Semimajor axis of the orbit [AU]
//   e        Eccentricity of the orbit
//   i        Inclination of the orbit to the ecliptic in [rad]
//   Omega    Longitude of the ascending node of the orbit in [rad]
//   omega    Argument of perihelion in [rad]
//   M        Mean anomaly in [rad]
//
//------------------------------------------------------------------------------
void Elements ( double GM, const Vec3D& r, const Vec3D& v,
                double& a, double& e, double& i,
                double& Omega, double& omega, double& M );


//------------------------------------------------------------------------------
//
// Elements:  computes orbital elements from two given positions 
//
// Input:
//
//   GM       Product of gravitational constant and centre mass [AU^3*d^-2]
//   Mjd_a    Time of the passage of position A [MJD]
//   Mjd_b    Time of the passage of position B [MJD]
//   r_a      Heliocentric ecliptical position A in [AU]
//   r_b      Heliocentric ecliptical position B in [AU]
//
// Output:
//
//   Mjd_p    Time of perihelion passage
//   q        Perihelion distance in [AU]
//   e        Eccentricity of the orbit
//   i        Inclination of the orbit to the ecliptic in [rad]
//   Omega    Longitude of the ascending node of the orbit in [rad]
//   omega    Argument of perihelion in [rad]
//
//------------------------------------------------------------------------------
void Elements ( double GM, double Mjd_a, double Mjd_b,
                const Vec3D& r_a, const Vec3D& r_b,
                double& Mjd_p, double& q, double& e,
                double& i, double& Omega, double& omega );


//------------------------------------------------------------------------------
//
// FindEta: computes the sector-triangle ratio for two given positions and 
//          the time between them
//
// Input:
//
//   r_a      Position at first instant in [AU]
//   r_a      Position at second instant in [AU]
//   tau      Time between the positions r_a and r_b (kGauss * dT in [d])
//
// <return>   Sector-triangle ratio
//
//------------------------------------------------------------------------------
double FindEta (const Vec3D& r_a, const Vec3D& r_b, double tau);


#endif  // include blocker
