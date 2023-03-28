//------------------------------------------------------------------------------
//
// APC_Planets.h
// 
//    Calculation of planetary positions
//
//------------------------------------------------------------------------------


#ifndef INC_APC_PLANETS_H
#define INC_APC_PLANETS_H

#include "APC_VecMat3D.h"


//
// Types
//
enum PlanetType { 
  Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto };


//------------------------------------------------------------------------------
//
// KepPosition: Planetary position from Keplerian elements
//
// Input:
//
//   Planet    Identifies the planet
//   T         Time in Julian centuries since J2000
//
// <return>:   Heliocentric position [AU], true-of-date ecliptic and equinox
//
//------------------------------------------------------------------------------
Vec3D KepPosition (PlanetType Planet, double T);


//------------------------------------------------------------------------------
//
// KepVelocity: Planetary velocity from Keplerian elements
//
// Input:
//
//   Planet    Identifies the planet
//   T         Time in Julian centuries since J2000
//
// <return>    Heliocentric velocity [AU/d], true-of-date ecliptic and equinox
//
//------------------------------------------------------------------------------
Vec3D KepVelocity (PlanetType Planet, double T);


//------------------------------------------------------------------------------
//
// <Planet>Pos: Planetary position from analytical series
//
// Input:
//
//   T         Time in Julian centuries since J2000
//
// <return>:   Heliocentric position [AU], true-of-date ecliptic and equinox
//
//------------------------------------------------------------------------------
Vec3D MercuryPos (double T);
Vec3D VenusPos   (double T);
Vec3D MarsPos    (double T);
Vec3D JupiterPos (double T);
Vec3D SaturnPos  (double T);
Vec3D UranusPos  (double T);
Vec3D NeptunePos (double T);
Vec3D PlutoPos   (double T);


//------------------------------------------------------------------------------
//
// PertPosition: Planetary position from analytical perturbation series
//
// Input:
//
//   Planet    Identifies the planet
//   T         Time in Julian centuries since J2000
//
// <return>:   Heliocentric position [AU], true-of-date ecliptic and equinox
//
//------------------------------------------------------------------------------
Vec3D PertPosition (PlanetType Planet, double T);

#endif   // include blocker