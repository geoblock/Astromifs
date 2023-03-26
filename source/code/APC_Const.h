//------------------------------------------------------------------------------
//
// APC_Const.h
// 
// Definition of constants used in APC library and applications
//
//------------------------------------------------------------------------------

#ifndef INC_APC_CONST_H
#define INC_APC_CONST_H


//
// Constants
//
const double pi        = 3.14159265358979324;
const double pi2       = 2.0*pi;
const double Rad       = pi / 180.0;
const double Deg       = 180.0 / pi;
const double Arcs      = 3600.0*180.0/pi;

// Radii of Earth, Sun and Moon
const double R_Earth   =   6378.137;     // [km]
const double R_Sun     = 696000.0;       // [km]
const double R_Moon    =   1738.0;       // [km]

const double MJD_J2000 = 51544.5;        // MJD of Epoch J2000.0
const double T_J2000   =  0.0;           // Epoch J2000.0
const double T_B1950   = -0.500002108;   // Epoch B1950

const double kGauss    = 0.01720209895;  // gravitational constant
const double GM_Sun    = kGauss*kGauss;  // [AU^3/d^2]

const double AU        = 149597870.0;    // Astronomical unit [km]

const double c_light   = 173.14;         // speed of light [AU/d]

#endif  // include blocker