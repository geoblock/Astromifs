//------------------------------------------------------------------------------
//
// Phys.cpp
//
// Physical ephemeris of the major planets and the Sun
//
//------------------------------------------------------------------------------

#include <cmath>
#include <iomanip>
#include <iostream>

#include "APC_Const.h"
#include "APC_Math.h"
#include "APC_Phys.h"
#include "APC_Planets.h"
#include "APC_PrecNut.h"
#include "APC_Spheric.h"
#include "APC_Sun.h"
#include "APC_Time.h"
#include "APC_VecMat3D.h"

#ifdef __GNUC__  // GNU C++ adaptation
#include "GNU_iomanip.h"
#endif

using namespace std;


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------
int main()
{
  //
  // Constants
  //
  
  // Body names
  const char* Name[] = 
    { "Sun", "Mercury", "Venus", "Earth", "Mars",
      "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto" };

  
  //
  // Variables
  //
  PlanetType    Planet;
  int           iPlanet, year, month, day;
  double        hour, T, T0;
  Vec3D         r, r_Earth, r_geoc;
  double        Delta;
  double        R_equ, D_equ, f;
  Mat3D         E_I, E_II, E_III;         
  Mat3D         P;                      
  Vec3D         d; 
  RotationType  Sense;
  double        long_I, long_II, long_III, lat_g, lat_c;
  double        long_Sun, latg_Sun, latc_Sun;
  double        Elong, phi, k, mag;
  double        PosAng_Ax, PosAng_Sun;


  // Title
  cout << endl
       << " PHYS: physical ephemeris of the planets and the Sun" << endl
       << "       (c) 1999 Oliver Montenbruck, Thomas Pfleger    " << endl
       << endl;

  // Prompt user for input
  cout << " Date (yyyy mm dd hh.hhh) ... ";
  cin >> year >> month >> day >> hour; cin.ignore(81,'\n');

  T = (Mjd(year,month,day)+hour/24.0-MJD_J2000)/36525.0;

  
  // Header
  cout << endl
       << "            D      V      i    PA(S)   PA(A)" 
       << "    L(I)  L(II) L(III)     B" << endl
       << "         '  \"     mag     o      o       o  "
       << "     o      o      o       o" << endl;


  // Precession (mean-of-date to EME2000)
  P = PrecMatrix_Equ(T,T_J2000);

  // Heliocentric position of the Earth, mean equator and equinox of date
  r_Earth = Ecl2EquMatrix(T) * PertPosition(Earth,T);


  // Loop over bodies
  for (iPlanet=Sun; iPlanet<=Pluto; iPlanet++) {

    Planet = (PlanetType) iPlanet;
    if (Planet==Earth) continue;  // Skip Earth

    // Geometric geocentric position and light time
    r = Ecl2EquMatrix(T) * PertPosition(Planet,T);
    Delta = Norm(r-r_Earth);
    T0 = T - (Delta/c_light)/36525.0,

    // Antedated planetary position at emission of light time
    r = Ecl2EquMatrix(T) * PertPosition(Planet,T0);

    // Geocentric position
    r_geoc = r - r_Earth;

    // Apparent semi-diameter (in ["])
    Shape ( Planet, R_equ,f );
    D_equ = Arcs*2.0*asin(R_equ/(Delta*AU));

    // Transformation from mean equator and equinox of date 
    // to body-fixed frame
    Orient ( Planet, Sys_I  ,T0, E_I  , Sense );  E_I   = E_I   * P;
    Orient ( Planet, Sys_II ,T0, E_II , Sense );  E_II  = E_II  * P;
    Orient ( Planet, Sys_III,T0, E_III, Sense );  E_III = E_III * P;

    // Planetocentric longitude and latitude of the Earth
    Rotation ( -r_geoc, E_I,   Sense, f, long_I,  lat_g,lat_c );
    Rotation ( -r_geoc, E_II,  Sense, f, long_II, lat_g,lat_c );
    Rotation ( -r_geoc, E_III, Sense, f, long_III,lat_g,lat_c );
    
    // Position angle of the rotation axis
    PosAng_Ax  = PosAng ( r_geoc, Row(E_I,z) );

    // Specific computations for Sun and planets
    if (Planet==Sun) {

      // Sun specific convention for longitude and position angle
      // of the rotation axis

      long_I = 2.0*pi-long_I;
      if (PosAng_Ax>pi) PosAng_Ax-=2.0*pi;
    }
    else  { 

      // Ilumination and apparent magnitude
      Illum ( r, r_Earth, Elong,phi,k );

      if (Planet==Saturn) {
        // Planetocentric longitude and latitude of the Sun
        Rotation ( -r, E_I, Sense, f, long_Sun, latg_Sun, latc_Sun );
        mag = Bright ( Planet, Norm(r), Norm(r_geoc), phi,
                       latc_Sun, long_Sun-long_I );
      }
      else
        mag = Bright ( Planet, Norm(r),Norm(r_geoc), phi);

      // Position angle of the Sun
      PosAng_Sun = PosAng ( r_geoc, -r );
    };


    // Output
    if (Planet==Sun) 
      cout << " " << left << setw(6) << Name[Planet] << right << fixed
           << setprecision(2) << setw(8) << Angle(D_equ/60.0,DMMm) 
           << setprecision(2) << setw(29) << Deg*PosAng_Ax; 
    else
      cout << " " << left << setw(9) << Name[Planet] << right << fixed 
           << setprecision(2) << setw(5) << D_equ
           << setprecision(1) << setw(6) << mag
           << setprecision(1) << setw(7) << Deg*phi
           << setprecision(2) << setw(8) << Deg*PosAng_Sun
           << setprecision(2) << setw(8) << Deg*PosAng_Ax;

    // Print planetographic longitude(s) of the Earth
    switch (Planet) {

      // Bodies with system I only  
      case Sun: case Mercury: case Venus: case Mars: case Pluto:
        cout << setprecision(2) << setw(8) << Deg*long_I 
             << setw(14) << " ";
        break;
      
      case Jupiter:
        cout << setprecision(2) << setw(8) << Deg*long_I 
             << setprecision(2) << setw(7) << Deg*long_II
             << setprecision(2) << setw(7) << Deg*long_III;
        break;
      
      case Saturn:
        cout << setprecision(2) << setw(8)  << Deg*long_I 
             << setprecision(2) << setw(14) << Deg*long_III;
        break;
      
      // Bodies with system III only  
      case Uranus: case Neptune:
        cout << setprecision(2) << setw(22) << Deg*long_III;
    }
    
    cout << setprecision(2) << setw(8)  << Deg*lat_g << endl;
 
   } // End Planet loop

}
