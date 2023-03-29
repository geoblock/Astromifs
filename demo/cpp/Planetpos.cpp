//-------------------------------------------------------------------------
//
// Planetpos.cpp
//
// Ephemeris of the major planets
//
//-------------------------------------------------------------------------

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "APC_Const.h"
#include "APC_Kepler.h"
#include "APC_Math.h"
#include "APC_Planets.h"
#include "APC_PrecNut.h"
#include "APC_Spheric.h"
#include "APC_Sun.h"
#include "APC_Time.h"
#include "APC_VecMat3D.h"

#ifdef __GNUC__  // GNU C++ adaptation
#include <ctype.h>
#include "GNU_iomanip.h"
#endif

using namespace std;


//-------------------------------------------------------------------------
//
// Main program
//
//-------------------------------------------------------------------------
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
  PlanetType Planet;
  char       Mode;
  int        iPlanet, last;
  int        year, month, day;
  double     Hour, MJD, T;
  double     dist, fac;
  Vec3D      R_Sun, r_helioc, r_geoc, r_equ;
  Mat3D      P;
  bool       End = false;


  // Title
  cout << endl
       << "   PlanetPos: geocentric and heliocentric planetary positions " << endl
       << "            (c) 1999 Oliver Montenbruck, Thomas Pfleger     " << endl
       << endl;


  // Command loop
  do {

    // Command input
    cout << endl
         << " (J) J2000 astrometric       (B) B1950 astrometric  " << endl
         << " (A) apparent coordinates    (E) exit               " << endl
         << endl
         << " Enter option ... "; 
    cin >> Mode; cin.ignore(81,'\n'); Mode = tolower(Mode);
    cout << endl;

    // Actions
    
    // Terminate, if Mode=Exit
    if (Mode=='e') {
      End=true; 
      break;
    }

    // Cycle, if Mode is invalid
    if (Mode!='j' && Mode!='b' && Mode!='a') continue;

    
    // Epoch date
    cout << " Date (yyyy mm dd hh.hhh) ... "; 
    cin >> year >> month >> day >> Hour; cin.ignore(81,'\n');
    cout << endl << endl << endl;

    MJD = Mjd(year,month,day)+Hour/24.0;
    T   = ( MJD - MJD_J2000 ) / 36525.0;

    // Header
    cout << " Date:  " << DateTime(MJD,HHh) << " (ET)"
         << "     JD: " << (MJD+2400000.5)
         << setw(18) << "equinox ";

    switch (Mode) {
      case 'a': cout << "of date" << endl; break;
      case 'j': cout << "J2000"   << endl; break;
      case 'b': cout << "B1950"   << endl;
    }

    cout << endl
         << endl
         << "               l            b          r          "
         << "Ra          Dec        delta"
         << endl
         << "            o  '  \"     o  '  \"       AU        "
         << "h  m  s      o  '  \"      AU"
         << endl;

   
    // Ecliptic coordinates of the Sun, equinox of date
    R_Sun = SunPos(T);

    // Loop over planets
    if ( (-1.1<T) && (T<1.0) ) 
      last=Pluto; 
    else last=Neptune;
    
    for (iPlanet=Sun; iPlanet<=last; iPlanet++) {
      
      // Heliocentric ecliptic coordinates of the planet; equinox of date
      Planet = (PlanetType) iPlanet;
      r_helioc = PertPosition(Planet, T);

      // Geocentric ecliptic coordinates (equinox of date)
      r_geoc = r_helioc + R_Sun;

      // First-order light-time/aberration correction
      dist = Norm(r_geoc);
      fac = dist/c_light;

      if (Mode=='a') 
        r_geoc -= fac*(KepVelocity(Planet,T)-KepVelocity(Earth,T)); 
      else        
        r_geoc -= fac*KepVelocity(Planet,T); 

      // Precession and equatorial coordinates;

      switch (Mode) {

        case 'a':  r_equ    = NutMatrix(T) * Ecl2EquMatrix(T) * r_geoc; 
                   break;
        case 'j':  P        = PrecMatrix_Ecl(T,T_J2000);
                   r_helioc = P * r_helioc;
                   r_geoc   = P * r_geoc;
                   r_equ    = Ecl2EquMatrix(T_J2000) * r_geoc;
                   break;
        case 'b':  P        = PrecMatrix_Ecl(T,T_B1950);
                   r_helioc = P * r_helioc;
                   r_geoc   = P * r_geoc;
                   r_equ    = Ecl2EquMatrix(T_B1950) * r_geoc;
      } 

      
      // Output
      cout << " " << setw(8) << left << Name[iPlanet] << right;
      cout << fixed << setprecision(1)
           << setw(11) << Angle(Deg*r_helioc[phi],DMMSSs) << " " 
           << setw(11) << Angle(Deg*r_helioc[theta],DMMSSs);
      if (Planet<=Earth) 
        cout << setprecision(6) << setw(11) << r_helioc[r];
      else
        cout << setprecision(5) << setw(10) << r_helioc[r] << " ";
      cout << setprecision(2) << setw(13) << Angle(Deg*r_equ[phi]/15.0,DMMSSs)
           << " " << showpos 
           << setprecision(1) << setw(11) << Angle(Deg*r_equ[theta],DMMSSs)
           << noshowpos;
      if (Planet<=Earth) 
        cout << setprecision(6) << setw(11) << dist;
      else
        cout << setprecision(5) << setw(10) << dist << " ";
      cout << endl;

    };


    // Trailer
    cout << endl
         << " l,b,r:   heliocentric ecliptic (geometric) " << endl
         << " Ra,Dec:  geocentric equatorial ";
    if (Mode=='a') 
      cout << "(apparent)" << endl;
    else 
      cout << "(astrometric)" << endl;
    cout << " delta:   geocentric distance   (geometric)" << endl
         << endl;

  }
  while (!End);
 
}
