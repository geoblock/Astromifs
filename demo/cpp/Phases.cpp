//------------------------------------------------------------------------------
//
// Phases program (Phases.cpp)
//
// Determines the times of lunar phases and checks for eclipses
//
//------------------------------------------------------------------------------

#include <cmath>
#include <iomanip>
#include <iostream>

#include "APC_Const.h"
#include "APC_Math.h"
#include "APC_Moon.h"
#include "APC_Sun.h"
#include "APC_Time.h"
#include "APC_VecMat3D.h"

#ifdef __GNUC__  // GNU C++ adaptation
#include "GNU_iomanip.h"
#endif

using namespace std;


//
// Types
//
enum enLunarPhase { NewMoon, FirstQuarter, FullMoon, LastQuarter };


//
// Global variables
//
enLunarPhase Phase = NewMoon;


//------------------------------------------------------------------------------
//
// PhasesFunc: Goal function for search of phase events [-pi, pi]
//
// Input:
//
//   T         Ephemeris Time (Julian centuries since J2000)
//
// <return>:   Difference between the longitude of the Moon from the Sun
//             and the nominal value for a given phase (New Moon 0, First 
//             Quarter pi/2, etc.) (in [rad])
//
// Global:
//
//   Phase     Flag for desired lunar phase
//
//------------------------------------------------------------------------------
double PhasesFunc (double T)
{
  //
  // Constants
  //
  const double tau_Sun  = 8.32 / (1440.0*36525.0);    // 8.32 min  [cy]

  
  //
  // Variables
  //
  const double LongDiff = MoonPos(T)[phi] - SunPos(T-tau_Sun)[phi];

  
  return  Modulo ( LongDiff - Phase*pi/2.0 + pi, pi2 ) - pi;
}


//------------------------------------------------------------------------------
//
// SolarEclipseFlag: Returns a 2 char. string indicating the possibility of 
//                   a solar eclipse
//
// Input:
//
//   beta      Ecliptic latitude of the Moon in [rad]
//
// <return>:   Eclipse flag
//
//------------------------------------------------------------------------------
char* SolarEclipseFlag (double beta)
{
  //
  // Constants
  //
  const double b = fabs(beta);

  
  if (b>0.027586) return "  ";   // No solar eclipse possible

  if (b<0.015223) return "c ";   // Central eclipse certain
  if (b<0.018209) return "c?";   // Possible central eclipse
  if (b<0.024594) return "p ";   // Partial solar eclipse certain
  
  return "p?";                   // Possible partial solar eclipse
}


//------------------------------------------------------------------------------
//
// LunarEclipseFlag: Returns a 2 char. string indicating the possibility of
//                   a lunar eclipse
//
// Input:
//
//   beta      Ecliptic latitude of the Moon [rad]
//
// <return>:   Eclipse flag
//
//------------------------------------------------------------------------------
char* LunarEclipseFlag (double beta)
{
  //
  // Constants
  //
  const double b = fabs(beta);


  if (b>0.028134) return "  ";   // No lunar eclipse possible

  if (b<0.006351) return "t ";   // Total lunar eclipse certain
  if (b<0.009376) return "t?";   // Possible total eclipse
  if (b<0.015533) return "p ";   // Partial lunar eclipse certain
  if (b<0.018568) return "p?";   // Possible partial eclipse
  if (b<0.025089) return "P ";   // Penumbral lunar eclipse certain
  
  return "P?";                   // Possible penumbral lunar eclipse
}


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
  const double dT  = 7.0 / 36525.0;           // Step (1 week)
  const double Acc = (0.5/1440.0) / 36525.0;  // Desired Accuracy (0.5 min)


  //
  // Variables
  //
  bool   valid   = false; // Do we have an approximate value for ET-UT ?
  bool   Success = false; // Flag for success of search for phase time
  int    iPhase, Lunation, Year;
  double MjdUT, ET_UT, T0, T1, TPhase, D0, D1, beta;


  // Title
  cout << endl
       << "        PHASES: Phases of the Moon and Check for Eclipses" << endl
       << "           (c) 1999 Oliver Montenbruck, Thomas Pfleger   " << endl
       << endl;


  // Input
  cout << " Dates of Lunar Phases for the year ... ";
  cin >> Year; cin.ignore(81,'\n');
  cout << endl;


  // Header
  cout << " " << left 
       << setw(19) << "New Moon"  << setw(19) << "First Quarter"
       << setw(19) << "Full Moon" << setw(19) << "Last Quarter" << endl;

  
  // Assure that we get the year's first New Moon
  // [T0, T1]: Interval bracketing the instant TPhase to be found
  T0 = ( Mjd(Year-1, 12, 1) - MJD_J2000 ) / 36525.0;
  T1 = T0 + dT;
  TPhase = T0;
  

  // Check 13 months
  for (Lunation=0; Lunation<=13; Lunation++) {
    
    // Search for phases
    for (iPhase = NewMoon; iPhase <= LastQuarter; iPhase++) {
  
      Phase = (enLunarPhase) iPhase;

      // Bracket desired phase event
      D0 = PhasesFunc(T0);
      D1 = PhasesFunc(T1);
      
      while ( (D0*D1>0.0) || (D1<D0) ) {
        T0=T1; D0=D1; T1+=dT; D1=PhasesFunc(T1);
      }

      // Iterate time of phase
      Pegasus ( PhasesFunc, T0, T1, Acc, TPhase, Success );

      // Correct for difference of ephemeris time and universal time
      ETminUT ( TPhase, ET_UT, valid );
      MjdUT = ( TPhase*36525.0+MJD_J2000) - ET_UT/86400.0;

      // Output date and time of event
      if ( Phase==NewMoon ) cout << ' ';
      cout << DateTime(MjdUT,HHMM); 
      
      // Eclipse check
      beta = MoonPos(TPhase)[theta];
      if      ( iPhase == NewMoon )
        cout << setw(3) << left << SolarEclipseFlag(beta);
      else if ( iPhase == FullMoon)
        cout << setw(3) << left << LunarEclipseFlag(beta);
      else 
        cout << setw(3) << " ";  // First or Last Quarter
      
      // Move search interval by one week
      T0 = TPhase;
      T1 = T0 + dT;
    
    }
    cout << endl;
  
  }
  
  // Footer
  if (valid) 
    cout << endl << " All times in Universal Time (UT)." << endl << endl;
  else
    cout << endl << " All times in Ephemeris Time (ET)." << endl << endl;
}
