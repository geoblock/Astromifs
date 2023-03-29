//------------------------------------------------------------------------------
//
// Luna.cpp
//
// Lunar ephemeris
//
//------------------------------------------------------------------------------

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "APC_Cheb.h"
#include "APC_Const.h"
#include "APC_Math.h"
#include "APC_Moon.h"
#include "APC_Time.h"
#include "APC_VecMat3D.h"

#ifdef __GNUC__  // GNU C++ adaptation
#include "GNU_iomanip.h"
#endif

using namespace std;


//------------------------------------------------------------------------------
//
// GetEph: Prompts user for date of start, step size and end of ephemeris
//
// Output:
//
//   MjdStart  Start time of ephemeris as Modified Julian Date
//   Step      Step size of ephemeris in days
//   MjdEnd    End time of ephemeris as Modified Julian Date
//
//------------------------------------------------------------------------------
void GetEph (double& MjdStart, double& Step, double& MjdEnd)
{
  //
  // Variables
  //
  int    year, month, day;
  double hour;


  cout << endl
       << " Begin and end of the ephemeris: " << endl
       << endl;

  cout << "  first date (yyyy mm dd hh.hhh) ... ";
  cin >> year >> month >> day >> hour; cin.ignore(81,'\n');
  MjdStart = Mjd(year,month,day) + hour/24.0 ;

  cout << "  final date (yyyy mm dd hh.hhh) ... ";
  cin >> year >> month >> day >> hour; cin.ignore(81,'\n');
  MjdEnd = Mjd(year,month,day) + hour/24.0 ;

  cout << "  step size (dd hh.hh)           ... ";
  cin >> day >> hour; cin.ignore(81,'\n');
  Step  = day + hour/24.0;
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
  const int    Degree    = 10;  
  const double Interval  = 10.0/36525.0;  // 10d in [cy]
  
  //
  // Variables
  //
  int    n_line = 0;
  double MjdStart, Step, MjdEnd;
  double Date, T;
  double Dist, Parallax;
  Vec3D  r_Moon;
  Cheb3D ChebMoonEqu (MoonEqu, Degree, Interval);


  // Title
  cout << endl
       << "                     LUNA: lunar ephemeris           " << endl
       << "        (c) 1999 Oliver Montenbruck, Thomas Pfleger  " << endl
       << endl;


  // Prompt user for date of start, step size and end of ephemeris 
  GetEph (MjdStart, Step, MjdEnd);

  
  // Header
  cout << endl << endl
       << "    Date      ET        Ra            Dec      Distance  "
       << " Parallax " << endl
       << "                      h  m  s       o  '  \"  Earth radii"
       << "    '  \"   " << endl;

  // Compute ephemeris
  Date = MjdStart;

  
  // Time loop
  while ( Date < MjdEnd + Step/2 ) {

    // Geocentric coordinates of the Moon from Chebyshev approximation
    // (true equator and equinox of date)
    T = ( Date - MJD_J2000 ) / 36525.0;
    
    r_Moon = ChebMoonEqu.Value(T);


    // Distance and parallax
    Dist     = r_Moon[r];
    Parallax = asin(R_Earth/Dist);


    // Output
    cout << " " << DateTime(Date,HHMM) 
         << fixed 
         << setprecision(2) << setw(14) << Angle(Deg*r_Moon[phi]/15.0,DMMSSs)
         << "  " << showpos 
         << setprecision(1) << setw(11) << Angle(Deg*r_Moon[theta],DMMSSs)
         << noshowpos 
         << setprecision(3) << setw(10) << Dist/R_Earth
         << setprecision(2) << setw(12) << Angle(60.0*Deg*Parallax,DMMm)
         << endl;
    
    ++n_line; if ( (n_line % 5) ==0 ) cout << endl;


    Date += Step; // Next time step

  };  // End time loop

}
