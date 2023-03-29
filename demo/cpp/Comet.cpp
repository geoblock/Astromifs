//------------------------------------------------------------------------------
//
// File:    Comet.cpp
//
// Calculation of unperturbed ephemeris with arbitrary eccentricity
//          for comets and minor planets
//
//------------------------------------------------------------------------------

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "APC_Const.h"
#include "APC_IO.h"
#include "APC_Kepler.h"
#include "APC_Math.h"
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
// GetElm: Reads orbital elements from input file
//
// Input:
//
//   Filename  Name of input file with orbital elements
//
// Output:
//
//   Mjd0      Time of perihelion as Modified Julian Date
//   q         Perihelion distance in [AU]
//   e         Eccentricity
//   PQR       Transformation matrix orbital plane -> ecliptic
//   T_eqx0    Equinox for the orbital elements
//
//------------------------------------------------------------------------------
void GetElm ( char* Filename,
              double& Mjd0, double& q, double& e, Mat3D& PQR, double& T_eqx0 )
{
  //
  // Variables
  //
  int       year, month;
  double    Omega, i, omega;
  double    Year,day;
  ifstream  inp;

  
  inp.open(Filename);   // Open orbital elements file


  // Read elements and dump to screen
  cout << endl
       << " Orbital elements from file " << Filename << endl
       << endl;

  inp >> year >> month >> day; inp.ignore(81,'\n');
  Mjd0  = Mjd(year,month,int(day))+(day-int(day));
  cout << "  perihelion time (y m d) "
       << "  " << setprecision(2) << DateTime(Mjd0,DDd) << endl;

  inp >> q; inp.ignore(81,'\n');
  cout << "  perihelion distance (q) "
       << setprecision(7) << setw(14) << q << " AU" << endl;

  inp >> e; inp.ignore(81,'\n');
  cout << "  eccentricity (e)        "
       << setprecision(7) << setw(14) << e << endl;

  inp >> i; inp.ignore(81,'\n');
  cout << "  inclination (i)         "
       << setprecision(5) << setw(12) << i << " deg" << endl;

  inp >> Omega; inp.ignore(81,'\n');
  cout << "  long. of ascending node "
       << setprecision(5) << setw(12) << Omega << " deg" << endl;

  inp >> omega; inp.ignore(81,'\n');
  cout << "  argument of perihelion  "
       << setprecision(5) << setw(12) << omega << " deg" << endl; 

  inp >> Year; inp.ignore(81,'\n');
  cout << "  equinox                 "
       << setprecision(2) << setw(9)  << Year << endl;
  
  inp.close();
  
  
  T_eqx0 = (Year-2000.0)/100.0;

  PQR = GaussVec(Rad*Omega, Rad*i, Rad*omega);
}


//------------------------------------------------------------------------------
//
// GetEph: Prompts user for date of start and end of ephemeris, step size and
//         desired equinox
//
// Output:
//
//   MjdStart  Start time of ephemeris as Modified Julian Date
//   Step      Step size of ephemeris in days
//   MjdEnd    End time of ephemeris as Modified Julian Date
//   T_eqx     Desired equinox for ephemeris in Jul. cent. since J2000
//
//------------------------------------------------------------------------------
void GetEph (double& MjdStart, double& Step, double& MjdEnd, double& T_eqx)
{
  //
  // Variables
  //
  int    year, month, day;
  double hour, Year;


  cout << endl
       << " Begin and end of the ephemeris: " << endl
       << endl;

  
  // Query user for start and end date, step and desired equinox
  cout << "  first date (yyyy mm dd hh.hhh)            ... ";
  cin >> year >> month >> day >> hour; cin.ignore(81,'\n');
  MjdStart = Mjd(year, month, day) + hour/24.0 ;

  cout << "  final date (yyyy mm dd hh.hhh)            ... ";
  cin >> year >> month >> day >> hour; cin.ignore(81,'\n');
  MjdEnd = Mjd(year, month, day) + hour/24.0 ;

  cout << "  step size (dd hh.hh)                      ... ";
  cin >> day >> hour; cin.ignore(81,'\n');
  Step  = day + hour/24.0;

  cout << endl
       << " Desired equinox of the ephemeris (yyyy.y)  ... ";
  cin >> Year; cin.ignore(81,'\n');
  
  T_eqx = (Year-2000.0)/100.0;
}


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------
void main(int argc, char* argv[]) {

  //
  // Variables
  //
  int       n_line = 0;
  double    MjdStart, Step, MjdEnd, T_eqx;
  double    Mjd0, q, e, T_eqx0;
  Mat3D     PQR;
  double    Date, T;
  double    dist, fac;
  Vec3D     R_Sun, r_helioc, v_helioc, r_geoc, r_equ;
  char      InputFile[APC_MaxFilename]  = "";
  char      OutputFile[APC_MaxFilename] = "";
  bool      FoundInputfile  = false;
  bool      FoundOutputfile = false;
  ofstream  OutFile;


  // Title
  cout << endl
       << " COMET: ephemeris calculation for comets and minor planets" << endl
       << "        (c) 1999 Oliver Montenbruck, Thomas Pfleger       " << endl
       << endl;

  
  // Find input file
  GetFilenames( argc, argv, "Comet.dat", InputFile, FoundInputfile,
                OutputFile, FoundOutputfile );
  
  // Terminate program if input file could not be found
  if ( !FoundInputfile ) {
    cerr << " Terminating program." << endl;
    exit(-1);
  }

  
  // Input orbital elements and prompt user
  GetElm ( InputFile, Mjd0, q, e, PQR, T_eqx0 );
  GetEph ( MjdStart, Step, MjdEnd, T_eqx );


  PQR = PrecMatrix_Ecl(T_eqx0, T_eqx) * PQR; 
  

  // Redirect output if output file shall be created
  if (FoundOutputfile) {
    OutFile.open(OutputFile);
    if (OutFile.is_open())
      cout = OutFile;
  }
  
  
  // Header
  cout << endl << endl
       << "    Date     ET   Sun     l      b     r"
       << "       Ra          Dec      Distance " << endl
       << setw(43) << " " << "   h  m  s      o  '  \"     (AU) " << endl;

  
  // Compute ephemeris
  Date = MjdStart;

  // Time loop
  while ( Date < MjdEnd + Step/2 ) {

    // Geocentric ecliptic coordinates of the Sun, equinox T_eqx
    T = ( Date - 51544.5 ) / 36525.0;
    R_Sun = PrecMatrix_Ecl(T, T_eqx) * SunPos(T);

    // Heliocentric ecliptic coordinates of the comet, equinox T_eqx
    Kepler (GM_Sun, Mjd0, Date, q, e, PQR, r_helioc, v_helioc);
    
    // Geometric geocentric coordinates of the comet
    r_geoc = r_helioc + R_Sun;

    // First-order light-time correction
    dist = Norm(r_geoc);
    fac  = 0.00578*dist;
    r_geoc = r_geoc - fac*v_helioc;

    // Equatorial coordinates
    r_equ = Ecl2EquMatrix(T_eqx) * r_geoc;

    // Output
    cout << DateTime(Date,HHh) 
         << fixed << setprecision(1)
         << setw(7) << Deg*R_Sun[phi] 
         << setw(7) << Deg*r_helioc[phi]
         << setw(6) << Deg*r_helioc[theta]
         << setprecision(3) << setw(7)  << r_helioc[r]
         << setprecision(1) << setw(12) << Angle(Deg*r_equ[phi]/15.0,DMMSSs)
         << "  " << showpos << setw(9)  << Angle(Deg*r_equ[theta],DMMSS)
         << noshowpos << setprecision(6) << setw(11) << dist
         << endl;
    
    ++n_line; 
    
    if ( (n_line % 5) == 0 ) cout << endl; // insert line feed every 5 lines

    Date += Step; // Next time step

  };  // End of time loop

  if (OutFile.is_open()) OutFile.close();
}
