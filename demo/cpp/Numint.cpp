//------------------------------------------------------------------------------
//
// Numint.cpp
//
// Numerical integration of perturbed minor planet orbits
//
//------------------------------------------------------------------------------

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "APC_Const.h"
#include "APC_DE.h"
#include "APC_IO.h"
#include "APC_Kepler.h"
#include "APC_Math.h"
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
// Accel: Calculates the acceleration vector for a minor body in 
//        the solar system
//
// Input:
//
//   Mjd      Time as Modified Julian Date
//   r        J2000 ecliptical coordinates of minor body
//
// <return>:  J2000 ecliptical acceleration in [AU/d**2]
//
//------------------------------------------------------------------------------
Vec3D Accel (double Mjd, const Vec3D& r)
{
  //
  // Constants
  //
  
  // Grav. constant * solar and planetary masses in AU**3/d**2
  static const double  GM[10] = { 
                           GM_Sun,                 // Sun
                           GM_Sun / 6023600.0,     // Mercury
                           GM_Sun /  408523.5,     // Venus
                           GM_Sun /  328900.5,     // Earth
                           GM_Sun / 3098710.0,     // Mars
                           GM_Sun /    1047.355,   // Jupiter
                           GM_Sun /    3498.5,     // Saturn
                           GM_Sun /   22869.0,     // Uranus
                           GM_Sun /   19314.0,     // Neptune
                           GM_Sun / 3000000.0      // Pluto
                         };

  //
  // Variables
  //
  int         iPlanet;
  PlanetType  Planet;
  Vec3D       a, r_p, d;
  Mat3D       P;
  double      T, D;

  
  // Solar attraction
  D = Norm(r);
  a = - GM_Sun * r / (D*D*D);

  
  // Planetary perturbation
  T = ( Mjd - MJD_J2000 ) / 36525.0;
  P = PrecMatrix_Ecl ( T, T_J2000 );    // Precession 

  for (iPlanet=Mercury; iPlanet<=Neptune; iPlanet++) {

    Planet = (PlanetType) iPlanet;

    // Planetary position (ecliptic and equinox of J2000)
//  r_p = P * KepPosition(Planet,T);    // Fast but low precison
    r_p = P * PertPosition(Planet,T);   // Accurate but slow
    d   = r - r_p;

    // Direct acceleration
    D = Norm(d);
    a += - GM[iPlanet] * d / (D*D*D);
  
    // Indirect acceleration
    D = Norm(r_p);
    a += - GM[iPlanet] * r_p / (D*D*D);
  }

  return a;
}


//------------------------------------------------------------------------------
//
// F: Equation of motion in a form compatible with the SolverDE integrator
//
// Input:
//
//   X        Independent variable (Time as Modified Julian Date)
//   Y[]      State vector (position and velocity)
//
// Output: 
//
//   dYdX[]   Derivative of state vector w.r.t. X (veloc. and acceleration)
//
//------------------------------------------------------------------------------
void F (double X, double Y[], double dYdX[])
{
  //
  // Variables
  //
  Vec3D a, r;


  r = Vec3D(Y[1], Y[2], Y[3]);
  a = Accel(X, r);

  dYdX[0] = 0.0;
  dYdX[1] = Y[4]; dYdX[2] = Y[5];  dYdX[3] = Y[6];  // velocity
  dYdX[4] = a[x]; dYdX[5] = a[y];  dYdX[6] = a[z];  // acceleration
}


//------------------------------------------------------------------------------
//
// WriteElm: Screen output of orbital elements
//
// Input:
//
//   MjdEpoch Epoch as Modified Julian Date
//   a        Semimayor axis of the orbit in [AU]
//   e        Eccentricity of the orbit
//   i        Inclination of the orbit to the ecliptic in [rad]
//   Omega    Longitude of the ascending node of the orbit in [rad]
//   omega    Argument of perihelion in [rad]
//   M        Mean anomaly in [rad]
//   T_eqx0   Equinox in Julian centuries since J2000
//
//------------------------------------------------------------------------------
void WriteElm ( double MjdEpoch, double a, double e, double i,
                double Omega, double omega, double M, double T_eqx0 )
{
  cout << "  Epoch (y m d)           "
       << "  " << setprecision(2) << DateTime(MjdEpoch,DDd) << endl;
  cout << "  Semi-major axis (a)     "
       << setprecision(7) << setw(14) << a << " AU" << endl;
  cout << "  Eccentricity (e)        "
       << setprecision(7) << setw(14) << e << endl;
  cout << "  Inclination (i)         "
       << setprecision(5) << setw(12) << i*Deg << " deg" << endl;
  cout << "  Long. of ascending node "
       << setprecision(5) << setw(12) << Omega*Deg << " deg" << endl;
  cout << "  Argument of perihelion  "
       << setprecision(5) << setw(12) << omega*Deg << " deg" << endl; 
  cout << "  Mean anomaly (M)        "
       << setprecision(5) << setw(12) << M*Deg     << " deg" << endl; 
  cout << "  Equinox                 "
       << setprecision(2) << setw(9)  << 2000.0+100*T_eqx0 << endl;
  cout << endl;
};


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
//   MjdEpoch Epoch as Modified Julian Date
//   a        Semimayor axis of the orbit [AU]
//   e        Eccentricity of the orbit
//   M        Mean anomaly in [rad]
//   PQR      Transformation matrix orbital plane -> ecliptic
//   T_eqx0   Equinox for the orbital elements
//
//------------------------------------------------------------------------------

void GetElm ( char* Filename,
              double& MjdEpoch,
              double& a, double& e, double& M, Mat3D& PQR, double& T_eqx0 )
{
  //
  // Variables
  //
  int       year, month;
  double    Omega, i, omega;
  double    Year, day;
  ifstream  inp;


  // Read elements
  inp.open(Filename);

  inp >> year >> month >> day; inp.ignore(81,'\n');
  inp >> a; inp.ignore(81,'\n');
  inp >> e; inp.ignore(81,'\n');
  inp >> i; inp.ignore(81,'\n');
  inp >> Omega; inp.ignore(81,'\n');
  inp >> omega; inp.ignore(81,'\n');
  inp >> M; inp.ignore(81,'\n');
  inp >> Year; inp.ignore(81,'\n');
  
  inp.close();
  

  // Compute derived data
  i*=Rad;
  Omega*=Rad;
  omega*=Rad;
  M*=Rad;

  MjdEpoch = Mjd(year,month,int(day))+(day-int(day));
  T_eqx0   = (Year-2000.0)/100.0;
  PQR      = GaussVec(Omega,i,omega);


  // Print elements
  cout << " Orbital elements from file " << Filename << ":" << endl << endl; 
  WriteElm ( MjdEpoch, a,e,i,Omega,omega,M, T_eqx0 );
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
//   T_eqx     Desired equinox for ephemeris in jul. cent. since J2000
//
//------------------------------------------------------------------------------
void GetEph (double& MjdStart, double& Step, double& MjdEnd, double& T_eqx)
{
  //
  // Variables
  //
  int    year, month, day;
  double hour, Year;


  cout << " Begin and end of the ephemeris: " << endl;

  cout << "  First date (yyyy mm dd hh.hhh)           ... ";
  cin >> year >> month >> day >> hour; cin.ignore(81,'\n');
  MjdStart = Mjd(year,month,day) + hour/24.0 ;

  cout << "  Final date (yyyy mm dd hh.hhh)           ... ";
  cin >> year >> month >> day >> hour; cin.ignore(81,'\n');
  MjdEnd = Mjd(year,month,day) + hour/24.0 ;

  cout << "  Step size (dd hh.hh)                     ... ";
  cin >> day >> hour; cin.ignore(81,'\n');
  Step  = day + hour/24.0;

  cout << endl
       << " Desired equinox of the ephemeris (yyyy.y) ... ";
  cin >> Year; cin.ignore(81,'\n');
  T_eqx = (Year-2000.0) /100.0;
}


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  //
  // Constants
  //
  const pol_index R    = pol_index(2);   // Index for vector length
  const int       Neqn = 6;              // Number of differential eqns.
  const double    eps  = 1.0e-10;        // Relative accuracy


  //
  // Variables
  //
  int       n_line;
  double    MjdStart, Step, MjdEnd, T_eqx;
  double    MjdEpoch, a, e, M, T_eqx0;
  double    Omega, i, omega;
  Mat3D     PQR, P;
  double    Mjd, MjdStep, T;
  SolverDE  Orbit(F, Neqn);
  double    Y[Neqn+1];
  double    relerr,abserr;
  DE_STATE  State;
  Vec3D     r, v;
  Vec3D     R_Sun, r_helioc, v_helioc, r_geoc, r_equ;
  double    dist, fac;
  char      InputFile[APC_MaxFilename] = "";
  char      OutputFile[APC_MaxFilename] = "";
  bool      FoundInputfile = false;
  bool      FoundOutputfile = false;
  ofstream  OutFile;


  // Title
  cout << endl
       << "  NUMINT: numerical integration of perturbed minor planet orbits" 
       << endl
       << "           (c) 1999 Oliver Montenbruck, Thomas Pfleger          " 
       << endl << endl;

  
  // Find input and optional output files
  GetFilenames( argc, argv, "Numint.dat", InputFile, FoundInputfile,
                OutputFile, FoundOutputfile );
  
  // Terminate program if input file could not be found
  if (!FoundInputfile) {
    cerr << " Terminating program." << endl;
    exit(-1);
  }

  
  // Read orbital elements from file and prompt user for ephemeris data
  GetElm ( InputFile, MjdEpoch, a, e, M, PQR, T_eqx0 );
  GetEph ( MjdStart, Step, MjdEnd, T_eqx );
  
  
  // Redirect output if output file shall be created
  if (FoundOutputfile) {
    OutFile.open(OutputFile);
    if (OutFile.is_open())
      cout = OutFile;
  }
  
  
  P = PrecMatrix_Ecl(T_J2000,T_eqx); 

  
  // Initial state vector (ecliptic and equinox of J2000)
  PQR = PrecMatrix_Ecl(T_eqx0,T_J2000) * PQR; 

  Ellip ( GM_Sun, M,a,e, r,v );
  
  r = PQR*r;
  v = PQR*v;

  Y[0]=0.0;  // Unused
  Y[1]=r[x]; Y[2]=r[y]; Y[3]=r[z];
  Y[4]=v[x]; Y[5]=v[y]; Y[6]=v[z];

  Mjd = MjdEpoch;
  
  
  // Start integration: propagate state vector from epoch 
  // to start of ephemeris                                
  State = DE_INIT; 
  relerr = eps;  abserr = 0.0;
    
  do {                                    
    
    Orbit.Integ(Y, Mjd, MjdStart, relerr, abserr, State);
    
    if ( State==DE_INVALID_PARAMS ) { 
      cerr << "Exit (invalid parameters)" << endl;
      exit(1); 
    }
  }
  while ( State > DE_DONE );


  // Orbital elements at start of ephemeris (equinox T_eqx)
  r = P*Vec3D(Y[1],Y[2],Y[3]); 
  v = P*Vec3D(Y[4],Y[5],Y[6]);

  Elements (GM_Sun, r, v, a, e, i, Omega, omega, M);


  // Print elements
  cout << endl << endl
       << " Orbital elements at start epoch:" << endl
       << endl;
  WriteElm ( Mjd, a,e,i,Omega,omega,M, T_eqx );


  // Ephemeris header
  cout << endl << endl
       << "    Date    ET    Sun     l      b     r"
       << "       Ra          Dec      Distance " << endl
       << setw(43) << " " << "   h  m  s      o  '  \"    (AU) " << endl;



  // Compute ephemeris
  n_line = 0;
  MjdStep = MjdStart;

  
  // Time loop
  while ( MjdStep < MjdEnd + Step/2 ) {

    // Integrate orbit to MjdStep
    do {                                    
      
      Orbit.Integ(Y, Mjd, MjdStep, relerr, abserr, State);
      
      if (State==DE_INVALID_PARAMS) { 
        cerr << "Exit (invalid parameters)" << endl;
        exit(1);
      }
    }
    while ( State > DE_DONE );


    // Heliocentric ecliptic coordinates, equinox T_eqx
    r_helioc = P*Vec3D(Y[1],Y[2],Y[3]); 
    v_helioc = P*Vec3D(Y[4],Y[5],Y[6]);

    
    // Geocentric ecliptic coordinates of the Sun, equinox T_eqx
    T = ( Mjd - MJD_J2000 ) / 36525.0;
    R_Sun = PrecMatrix_Ecl(T,T_eqx) * SunPos(T);


    // Geometric geocentric coordinates
    r_geoc = r_helioc + R_Sun;

    
    // First-order light-time correction
    dist = Norm(r_geoc);
    fac  = 0.00578*dist;
    r_geoc = r_geoc - fac*v_helioc;

    
    // Equatorial coordinates
    r_equ = Ecl2EquMatrix(T_eqx) * r_geoc;

    
    // Output
    cout << DateTime(Mjd,HHh) 
         << fixed << setprecision(1)
         << setw(7) << Deg*R_Sun[phi] 
         << setw(7) << Deg*r_helioc[phi]
         << setw(6) << Deg*r_helioc[theta]
         << setprecision(3) << setw(7) << r_helioc[R]
         << setprecision(1) << setw(12) << Angle(Deg*r_equ[phi]/15.0,DMMSSs)
         << "  " << showpos << setw(9) << Angle(Deg*r_equ[theta],DMMSS)
         << noshowpos << setprecision(6) << setw(11) << dist
         << endl;
    
    ++n_line; 
    
    if ( (n_line % 5) ==0 ) cout << endl; // insert line feed every 5 lines

    MjdStep += Step;  // Next time step

  } // End time loop

  if (OutFile.is_open()) OutFile.close();
}
