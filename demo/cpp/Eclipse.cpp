//------------------------------------------------------------------------------
//
// Eclipse program (Eclipse.cpp)
//
// Computes the center line of a solar eclipse
//
//------------------------------------------------------------------------------

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "APC_Cheb.h"
#include "APC_Const.h"
#include "APC_IO.h"
#include "APC_Math.h"
#include "APC_Moon.h"
#include "APC_PrecNut.h"
#include "APC_Spheric.h"
#include "APC_Sun.h"
#include "APC_Time.h"
#include "APC_VecMat3D.h"

#ifdef __GNUC__  // GNU C++ Adaptations
#include <memory>
#include "GNU_iomanip.h"
#endif

using namespace std;


//
// Types
//
enum enPhase { NoEclipse, partial, NonCenAnn, NonCenTot, annular, total };


//
// Constants
//
const double Interval = 1.0/36525.0;   // 1d [cy]


//
// Global lunar and solar ephemeris data
//
Cheb3D ChebMoonEqu(MoonEqu,10,Interval);
Cheb3D ChebSunEqu (SunEqu,  5,Interval);


//------------------------------------------------------------------------------
//
// Intersect: Computes intersection point of shadow's axis and the Earth,
//            the shadow axis unit vector, umbral diameter and phase
//
// Input:
//
//   T_ET      Ephemeris time in Julian centuries since J2000
//
// Output:
//
//   r         Intersection point of shadow's axis and surface of the Earth;
//             geocentric coordinates in [km]
//   e         Shadow axis unit vector
//   D_umbra   Shadow diameter on the fundamental plane in [km]
//   Phase     Phase of the eclipse
//
//------------------------------------------------------------------------------
void Intersect ( double T_ET,
                 Vec3D& r, Vec3D& e, double& D_umbra, enPhase& Phase )
{
  //
  // Constants
  //
  const double tau_Sun  = 8.32 / ( 1440.0*36525.0);    // 8.32 min  [cy]
  const double fac      = 0.996647;     // Ratio polar/equat. Earth radius

  
  //
  // Variables
  //
  double r_MS,s0,Delta,r0,s;
  double D_penumbra;
  Vec3D  r_Moon; 
  Vec3D  r_Sun;


  // Solar and lunar position w.r.t. true equator and equinox of date
  // (solar coordinates corrected for light time)
  r_Moon = ChebMoonEqu.Value(T_ET);             // [km]
  r_Sun  = ChebSunEqu.Value (T_ET-tau_Sun)*AU;  // [km]
  
  // Scale z-coordinate to compensate Earth flattening
  r_Moon = Vec3D ( r_Moon[x], r_Moon[y], r_Moon[z]/fac );
  r_Sun  = Vec3D ( r_Sun[x] , r_Sun[y] , r_Sun[z]/fac  );

  // Shadow axis (Sun->Moon unit vector)
  r_MS = Norm(r_Moon-r_Sun);
  e    = (r_Moon-r_Sun)/r_MS;

  // Distance of the Moon from the fundamental plane
  s0 = -Dot(r_Moon,e);

  // Distance of the shadow axis from the centre of the Earth
  Delta = s0*s0 + R_Earth*R_Earth - Dot(r_Moon,r_Moon);
  r0 = sqrt(R_Earth*R_Earth-Delta);

  // Umbra and penumbra diameter on the fundamental plane
  D_umbra    = 2.0*((R_Sun-R_Moon)*(s0/r_MS)-R_Moon);
  D_penumbra = 2.0*((R_Sun-R_Moon)*(s0/r_MS)+R_Moon);
  
  // Determine phase and shadow coordinates if required
  if (r0<R_Earth) {

    // Shadow axis intersects the Earth; total or annular eclipse
     
    // Intersection of the shadow axis and the surface of the Earth
    s = s0-sqrt(Delta);
    r = r_Moon + s*e;                
    
    // Re-scale z-component
    r = Vec3D(r[x],r[y],fac*r[z]);   
    
    // Umbra diameter at the surface of the Earth
    D_umbra = 2.0*((R_Sun-R_Moon)*(s/r_MS)-R_Moon);
    
    Phase = ( (D_umbra>0.0)? annular : total );  
  }
  else {

    if ( r0 < R_Earth+0.5*fabs(D_umbra) ) {
      // Non-central eclipse
      Phase = ( (D_umbra>0.0)? NonCenAnn : NonCenTot );  
    }
    else {
      // Partial or no eclipse
      Phase = ( (r0<R_Earth+0.5*fabs(D_penumbra))? partial : NoEclipse );
    };
    
    r = Vec3D(); // Dummy value (Null vector)
  };

};


//------------------------------------------------------------------------------
//
// Central: Computes centre-line coordinates, duration and phase
//
// Input:
//
//   MjdUT     Universal time as Modified Julian Date 
//   ET_UT     Difference of ephemeris time and universal time in [s]
//
// Output:
//
//   Lambda    Geographic longitude (east pos.) of centre of shadow [rad]
//   Phi       Geographic latitude of centre of shadow [rad]
//   Durat     Duration of central phase at location (Lambda, Phi) in [min]
//   Phase     Phase of the eclipse
//
//------------------------------------------------------------------------------
void Central ( double MjdUT, double ET_UT,
               double &Lambda, double& Phi, double& Durat, enPhase& Phase )
{
  //
  // Constants
  //
  const double dt    = 0.1  / (1440.0*36525.0); // 0.1 min [cy]
  const double omega = pi2 * 1.002738*36525.0;  // [rad/cy]
  

  //
  // Variables
  //
  double    T_ET;             // Ephemeris time
  Vec3D     r,r_G,rr,dr;      // Shadow coordinates
  Vec3D     e,ee;             // Shadow axis
  double    D_umbra, DU;      // Umbra diameter
  enPhase   Ph=NoEclipse;     // Phase 
  double    w;                // Earth rotation angle                                       
  double    drp;              // Shadow displacement

  
  // Ephemeris Time
  T_ET = ( MjdUT + ET_UT/86400.0 - MJD_J2000 ) / 36525.0;

  // Intersection of shadow axis with the surface of the Earth
  Intersect ( T_ET, r, e, D_umbra, Phase ); 

  // For central phase only: geogr. coordinates and duration of central phase
  Lambda = Phi = Durat = 0.0;
  
  if ( Phase >= annular ) {

    // Geographic shadow coordinates
    r_G    = R_z(GMST(MjdUT))*r;               // Greenwich coordinates
    Lambda = Modulo(r_G[phi]+pi,2*pi)-pi;      // East longitude
    Phi    = r_G[theta];                       // Geocentric latitude
    Phi    = Phi + 0.1924*Rad*sin(2*Phi);      // Geographic latitude
    
    // Duration of central phase for this place
   
    // (a) shadow coordinates at time T+dt (or T-dt)
    if (Ph<annular) {
      Intersect ( T_ET-dt, rr, ee, DU, Ph ); w = -dt*omega;
    }
    else {
      Intersect ( T_ET+dt, rr, ee, DU, Ph ); w = +dt*omega;
    }

    // (b) displacement dr of the shadow on Earth and 
    //     fraction drp perpendicular to the shadow axis
    dr  = Vec3D ( rr[x]-r[x]+w*r[y], rr[y]-r[y]-w*r[x], rr[z]-r[z] );
    drp = Norm ( dr - Dot(dr,e)*e );
    
    Durat = (36525.0*1440.0)*dt * fabs(D_umbra) / drp;  // [min]
  };  

};


//------------------------------------------------------------------------------
//
// GetInput: Prompts user for approximate date and time of New Moon,
//           desired step size and ET-UT
//
// Output:
//
//   MjdStart  Start time for calculation of eclipse
//   Step      Step size in [d]
//   MjdEnd    End time for calculation of eclipse
//   ET_UT     Difference ephemeris time - universal time in [s]
//
// Note: MjdStart and MjdEnd as Modified Julian Dates
//
//------------------------------------------------------------------------------
void GetInput(double& MjdStart, double& Step, double& MjdEnd, double& ET_UT)
{
  //
  // Variables
  //
  int    year, month, day, step;
  double hour, MJD;
  bool   valid;


  // Date and step size
  cout << endl
       << " Date of New Moon (yyyy mm dd hh.h)     ... ";
  cin >> year >> month >> day >> hour; cin.ignore(81,'\n');

  cout << " Output step size (min)" << setw(21) << right << "... ";
  cin >> step; cin.ignore(81,'\n');


  // Round date to step size
  Step = step/1440.0;  // [d]

  MJD = Mjd(year,month,day) + floor((hour/24.0)/Step+0.5)*Step;


  // Start and end date (date of New Moon +/- 6 hours)
  MjdStart = MJD - 0.25;
  MjdEnd   = MJD + 0.25;


  // Difference between ephemeris time and universal time
  ETminUT ( (MJD-MJD_J2000)/36525.0, ET_UT, valid );
  if ( valid ) {
    cout << " Difference ET-UT (proposal:" << fixed << setw(6) 
         << setprecision (1) << ET_UT << " sec) ... ";
    cin >> ET_UT; cin.ignore(81,'\n');
  }
  else {
    cout << " Difference ET-UT (sec)" << setw(21) << right << "... "; 
    cin >> ET_UT; cin.ignore(81,'\n');
  }
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
  double    MjdStart, Step, MjdEnd;
  double    MjdUT, ET_UT;
  double    Lambda, Phi, Durat;
  enPhase   Phase;
  char      OutputFile[APC_MaxFilename] = "";
  ofstream  OutFile;
    

  // Title
  cout << endl
       << "     ECLIPSE: central line and duration of solar eclipses " << endl
       << "          (c) 1999 Oliver Montenbruck, Thomas Pfleger     " << endl
       << endl;

  
  // Prompt user for input
  GetInput ( MjdStart, Step, MjdEnd, ET_UT );


  // Check for output file and eventually redirect output
  if (argc==2) {
    strncpy(OutputFile, argv[1], APC_MaxFilename);
    OutFile.open(OutputFile);
    if (OutFile.is_open())
      cout = OutFile;
  }


  // Header
  cout << endl << endl
       << "    Date      UT         Phi     Lambda   Durat  Phase " << endl
       << "              h  m        o  '      o  '   min         " << endl
       << endl;

  // Phase and central line of the eclipse
  MjdUT = MjdStart;

  
  // Time loop
  while ( MjdUT < MjdEnd + Step/2 ) {
  
    // Phase and location of shadow point     
    Central ( MjdUT, ET_UT, Lambda, Phi, Durat, Phase );

    // Output
    if ( Phase != NoEclipse ) {

      cout << " " << DateTime(MjdUT,HHMM);

      if ( Phase < annular )
        cout << "       -- --     -- --    ---";
      else {
        cout << showpos
             << "      " << setw(6) << Angle(Deg*Phi,DMM)
             << "   "    << setw(7) << Angle(Deg*Lambda,DMM)
             << noshowpos 
             << fixed << setprecision(1) << setw(7) << Durat;
      };
      
      switch(Phase) {
      
        case partial  : cout << "   partial               "; break;
        case NonCenAnn: cout << "   annular (non-central) "; break;
        case NonCenTot: cout << "   total (non-central)   "; break;
        case annular  : cout << "   annular               "; break;
        case total    : cout << "   total                 ";
      };
      
      cout << endl;
    };

    MjdUT += Step;  // Next time step

  };  // End time loop

  if (OutFile.is_open()) OutFile.close();
}
