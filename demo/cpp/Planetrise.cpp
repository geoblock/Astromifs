//-------------------------------------------------------------------------
//
// Planetrise program (Planetrise.cpp)
//
// Calculates rising and setting times of the major planets
//             and the Sun
//
//-------------------------------------------------------------------------

#include <cmath>
#include <iomanip>
#include <iostream>

#include "APC_Const.h"
#include "APC_Math.h"
#include "APC_Planets.h"
#include "APC_PrecNut.h"
#include "APC_Spheric.h"
#include "APC_Time.h"
#include "APC_VecMat3D.h"

#ifdef __GNUC__  // GNU C++ adaptation
#include "GNU_iomanip.h"
#endif

using namespace std;


//------------------------------------------------------------------------------
//
// GetInput: Prompts user for input
//
// Output:
//
//   MJD       Time as Modified Julian Date
//   Lambda    Geographical longitude of the observer (east pos.)in [rad]
//   Phi       Geographical latitude of the observer in [rad]
//   Zone      Difference local time - universal time in [d]
//
//------------------------------------------------------------------------------
void GetInput ( double& MJD, double& lambda,
				double& phi,  double& zone )
{
  //
  // Variables
  //
  int    year, month, day;


  cout << " Date (yyyy mm dd)                     ... ";
  cin >> year >> month >> day; cin.ignore(81,'\n');

  cout << endl;
  cout << " Observing site:  east longitude [deg] ... ";
  cin >> lambda; cin.ignore(81,'\n');
  cout << "                  latitude [deg]       ... ";
  cin >> phi; cin.ignore(81,'\n');
  cout << "                  local time - UT [h]  ... ";
  cin >> zone; cin.ignore(81,'\n');

  lambda*= Rad;
  phi*=Rad;
  zone/=24.0;

  MJD = Mjd(year, month, day) - zone; // MJD at 0h LT
}



//------------------------------------------------------------------------------
//
// GeocPosition: Geocentric equatorial coordinates of the Sun and planets 
//
// Input:
//
//   Planet    Identifies planet (or Sun)
//   T         Time in Julian centuries since J2000
//
// Output:
//
//   Ra        Right ascension in [rad]
//   Dec       Declination in [rad]
//
//------------------------------------------------------------------------------
void GeocPosition ( PlanetType Planet, double T,
                    double& Ra, double& Dec )
{
  //
  // Variables
  //
  Vec3D r_plan = Ecl2EquMatrix(T) * 
               ( PertPosition(Planet,T) - PertPosition(Earth,T) );


  Ra  = r_plan[phi];
  Dec = r_plan[theta];
}


//-------------------------------------------------------------------------
//
// Main program
//
//-------------------------------------------------------------------------
int main()
{
  //
  // Types
  //
  enum enEvent { Rising, Transit, Setting };
  enum enState { Never_Rises, Ok, Circumpolar };

  
  //
  // Constants
  //
  const double Sid     =  0.9972696;   // Conversion sidereal/solar time 
  const double Sin_h0p = -9.890038e-3; // sin(-34'); altitude for planets
  const double Sin_h0s = -1.45439e-2;  // sin(-50'); altitude for the Sun
  const int    MaxIter =  10;          // Max. number of iteration steps

  const char* Name[] = 
   { "Sun", "Mercury", "Venus", "Earth", "Mars",
     "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto" };

  
  //
  // Variables
  //
  PlanetType  Planet;
  enState     State;
  enEvent     Event;
  int         iPlanet, iEvent, cnt;
  double      lambda, zone, phi, Cphi, Sphi;
  double      Date, T, LT, D_LT;
  double      LST_0h, LST;
  double      Ra, RA_0h, RA_24h, Dec, Dec_0h, Dec_24h;
  double      Sin_h0, SDA, Dtau;

  
  // Title
  cout << endl
       << "  PLANRISE: planetary and solar rising and setting times " << endl
       << "            (c) 1999 Oliver Montenbruck, Thomas Pfleger  " << endl
       << endl;


  GetInput (Date, lambda, phi, zone );
  
  Cphi = cos(phi);
  Sphi = sin(phi);

  
  // Local sidereal time at 0h local time 
  LST_0h = GMST(Date) + lambda;


  // Header
  cout << endl
       << setw(18) << "rise" 
       << setw(15) << "culmination"
       << setw(7)  << "set"
       << endl << endl;


  // Loop over bodies
  for (iPlanet=Sun; iPlanet<=Pluto; iPlanet++) {
    
    Planet = (PlanetType) iPlanet;

    if (Planet==Earth) continue;  // Skip Earth

    // Compute geocentr. planetary position at 0h and 24h local time
    T = (Date-MJD_J2000)/36525.0;
    
    GeocPosition(Planet, T            , RA_0h , Dec_0h );
    GeocPosition(Planet, T+1.0/36525.0, RA_24h, Dec_24h);

    // Generate continuous right ascension values in case of jumps
    // between 0h and 24h                                          
    if ( RA_0h-RA_24h > pi ) RA_24h = RA_24h + 2.0*pi;
    if ( RA_0h-RA_24h <-pi ) RA_0h  = RA_0h  + 2.0*pi;

    // Print name
    cout << " " << setw(10) << left << Name[Planet] << right << "   ";
  
    // Compute rising, transit and setting times
    State = Ok;

    for (iEvent=Rising; iEvent<=Setting; iEvent++) {

      Event = (enEvent) iEvent;
      
      // Starting value 12h local time 
      LT  = 12.0;  
      cnt = 0;

      
	  // Iteration
      do {

        // Linear interpolation of planetary position 
        Ra  = RA_0h  + (LT/24.0) * (RA_24h  - RA_0h );
        Dec = Dec_0h + (LT/24.0) * (Dec_24h - Dec_0h);

        // Compute semi-diurnal arc (in rad)
        if (Planet==Sun) Sin_h0=Sin_h0s; else Sin_h0=Sin_h0p;

        SDA = (Sin_h0 - sin(Dec)*Sphi ) / (cos(Dec)*Cphi);

        if (fabs(SDA)<1.0) {
           SDA = acos(SDA); State=Ok; 
        }
        else {
         
          // Test for circumpolar motion or invisibility
          if (phi*Dec>=0.0)  
            State = Circumpolar; 
          else 
            State = Never_Rises;
        
          break;  // Terminate iteration loop
        }
        
        // Local sidereal time 
        LST = LST_0h + LT/(Sid*12.0/pi); 

        // Difference in hour angle
        switch (Event) {
          case Rising : Dtau = (LST-Ra) + SDA; break;
          case Transit: Dtau = (LST-Ra);       break;
          case Setting: Dtau = (LST-Ra) - SDA;
        }

        // Improved times for rising, culmination and setting
        D_LT = Sid*(12.0/pi)*(Modulo(Dtau+pi,2*pi)-pi);
        LT   = LT-D_LT;
        cnt++;

      }
      while ( (fabs(D_LT)>0.008) && (cnt<=MaxIter) );

      // Print result 
      switch (State) {
        case Ok          : cout << Time(LT,HHMM) << setw(6) << " ";      break;
        case Never_Rises : cout << "-------- always invisible --------"; break;
        case Circumpolar : cout << "--------- always visible ---------";
      }
    
    } // End Event loop
    
    cout << endl;
  
  } // End Planet loop

  
  // Trailer
  cout << endl;

  cout << " all times in local time (=UT"
       << showpos << 24.0*zone << noshowpos << "h)" 
       << endl;
}
