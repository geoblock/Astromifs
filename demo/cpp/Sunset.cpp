//------------------------------------------------------------------------------
//
// Sunset.cpp
//
// Rising and setting times of Sun and Moon and twilight times
//
//------------------------------------------------------------------------------

#include <cmath>
#include <iomanip>
#include <iostream>

#include "APC_Const.h"
#include "APC_Math.h"
#include "APC_Moon.h"
#include "APC_Spheric.h"
#include "APC_Sun.h"
#include "APC_Time.h"

#ifdef __GNUC__  // GNU C++ adaptation
#include <ctype.h>
#include "GNU_iomanip.h"
#endif

using namespace std;


//
// Types
// 

// Events to search for
enum enEvent { 
  Moon,      // indicates moonrise or moonset
  Sun,       // indicates sunrise or sunset
  CivilTwi,  // indicates Civil twilight
  NautiTwi,  // indicates nautical twilight
  AstroTwi   // indicates astronomical twilight
};


//------------------------------------------------------------------------------
//
// SinAlt: Sine of the altitude of Sun or Moon
//
// Input:
//
//   Event     Indicates event to find
//   MJD0      0h at date to investigate (as Modified Julian Date)
//   Hour      Hour
//   lambda    Geographic east longitude in [rad]
//   Cphi      Cosine of geographic latitude
//   Sphi      Sine of geographic latitude
//
// <return>:   Sine of the altitude of Sun or Moon at instant of Event
//
//------------------------------------------------------------------------------
double SinAlt ( enEvent Event, double MJD0, double Hour,
                double lambda, double Cphi, double Sphi )
{
  //
  // Variables
  //  
  double  MJD, T, Ra, Dec, tau;
  

  MJD = MJD0 + Hour/24.0;         
  T   = (MJD-51544.5)/36525.0;

  if ( Event == Moon )
    MiniMoon ( T, Ra, Dec);
  else
    MiniSun ( T, Ra, Dec );

  tau = GMST(MJD) + lambda - Ra;
  
  return ( Sphi*sin(Dec)+Cphi*cos(Dec)*cos(tau) );
}


//------------------------------------------------------------------------------
//
// GetInput: Prompts user for input
//
// Output:
//
//   MJD       Time as Modified Julian Date
//   lambda    Geographic east longitude of the observer in [rad]
//   phi       Geographic latitude of the observer in [rad]
//   zone      Difference local time - universal time in [d]
//   twilight  Indicates civil, nautical or astronomical twilight
//
//------------------------------------------------------------------------------
void GetInput ( double& MJD, double& lambda,
                double& phi, double& zone, enEvent& twilight )
{
  //
  // Variables
  //  
  int  year, month, day;
  char cTwilight;


  cout << " First date (yyyy mm dd)               ... ";
  cin >> year >> month >> day; cin.ignore(81,'\n');

  cout << endl;
  cout << " Observing site:  east longitude [deg] ... ";
  cin  >> lambda; cin.ignore(81,'\n');
  cout << "                  latitude [deg]       ... ";
  cin  >> phi; cin.ignore(81,'\n');
  cout << "                  local time - UT [h]  ... ";
  cin  >> zone; cin.ignore(81,'\n');
  cout << " Twilight definition (c, n or a)       ... ";
  cin  >> cTwilight; cin.ignore(81,'\n');
  
  switch (tolower(cTwilight)) {
    case 'c': twilight = CivilTwi; break;
    case 'a': twilight = AstroTwi; break;
    case 'n': 
    default:  twilight = NautiTwi;
  }

  lambda*=Rad;
  phi*=Rad;
  zone/=24.0;

  MJD = Mjd(year,month,day) - zone;
}


//------------------------------------------------------------------------------
//
// FindEvents: Search for rise/set/twilight events of Sun or Moon
//
// Input:
//
//   Event     Indicates event to search for
//   MJD0h     0h at desired date as Modified Julian Date
//   lambda    Geographic east longitude of the observer in [rad]
//   phi       Geographic latitude of the observer in [rad]
//
// Output:
//
//   LT_Rise   Local time of rising or beginning of twilight
//   LT_Set    Local time of setting or end of twilight
//   rises     Event takes place
//   sets      Event takes place
//   above     Sun or Moon is circumpolar
//
//------------------------------------------------------------------------------
void FindEvents ( enEvent Event, double MJD0h, double lambda, double phi,              
                  double& LT_Rise, double& LT_Set, 
                  bool& rises, bool& sets, bool& above )
{
  //
  // Constants
  //
  static const double sinh0[5] = {

    sin(Rad*( +8.0/60.0)), // Moonrise              at h= +8' 
    sin(Rad*(-50.0/60.0)), // Sunrise               at h=-50'
    sin(Rad*(   - 6.0  )), // Civil twilight        at h=-6 deg
    sin(Rad*(   -12.0  )), // Nautical twilight     at h=-12deg
    sin(Rad*(   -18.0  )), // Astronomical twilight at h=-18deg
  };
  
  const double Cphi = cos(phi);
  const double Sphi = sin(phi);

  
  //
  // Variables
  //
  double hour = 1.0;
  double y_minus, y_0, y_plus;
  double xe, ye, root1, root2;
  int nRoot;
  
  
  // Initialize for search
  y_minus = SinAlt(Event, MJD0h, hour-1.0, lambda, Cphi, Sphi) - sinh0[Event];

  above = (y_minus>0.0);
  rises = false;
  sets  = false;

  
  // loop over search intervals from [0h-2h] to [22h-24h]
  do {

    y_0    = SinAlt( Event, MJD0h, hour    , lambda, Cphi, Sphi )-sinh0[Event];
    y_plus = SinAlt( Event, MJD0h, hour+1.0, lambda, Cphi, Sphi )-sinh0[Event];

    // find parabola through three values y_minus,y_0,y_plus
    Quad ( y_minus, y_0, y_plus, xe, ye, root1, root2, nRoot );

    if ( nRoot==1 ) {
      if ( y_minus < 0.0 ) 
        { LT_Rise = hour+root1;  rises = true; }
      else  
        { LT_Set  = hour+root1;  sets  = true; }
    }                   

    if ( nRoot == 2 ) {
      if ( ye < 0.0 )
        { LT_Rise = hour+root2;  LT_Set = hour+root1; }
      else
        { LT_Rise = hour+root1;  LT_Set = hour+root2; }
      rises = true; sets = true;
    }          
    
    y_minus = y_plus;     // prepare for next interval
    hour += 2.0;         

  }
  while ( !( ( hour == 25.0 ) || ( rises && sets ) ) );
}


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------
int main()
{
  //
  // Variables
  //
  bool     above, rise, sett;
  int      iEvent, day;
  double   lambda, zone, phi;
  double   date, start_date, LT_Rise, LT_Set;
  enEvent  Event, Twilight;

  
  // Title
  cout << endl
       << "      SUNSET: solar and lunar rising and setting times  " << endl
       << "        (c) 1999 Oliver Montenbruck, Thomas Pfleger     " << endl
       << endl;


  // Prompt user for input
  GetInput (start_date, lambda, phi, zone, Twilight );
  

  // Header
  cout << endl
       << "    Date"
       << "             Moon              Sun             Twilight " << endl
       << "        "
       << "           rise/set          rise/set        beginning/end"
       << endl << endl;


  // loop over 10 subsequent days
  for (day=0; day<10; day++) {
    
    // current date
    date = start_date + day; 
    cout << " " << DateTime(date+zone) << "  ";        

    // loop over cases (Moon, Sun, Twilight)
    for (iEvent=0; iEvent<=2; iEvent++)  {

      // After events for Moon and Sun: specify kind of twilight event
      Event = (iEvent<2) ? (enEvent) iEvent : Twilight;

      // Now try to find times of events
      FindEvents ( Event, date, lambda, phi,
                   LT_Rise, LT_Set, rise, sett, above );

      
      // Output
      if ( rise || sett )  {
        if ( rise )
          cout << "   " << Time(LT_Rise,HHMM) << " ";
        else
          cout << "   ----- ";
        if ( sett )
          cout << "   " << Time(LT_Set,HHMM) << " ";
        else
          cout << "   ----- ";
      }
      else 
        if ( above ) {
          if ( Event >= CivilTwi )
            cout << "    always bright ";
          else 
            cout << "   always visible ";
        }
        else {
          if ( Event >= CivilTwi )
            cout << "     always dark  ";
          else 
            cout << "  always invisible";
        }

    }

    cout << endl;
  }

  
  // Trailer
  cout << endl;

  cout << " all times in local time (=UT"
       << showpos << 24.0*zone << noshowpos << "h)" 
       << endl;
}
