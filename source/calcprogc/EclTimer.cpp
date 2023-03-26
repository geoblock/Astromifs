//------------------------------------------------------------------------------
//
// EclTimer program (EclTimer.cpp)
//
// Computes local circumstances of a solar eclipse
//                                                                           
//------------------------------------------------------------------------------

#include <cmath>
#include <iomanip>
#include <iostream>

#include "APC_Cheb.h"
#include "APC_Const.h"
#include "APC_Math.h"
#include "APC_Moon.h"
#include "APC_PrecNut.h"
#include "APC_Spheric.h"
#include "APC_Sun.h"
#include "APC_Time.h"
#include "APC_VecMat3D.h"

#ifdef __GNUC__  // GNU C++ adaptation
#include "GNU_iomanip.h"
#endif

using namespace std;


namespace // unnamed namespace
{
  //
  // Constants
  //
  const int    Degree   = 8;             // Degree of Chebyshev polynomials
  const double Interval = 2.0/36525.0;   // Interval for Chebyshev polyn. [cy]  

  const double k_Moon   = 0.2725076;     // Ratio Moon/Earth radii (penumbral)
  const double k_Sun    = R_Sun/R_Earth; // Ratio Sun/Earth radii
}


//
// Types
//
enum enPhase { noEclipse, partial, annular, total };
enum enShadow { Penumbra, Umbra };

// Structure holding the local circumstances of an eclipse
struct LocCircs {
   
  // Times in Julian cent. since J2000:
  double Times[4];  // Contact times (1st to 4th contact) 
  double TMax;      // Time of maximum eclipse

  double Mag;       // Magnitude of the eclipse
  double Obsc;      // Maximum obscuration
  enPhase Phase;    // Maximal phase
};


//
// Global variables and objects
//
enShadow ShadowType;

double ET_UT = 0.0; // Difference between ephemeris time 
                    // and universal time in [s]

// Chebyshev approximation of solar and lunar coordinates
Cheb3D LunarEph(MoonEqu, Degree, Interval);
Cheb3D SolarEph(SunEqu,  Degree, Interval);

Vec3D r_Obs_gc; // Observer's geocentric position [Earth radii]


//------------------------------------------------------------------------------
//
// GetInput: Prompts user for input
//
// Output:
//
//   TNewMoon  Approximate time of New Moon in jul. cent. since J2000
//   ET_UT     Difference betweeen ephemeris time and universal time in [s]
//   ObsLambda Geographical longitude of the observer (east pos.) in [rad]
//   ObsPhi    Geographical latitude of the observer in [rad]
//
//------------------------------------------------------------------------------
void GetInput( double& TNewMoon, double& ET_UT, 
               double& ObsLambda, double&ObsPhi )
{
  //
  // Variables
  //
  int    year, month, day;
  double hour;
  bool   valid;


  cout << " Date of new Moon (yyyy mm dd hh.hhh)         ... ";
  cin >> year >> month >> day >> hour; cin.ignore(81,'\n');

  TNewMoon = ( Mjd(year, month, day) + hour/24.0 - MJD_J2000) / 36525.0;

  ETminUT (TNewMoon, ET_UT, valid);
  if ( valid ) {
    cout << " Difference ET-UT (proposal:" << fixed << setw(6) 
         << setprecision (1) << ET_UT << " sec)       ... ";
    cin >> ET_UT; cin.ignore(81,'\n');
  }
  else {
    cout << " Difference ET-UT (sec)" << setw(27) << right << "... "; 
    cin >> ET_UT; cin.ignore(81,'\n');
  }

  cout << " Observer's coordinates: East longitude [deg] ... ";
  cin >> ObsLambda; cin.ignore(81,'\n');

  ObsLambda *= Rad;
  
  cout << "                         latitude [deg]       ... ";
  cin >> ObsPhi; cin.ignore(81,'\n');
  
  ObsPhi *= Rad;
}


//------------------------------------------------------------------------------
//
// Bessel: Computes the orientation of the fundamental plane, the 
//         coordinates of the shadow center and the parameters of the 
//         shadow cone
//
// Input:
//
//   T_UT        Time in Julian centuries UT since J2000
//
// Output:
//
//   IJK         Unit vectors of the fundamental plane (equatorial)
//   XSh         Coordinates of the shadow center on the 
//   YSh           fundamental plane in [Earth radii]
//   F1          Half angle of the penumbra cone in [rad]
//   L1          Radius of penumbra cone on the fund. plane [Earth radii]
//   F2          Half angle of the umbra cone [rad]
//   L2          Radius of umbra cone on the fund. plane [Earth radii]
//
// Globals used:
//
//   ET_UT       Difference between ephemeris time and universal time in [s]
//   SolarEph    Chebyshev polynomials of solar coordinates
//   LunarEph    Chebyshev polynomials of lunar coordinates
//
//------------------------------------------------------------------------------
void Bessel ( double T_UT,
              Mat3D& IJK, double& XSh, double& YSh,
              double& F1, double& L1, double& F2, double& L2 )
{
  //
  // Constants
  //
  const double tau_Sun  = 8.32 / ( 1440.0*36525.0); // 8.32 min  [cy]


  //
  // Variables
  //
  double T_ET, DistMS, ZSh;
  Vec3D  r_Sun, r_Moon, r_MS, i, j, k;


  T_ET = T_UT + ET_UT/(86400.0*36525.0);
  
  // Calculate solar and lunar equatorial coordinates in Earth radii
  r_Sun  = SolarEph.Value(T_ET-tau_Sun) * AU / R_Earth;
  r_Moon = LunarEph.Value(T_ET)              / R_Earth;

  // The Sun-Moon direction unit vector becomes vector k
  // of the triade defining the fundamental plane
  r_MS = r_Sun - r_Moon;
  DistMS = Norm(r_MS);
  k = r_MS / DistMS;
  
  // The unit vector i lies in the equatorial
  // plane and is perpendicular to k
  i = Polar( k[phi] + pi/2.0, 0 );
  
  // The unit vector j is perpendicular to k and i
  j = Cross(k, i);

  IJK = Mat3D(i, j, k);

  // Lunar coordinates in the fundamental plane system
  XSh = Dot(r_Moon, i);
  YSh = Dot(r_Moon, j);
  ZSh = Dot(r_Moon, k);

  // Shadow cones
  F1 = asin ( (k_Sun + k_Moon) / DistMS );
  F2 = asin ( (k_Sun - k_Moon) / DistMS );

  L1 = ZSh * (k_Sun + k_Moon) / DistMS + k_Moon;
  L2 = ZSh * (k_Sun - k_Moon) / DistMS - k_Moon;
}
   

//------------------------------------------------------------------------------
//
// Observer: Computes the observer's location on the fundamental plane
//
// Input:
//
//   T_UT        Time in Julian centuries UT since J2000
//   IJK         Unit vectors of the fundamental plane (equatorial)
//
// <return>:     Observer's location on the fundamental plane [Earth radii]
//
//
// Globals used:
//
//   r_Obs_gc    Geocentric coordinates of the observer           
//
//------------------------------------------------------------------------------
Vec3D Observer (double T_UT, const Mat3D& IJK)
{
  //
  // Variables
  //
  Vec3D r_Obs_eq;

  
  // Compute equatorial coordinates of observer
  r_Obs_eq = R_z(-GMST(T_UT*36525.0+MJD_J2000)) * r_Obs_gc;
  
  return Transp(IJK)*r_Obs_eq; 
}


//------------------------------------------------------------------------------
//
// ShadowDist:
//
//   Computes the shadow distance function f(t)=D(t)**2-L(t)**2, where D is  
//   the distance of the observer from the shadow axis, while L is the       
//   radius of the shadow cone. For f(t)=0 the observer touches the shadow   
//   cone (unit Earth radii squared).                                        
//
// Input:
//
//   T_UT        Time in Julian centuries UT since J2000
//
// <return>:     Distance observer - center of shadow in [Earth radii]
//
// Globals used:
//
//   ShadowType  Umbra or penumbra                                
//   ET_UT       Difference between ephemeris time and universal time [s]
//   r_Obs_gc    Geocentric coordinates of the observer           
//   SolarEph    Chebyshev polynomials of solar coordinates
//   LunarEph    Chebyshev polynomials of lunar coordinates
//
//------------------------------------------------------------------------------
double ShadowDist (double T_UT)
{
  //
  // Variables
  //
  Mat3D  IJK;
  double XSh, YSh, F1, L1, F2, L2, LL;


  // Fundamental plane coordinates of the shadow and the observer
  Bessel ( T_UT, IJK, XSh, YSh, F1, L1, F2, L2 );
  Vec3D r_Obs = Observer(T_UT, IJK);

  // Shadow radius at the observer's place (LL<0 for a total eclipse)
  switch ( ShadowType ) {
    case Penumbra: LL = L1 - r_Obs[z] * tan(F1); break;
    case Umbra:    LL = L2 - r_Obs[z] * tan(F2);
  }

  return (XSh-r_Obs[x])*(XSh-r_Obs[x])+(YSh-r_Obs[y])*(YSh-r_Obs[y])-LL*LL;
}



//------------------------------------------------------------------------------
//
// Contacts:
//
//   Determines the phase, magnitude and times of contacts of a solar        
//   eclipse for a specified location of the observer based on the time      
//   of new Moon.                                                            
//
// Description:                                                              
//                                                                           
//   For computing the times of contact the roots of f(t)=D(t)**2-L(t)**2
//   are determined, where D is the observer's distance from the shadow 
//   axis, L is the radius of the shadow cone. 
//   For f(t)=0 the observer lies on the surface of the shadow cone and 
//   sees the touching rims of Sun and Moon.    
//                                                                           
//   Using quadratic interpolation, the times of 1st and 4th contacts are 
//   first computed, at which the observer touches the penumbra cone. 
//   Simultaneously the time of maximum eclipse is obtained, at which the
//   distance between observer and shadow axis attains a minimum.                               
//                                                                           
//   If the eclipse is at least a partial one, Contacts computes the 
//   magnitude of the eclipse and the obscuration of the Sun by the Moon
//   for the time of maximum eclipse.                                                          
//                                                                           
//   In case the observer's distance from the shadow axis and the umbra
//   cone diameter at the time of maximum eclipse result in a total or
//   annular eclipse, the times of 2nd and 3rd contact may subsequently be
//   derived. Since the duration of the total or annular phase never 
//   exceeds a value of approximately 13 minutes, the 2nd and 3rd contacts
//   take place in a corresponding interval before or after the maximum. 
//   Making use of this fact, one may use a method like the regula falsi 
//   or the Pegasus method, which converge more rapidly than quadratic 
//   interpolation, for computing the contact times.                                                        
//
// Input:
//
//   TNewMoon    Approximate time of New Moon
//
// <return>:     Struct holding the local circumstances
//
// Globals used:
//
//   ET_UT       Difference between ephemeris time and universal time [s]
//   r_Obs_gc    Geocentric coordinates of the observer           
//   SolarEph    Chebyshev polynomials of solar coordinates
//   LunarEph    Chebyshev polynomials of lunar coordinates
//
//------------------------------------------------------------------------------

LocCircs Contacts(double TNewMoon)
{
  //
  // Constants
  //
  const double Range = 5.0;       // Search interval = +/-(Range+dTab) [h]
  const double dTab  = 0.1;       // Step size [h]
  const double hCent  = 876600.0; // 24*36525 hours per century
  const double dT  = 0.25/hCent;  // ~15 min. in Julian centuries
  const double eps = 1.0E-10;     // Accuracy for contact times ( ~0.3 s )

    
  //
  // Variables
  //
  LocCircs  LC;
  int       nContactsFound, nRoots;
  bool      Success;
  double    T_UT, PDminus, PD0, PDplus;
  double    xe, ye, Root1, Root2;
  double    XSh, YSh, F1, L1, F2, L2, LL1, LL2;
  double    m, A, B, C, S;
  Mat3D     IJK;


  // Prepare for search of outer contacts
  for (int i = 0; i < 4; i++)
    LC.Times[i] = TNewMoon;
  
  nContactsFound = 0;
  nRoots = 0;
  
  T_UT = TNewMoon + (-Range-dTab) / hCent;
  ShadowType = Penumbra;
  PDplus = ShadowDist(T_UT);
  
  T_UT = TNewMoon + (-Range-2.0*dTab) / hCent;

  
  do {  // try to find outer contacts

    T_UT += 2.0 * dTab/hCent; // Compute next time step

    
    // Compute square of the shadow distance at times T_UT-dTab,
    // T_UT and T_UT+dTab and interpolate
    PDminus = PDplus;
    PD0     = ShadowDist (T_UT);
    PDplus  = ShadowDist (T_UT + dTab/hCent);
    Quad (PDminus, PD0, PDplus, xe, ye, Root1, Root2, nRoots);


    // Store number of contacts found so far and compute contact times
    nContactsFound += nRoots;

    switch (nRoots) {
      case 1:
        if ( nContactsFound == 1 )
          LC.Times[0] = T_UT + (dTab*Root1)/hCent; // 1st contact
        else
          LC.Times[3] = T_UT + (dTab*Root1)/hCent; // 4th contact
      break;

      case 2:
        LC.Times[0] = T_UT + (dTab*Root1)/hCent; // 1st contact
        LC.Times[3] = T_UT + (dTab*Root2)/hCent; // 4th contact
      break;
    }


    // If the minimum of the shadow distance lies within the present
    // interval, a sufficiently accurate value of the time of maximum
    // eclipse is determined from its position without further iteration.
    if ( (-1.0 < xe) && (xe < 1.0) )
      LC.TMax = T_UT + dTab*xe/hCent;

  } while ( (TNewMoon+Range/hCent >= T_UT) && (nContactsFound != 2));


  // Determine type of eclipse
  if ( nContactsFound == 0 )
    LC.Phase = noEclipse; // No eclipse found
  else
    LC.Phase = partial;   // Eclipse is at least partial

  
  //----------------------------------------------------------------------------
  //
  // Degree of obscuration and magnitude of the eclipse
  //
  //----------------------------------------------------------------------------
  if ( LC.Phase > noEclipse ) {

    // Coordinates of shadow and observer on the fundamental plane
    // at maximum eclipse
    Bessel ( LC.TMax, IJK, XSh, YSh, F1, L1, F2, L2 );
    Vec3D r_Obs = Observer (LC.TMax, IJK);

    // Distance between observer and shadow axis
    m = sqrt ( 
      (r_Obs[x] - XSh) * (r_Obs[x] - XSh) + 
      (r_Obs[y] - YSh) * (r_Obs[y] - YSh) );
    
    // Penumbral and umbral radius at the observer's place
    // (L2<0 for a total eclipse!)
    LL1 = L1 - r_Obs[z] * tan(F1);
    LL2 = L2 - r_Obs[z] * tan(F2);

    // Eclipse type
    if ( m <  LL2 )
      LC.Phase = annular;
    if ( m < -LL2 )
      LC.Phase = total;

    // Eclipse magnitude
    if ( LC.Phase == partial )
      LC.Mag = (LL1-  m)/(LL1+LL2); // Penumbra region
    else
      LC.Mag = (LL1-LL2)/(LL1+LL2); // Umbra region

    // Compute degree of obscuration
    switch ( LC.Phase ) {

      case noEclipse: LC.Obsc = 0.0; break;

      case partial:       
        
        B = acos( (LL1*LL2+m*m)/(m*(LL1+LL2)) );
        C = acos( (LL1*LL1+LL2*LL2-2.0*m*m) /
                  (LL1*LL1-LL2*LL2) );
        A = pi - B - C;
        S = (LL1-LL2)/(LL1+LL2);
        LC.Obsc = (S*S*A + B - S*sin(C)) / pi;
        break;

      case annular:

        S = (LL1-LL2)/(LL1+LL2);
        LC.Obsc = S*S;
        break;

      case total:
        LC.Obsc = 1.0;
    }
  
  } // if ( LC.Phase > noEclipse )


  
  //----------------------------------------------------------------------------
  // Use the Pegasus method to locate the times of 2nd and 3rd contact
  // starting at the time of maximum eclipse. Search interval
  // [LC.TMaximum-DT,LC.TMaximum] for 2nd contact,
  // [LC.TMaximum,LC.TMaximum+DT] for 3rd contact.
  //----------------------------------------------------------------------------
  if ( LC.Phase > partial ) {

    ShadowType = Umbra;

    Pegasus(ShadowDist, LC.TMax - dT, LC.TMax, eps, LC.Times[1], Success);
    Pegasus(ShadowDist, LC.TMax, LC.TMax + dT, eps, LC.Times[2], Success);
  }

  return LC;
}


//------------------------------------------------------------------------------
//
// PosAngles:
//
//   Computes the position angles w.r.t. North (standard definition) and     
//   with respect to the local vertical.                                     
//
// Input:
//
//   TContacts[] Contact times (in Julian cent. since J2000 UT)
//   Phase       Maximum phase of the eclipse
//
// Output:
//
//   P[]         Position angle measured from north [rad]
//   V[]         Position angle with respect to the vertical [rad]
//
// Globals used:
//
//   ET_UT       Difference between ephemeris time and universal time [s]
//   r_Obs_gc    Geocentric coordinates of the observer           
//   SolarEph    Chebyshev polynomials of solar coordinates
//   LunarEph    Chebyshev polynomials of lunar coordinates
//
//------------------------------------------------------------------------------
void PosAngles ( double TContacts[], enPhase Phase,
                 double P[], double V[] )
{
  //
  // Variables
  //
  int    Contact;
  double T_UT, XSh, YSh, F1, L1, F2, L2, LL;
  Mat3D  IJK;


  for ( Contact = 0; Contact <= 3; Contact++ ) {

    if (   (Phase > partial) || 
         ( (Phase==partial) && ( (Contact==0) || (Contact==3))) ) {

      T_UT = TContacts[Contact];

      // Coordinates of shadow and observer on the fundamental plane
      Bessel ( T_UT, IJK, XSh, YSh, F1, L1, F2, L2 );
      Vec3D r_Obs = Observer (T_UT, IJK);

      // Shadow radius at the observer's place
      switch (Contact) {
        case 0: 
        case 3: LL = L1 - r_Obs[z] * tan(F1); break;
        case 1: 
        case 2: LL = L2 - r_Obs[z] * tan(F2);
      }

      // Position angles w.r.t North and w.r.t the vertical
      P[Contact] = Modulo(atan2((XSh  - r_Obs[x])/LL, (YSh - r_Obs[y])/LL), pi2);
      V[Contact] = Modulo(P[Contact] - atan2(r_Obs[x], r_Obs[y]), pi2);
    }
    else {
      P[Contact] = 0.0;
      V[Contact] = 0.0;
    }
  }
}


//------------------------------------------------------------------------------
//
// PrintResults: Screen output of local circumstances
//
// Input:
//
//   LC          Structure holding local circumstances information
//   P[]         Position Angles w.r.t. North
//   V[]         Position Angles w.r.t. the vertical
//
//------------------------------------------------------------------------------
void PrintResults (const LocCircs& LC, const double P[], const double V[])
{
  cout << endl;

  // Eclipse type, magnitude and degree of ocultation
  switch ( LC.Phase ) {

    case noEclipse:
      cout << " The eclipse is not visible at the selected site." << endl;
    break;

    case partial:
      cout << " Partial eclipse with M=" << fixed << setw(5) 
           << setprecision(3) << LC.Mag  << " ("  << setw(4) 
           << setprecision(2) << LC.Obsc << ")."  << endl;
    break;

    case annular:
      cout << " Annular eclipse with M=" << fixed << setw(5) 
           << setprecision(3) << LC.Mag  << " ("  << setw(4) 
           << setprecision(2) << LC.Obsc << ")."  << endl;
    break;

    case total:
      cout << " Total eclipse with M=" << fixed << setw(5) 
           << setprecision(3) << LC.Mag << "."  << endl;
  }

  // Print contact times and position angles, if the eclipse
  // is at least partial at the place of the observer
  if ( LC.Phase > noEclipse ) {

    double MJD_UT = LC.TMax*36525.0+MJD_J2000;
    DateTime DTMax (MJD_UT, HHMMSS);
    cout << " Maximum at " << DTMax << " UT." << endl << endl;

    cout << setw(58) << right << "h  m  s [UT]    P [o]     V [o]" << endl;

    for ( int Contact=0; Contact <= 3; Contact++ ) {

      if ( (Contact==0) || (Contact==3) || (LC.Phase > partial) ) {

        switch (Contact) {
          case 0: cout << " 1st"; break;
          case 1: cout << " 2nd"; break;
          case 2: cout << " 3rd"; break;
          case 3: cout << " 4th";
        }
        cout << " contact: ";

        MJD_UT = LC.Times[Contact] * 36525.0 + MJD_J2000;
        DateTime DTContact(MJD_UT, HHMMSS);
        cout << DTContact;

        cout << setw(15) << setprecision(0) << right << Deg*P[Contact];
        cout << setw(10) << setprecision(0) << right << Deg*V[Contact];
        cout << endl;
      }
    }

  }
  cout << endl;
}


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------
void main() 
{
  //
  // Variables
  //
  double TNewMoon  = 0.0;
  double ObsLambda = 0.0; // Observers geographic longitude; east pos. [rad]
  double ObsPhi    = 0.0; // Observers geographic latitude [rad]
  double P[4];            // Arrays for position angles of...
  double V[4];            // ...contacts (w.r.t. North / Vertical) in [rad]


  // Title
  cout << endl
       << "      ECLTIMER: local circumstances of solar eclipses" << endl
       << "        (c) 1999 Oliver Montenbruck, Thomas Pfleger  " << endl
       << endl;

  GetInput (TNewMoon, ET_UT, ObsLambda, ObsPhi); // Input data

  // Geocentric coordinates of the observer (in Earth radii)
  r_Obs_gc = Site (ObsLambda, ObsPhi) / R_Earth;

  LocCircs LC = Contacts (TNewMoon);    // Compute contact times...
  
  PosAngles (LC.Times, LC.Phase, P, V); // ...and position angles

  PrintResults (LC, P, V);
}
