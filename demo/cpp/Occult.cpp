//------------------------------------------------------------------------------
//
// Occult program (Occult.cpp)
//
// Prediction of stellar occultations by the Moon
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
#include "APC_Phys.h"
#include "APC_PrecNut.h"
#include "APC_Spheric.h"
#include "APC_Sun.h"
#include "APC_Time.h"
#include "APC_VecMat3D.h"

#ifdef __GNUC__  // GNU C++ adaptation
#include "GNU_iomanip.h"
#endif

using namespace std;

#ifdef _MSVC60  // Bugfix (forward-declarations)
class Star;
istream& operator >> (istream& is, Star& aStar);
class Event;
ostream& operator << (ostream& os, Event& anEvent);
#endif

//
// Constants
//
const double Interval  = 1.0/36525.0;    // 1d [cy]


//
// Global objects and variables
//
Cheb3D ChebMoonEqu(MoonEqu, 10, Interval);
double RA_Star;


//
// Class Star
//
class Star {

  public:

    static void SetCatEquinox( double T );
    static void SetCatEpoch( double T );
    static void SetEpoch( double T );

    Vec3D  Apparent();  // Calculate apparent place of star

    inline char*  Name() { return name; };
    inline double Mag() { return mag; };

    friend istream& operator >> (istream& is, Star& aStar);

  private:

    double RA_Cat;               // Catalogue right ascension in [rad]
    double Dec_Cat;              // Catalogue declination in [rad]
    double mu_RA;                // Proper motion Ra in [rad/cy]
    double mu_Dec;               // Proper motion Dec in [rad/cy]
    double mag;                  // Visual brightness in [mag]
    char   name[12];             // Name

    static double T_CatEquinox;  // Catalogue equinox
    static double T_CatEpoch;    // Catalogue epoch
    static double T_Epoch;       // Current epoch

    static Mat3D  PN;            // Precession/nutation matrix
    static Vec3D  v_Earth;       // Earth velocity vector 

};


// Static class variables 
double Star::T_CatEquinox;  // Catalogue equinox
double Star::T_CatEpoch;    // Catalogue epoch
double Star::T_Epoch;       // Current epoch
Mat3D  Star::PN;            // Precession/nutation matrix
Vec3D  Star::v_Earth;       // Earth velocity vector 


void Star::SetCatEquinox (double T)
{
  T_CatEquinox = T;
}


void Star::SetCatEpoch (double T)
{
  T_CatEpoch = T;
}


void Star::SetEpoch (double T)
{
  T_Epoch = T;  // Save epoch
  
  // Transformation from mean equator and equinox of catalogue
  // to true equator and equinox of date
  PN = NutMatrix(T)*PrecMatrix_Equ(T_CatEquinox, T);

  // Earth velocity vector w.r.t. mean equator and equinox of date
  // in units of speed of light 
  v_Earth = Ecl2EquMatrix(T)*KepVelocity(Earth, T) / c_light;
}


Vec3D Star::Apparent()
{
  // Proper motion
  double Ra  = RA_Cat  + mu_RA *(T_Epoch-T_CatEpoch);
  double Dec = Dec_Cat + mu_Dec*(T_Epoch-T_CatEpoch);

  // Precession, nutation
  Vec3D e = PN * Vec3D(Polar(Ra, Dec));

  // Aberration
  e = e + v_Earth;
  e = e / Norm(e);

  return e;
};


// Input operator
istream& operator >> (istream& is, Star& aStar)
{
  int     h,m,deg,min;
  double  s,sec,mu_a,mu_d;
  char    c,sign;

  is >> h >> m >> s >> mu_a;
  is.get(c).get(c).get(c).get(sign);
  is >> deg >> min >> sec >> mu_d;
  
  is >> aStar.mag; 
  is.get(c).get(c).get(c).get(aStar.name,12);
  is.ignore(81,'\n'); 

  aStar.RA_Cat  = Rad*15.0*Ddd(h,m,s);
  aStar.Dec_Cat = Rad*Ddd(deg,min,sec);

  if (sign=='-') aStar.Dec_Cat = -aStar.Dec_Cat;

  aStar.mu_RA   = Rad*15.0*mu_a/3600.0;
  aStar.mu_Dec  = Rad*mu_d/3600.0;

  c=is.peek(); // Required to enable EOF check

  return is;
}


//
// Class Event
//

enum enEvent {dummy=-1, in=0, out=1};  // Identifies disappearance or reapperance

class Event {

  public:

    Event();
    Event( const Star& aStar, double MjdUT, enEvent InOut, bool visible,
           double PosAng, double a, double b );
    
    inline bool IsVisible() { return m_visible; };

    friend ostream& operator << (ostream& os, Event& anEvent);

  private:

    bool      m_visible;
    Star      m_Star;
    double    m_MjdUT;
    enEvent   m_InOut;
    double    m_PosAng;
    double    m_a,m_b;
};


//
// Constructors
//
Event::Event()
: m_visible(false),
  m_InOut(dummy)
{ 
}

Event::Event( const Star& aStar, double MjdUT, enEvent InOut, bool visible,
              double PosAng, double a, double b )
 : m_visible(visible),
   m_Star(aStar),
   m_MjdUT(MjdUT),
   m_InOut(InOut),
   m_PosAng(PosAng),
   m_a(a),
   m_b(b)
{
}


//
// Output operator
//
ostream& operator << (ostream& os, Event& anEvent)
{
  if ( anEvent.m_InOut != dummy ) {   
    
    // Output true events only
    os << " " << DateTime(anEvent.m_MjdUT)
       << "  " << anEvent.m_Star.Name()
       << fixed << setprecision(1)
       << setw(5) << anEvent.m_Star.Mag()
       << "    " << ( (anEvent.m_InOut==in)? "D" : "R" ) 
       << "    " << Time(24.0*Modulo(anEvent.m_MjdUT, 1.0),HHMMSS)
       << fixed << setprecision(1)
       << setw(7) << int(Deg*anEvent.m_PosAng + 0.5)
       << setw(8) << anEvent.m_a/Deg
       << setw(7) << anEvent.m_b/Deg
       << endl;
  }
  
  return os;
}


//------------------------------------------------------------------------------
//
// ReadCatalogue: Reads star catalogue from input file
//
// Input:
//
//   Filename     Name of catalogue file to read
//
// Output:  
//
//   Stars        Array of Star objects to store catalogue
//   StarsRead    Number of stars read from catalogue file
//
// Notes: Memory allocated for Stars has to be deleted by the caller
//
//------------------------------------------------------------------------------
void ReadCatalogue (char* Filename, Star*& Stars, int& StarsRead)
{
  //
  // Variables
  //
  ifstream inp;
  double   Y_Epoch, Y_Eqx;
  int      i;

  inp.open(Filename); 

  
  // Read and store epoch and equinox
  inp >> Y_Epoch >> Y_Eqx;

  Star::SetCatEpoch((Y_Epoch-2000.0)/100.0);
  Star::SetCatEquinox((Y_Eqx-2000.0)/100.0);


  // Read stars
  StarsRead = 0;
  inp >> StarsRead; inp.ignore(81,'\n');
  
  Stars = new Star[StarsRead];  // Allocate memory for catalogue array

  if (Stars != NULL) {
    for (i=0; i<StarsRead; i++)  inp >> Stars[i]; 
  }
  else
    cerr << "Insufficient memory to read star catalogue." << endl;
  
  inp.close();

  cout << endl << " " << StarsRead 
       << " stars read from catalogue" << endl << endl;
}


//------------------------------------------------------------------------------
//
// GetInput: Prompts user for input
//
// Output:
//
//   MjdStart  Start time for search for occultations
//   MjdEnd    End time for search for occultations
//   ET_UT     Difference ephemeris time - universal time in [s]
//   R_Obs     Geocentric position of the observer in [km]
//
// Note: MjdStart and MjdEnd as Modified Julian Dates
//
//------------------------------------------------------------------------------
void GetInput ( double& MjdStart, double& MjdEnd,  double& ET_UT,
                Vec3D& R_Obs )
{
  //
  // Variables
  //
  int     year,month,day;
  bool    valid;
  double  Lambda, Phi;

  // Query prediction interval
  cout << " Period of time for prediction of occultations" << endl;

  cout << "  first date (yyyy mm dd)                     ... ";
  cin >> year >> month >> day; cin.ignore(81,'\n');
  MjdStart = Mjd(year,month,day);

  cout << "  last  date (yyyy mm dd)                     ... ";
  cin >> year >> month >> day; cin.ignore(81,'\n');
  MjdEnd = Mjd(year,month,day);

  
  // Query difference between ephemeris time and universal time
  ETminUT ( (0.5*(MjdStart+MjdEnd)-MJD_J2000)/36525.0, ET_UT, valid );

  if ( valid ) {
    cout << " Difference ET-UT (proposal:" << fixed << setw(6) 
         << setprecision (1) << ET_UT << " sec)       ... ";
    cin >> ET_UT; cin.ignore(81,'\n');
  }
  else {
    cout << " Difference ET-UT (sec)" << setw(27) << right << "... "; 
    cin >> ET_UT; cin.ignore(81,'\n');
  }

  
  // Query geographic coordinates
  cout << " Observer's coordinates: East longitude [deg] ... ";
  cin >> Lambda; cin.ignore(81,'\n');
  cout << "                         latitude       [deg] ... ";
  cin >> Phi; cin.ignore(81,'\n'); 
 
  Lambda *= Rad; Phi *= Rad;
  R_Obs = Site (Lambda, Phi);
}


//------------------------------------------------------------------------------
//
// RA_Diff: Goal function for iteration of conjunction instant
//
// Input:
//
//   T        Time in Julian centuries ET since J2000
//
// <return>:  Difference in Ra between Moon and star in [rad]
//
// Globals used:
//
//   ChebMoonEqu  Chebyshev approximation of lunar coordinates
//   RA_Star      Right ascension of star currently examined in [rad]
//
//------------------------------------------------------------------------------
double RA_Diff (double T)
{
  Vec3D r_Moon = ChebMoonEqu.Value(T);   // Lunar position at time T
  return Modulo ( r_Moon[phi] - RA_Star + pi, pi2 ) - pi;
}


//------------------------------------------------------------------------------
//
// Conjunct: Checks for possibility of a lunar occultation 
//
// Input:
//
//   T1       Begin of search interval
//   T2       End of search interval
//   RA1      Right ascension of the Moon at T1 in [rad]
//   RA2      Right ascension of the Moon at T2 in [rad]
//   e        Direction of star to check (equatorial coordinates)
//
// Output:    
// 
//   Conj     Flag indicating whether an occultation is at least likely
//   T_Conj   Time of conjunction
//
// Globals used:
//
//   ChebMoonEqu  Chebyshev approximation of lunar coordinates
//   RA_Star      Right ascension of star in [rad]
//
//------------------------------------------------------------------------------
void Conjunct ( double T1,  double T2, double RA1, double RA2, Vec3D& e,
                bool& Conj, double& T_Conj )
{
  //
  // Constants
  //
  const double eps = Rad*1.0e-4;      // Required accuracy in Ra [rad]

  //
  // Variables
  //
  Vec3D r_Moon;
  T_Conj = 0.0;                       // Default value
  bool Success = false;

  
  if ( RA2<RA1 ) RA2+=pi2;            // Map to [RA1,RA1+2*pi[

  RA_Star = e[phi];
  if ( RA_Star<RA1 ) RA_Star+=pi2;    // Map to [RA1,RA1+2*pi[

  Conj = ( (RA1<=RA_Star) && (RA_Star<=RA2) );

  if ( Conj ) {
    Pegasus(RA_Diff, T1, T2, eps, T_Conj, Success);
    r_Moon = ChebMoonEqu.Value(T_Conj); 
    Conj = ( Norm( r_Moon-e*Dot(r_Moon,e) ) < 1.5*R_Earth );
  }
}


//------------------------------------------------------------------------------
//
// FG: Compute relative fundamental plane coordinates
//
// Input:
//
//   T_ET     Time ET in Julian centuries since J2000
//   ET_UT    Difference betweeen ephemeris time and universal time in [s]
//   e        Direction of star to check (equatorial coordinates)
//   R_Obs    Geocentric position of the observer in [km]
//
// Output:    
// 
//   f, g     Relative fundamental plane coordinates in [km]
//   s        Square of distance between shadow axis and observer in [km]
//
// Globals used:
//
//   ChebMoonEqu  Chebyshev approximation of lunar coordinates
//
//------------------------------------------------------------------------------
void FG ( double T_ET, double ET_UT, const Vec3D& e, const Vec3D& R_Obs, 
          double& f, double& g, double& s ) 
{
  //
  // Variables
  //
  double  MjdUT; 
  Vec3D   e_x, e_y, r_Moon, r_Obs;


  // Fundamental plane unit vectors
  e_x = Cross ( Vec3D(0,0,1), e );
  e_x = e_x / Norm(e_x);
  e_y = Cross ( e, e_x );

  // Lunar position
  r_Moon = ChebMoonEqu.Value(T_ET);

  // Observer's position
  MjdUT = (36525.0*T_ET+MJD_J2000) - ET_UT/86400;
  r_Obs = R_z(-GMST(MjdUT)) * R_Obs;

  // Relative fundamental plane coordinates
  f = Dot ( e_x, r_Moon - r_Obs );
  g = Dot ( e_y, r_Moon - r_Obs );

  s = f*f + g*g - R_Moon*R_Moon;
}


//------------------------------------------------------------------------------
//
// Contact: Compute position angle, longitude and latitude coefficients and
//          check for visibility (concerning altitude and sky brightness)
//
// Input:
//
//   T_ET       Time ET in Julian centuries since J2000
//   ET_UT      Difference betweeen ephemeris time and universal time in [s]
//   e          Direction of star to check (equatorial coordinates)
//   R_Obs      Geocentric position of the observer in [km]
//   InOut      Indicates disappearance or reappearance
//   mag        Magnitude of the star in [mag]
//
// Output:    
// 
//   PosAngle   Position angle w.r.t. North in [rad]
//   a          Longitude coefficient in [min/']
//   b          Latitude coefficient in [min/']
//   Visible    Flag indicating visibility of the occultation
//
// Globals used:
//
//   ChebMoonEqu  Chebyshev approximation of lunar coordinates
//
//------------------------------------------------------------------------------
void Contact ( double T_ET, double ET_UT, const Vec3D& e, const Vec3D& R_Obs,
               enEvent InOut, double mag,
               double& PosAngle, double& a, double& b, bool& Visible ) 
{
  //
  // Constants
  //
  const double cy = 36525.0*1440.0;   // Minutes per century
  const double dt = 1.0;              // Time step [min]


  //
  // Variables
  //
  double  MjdUT; 
  Vec3D   e_x, e_y, r_Moon, r_Obs, e_Obs, r_Sun;
  Vec3D   dt_dr, dr_dlambda, dr_dphi;
  double  s_0, f, g, ff, gg, fdot, gdot;
  double  Elev, ElevSun, PosAngSun, DPosAng, i;
  bool    BrightLimb;


  // Fundamental plane unit vectors
  e_x = Cross ( Vec3D(0,0,1), e );
  e_x = e_x / Norm(e_x);
  e_y = Cross ( e, e_x );


  // Lunar position
  r_Moon = ChebMoonEqu.Value(T_ET);

  
  // Observer's geocentric position
  MjdUT = (36525.0*T_ET+MJD_J2000) - ET_UT/86400;

  r_Obs = R_z(-GMST(MjdUT)) * R_Obs;
  e_Obs = r_Obs / Norm(r_Obs);

  
  // Altitude of star above the horizon
  Elev = asin ( Dot(e_Obs,e) ); 

  
  // Relative fundamental plane coordinates
  f = Dot ( e_x, r_Moon - r_Obs );
  g = Dot ( e_y, r_Moon - r_Obs );


  // Position angle
  PosAngle = Modulo(atan2(-f, -g), pi2);


  // Time derivative of fundamental plane coordinates
  FG ( T_ET+dt/cy, ET_UT, e, R_Obs, ff, gg, s_0 );
  fdot = (ff-f)/dt;  // [km/min]
  gdot = (gg-g)/dt;

  
  // Derivatives of observer position w.r.t. geocentric coordinates
  dr_dlambda = Cross ( Vec3D(0,0,1), r_Obs );
  dr_dphi    = Cross ( r_Obs, dr_dlambda / Norm(dr_dlambda) );


  // Compute longitude and latitude coefficients
  dt_dr = (f*e_x+g*e_y) / (f*fdot+g*gdot);

  a = Dot ( dt_dr, dr_dlambda );
  b = Dot ( dt_dr, dr_dphi    );
  

  // Compute solar coordinates, altitude and position angle 
  r_Sun = AU*SunEqu(T_ET);
  ElevSun = asin ( Dot ( e_Obs, r_Sun/Norm(r_Sun) ) );
  PosAngSun = PosAng(r_Moon-r_Obs, r_Sun-r_Moon);
  DPosAng = Modulo(PosAngle-PosAngSun+pi, pi2) - pi;
  BrightLimb = ( fabs(DPosAng) < pi/2.0 );


  // Phase angle of Moon
  i = acos ( Dot(-r_Moon, r_Sun-r_Moon) / (Norm(r_Sun-r_Moon)*Norm(r_Moon)) );


  // Visibility conditions
  Visible = ( // Minimum altitude of star 
              ( Elev > +5.0*Rad )                    ||
              ( Elev > +2.0*Rad ) && ( mag <= 1.9 ) 
            ) &&
            ( // Minimum depression of the Sun
              ( ElevSun < -6.0*Rad ) && ( mag <= 7.5 ) ||
              ( ElevSun < -3.0*Rad ) && ( mag <= 5.5 ) ||
              ( ElevSun < -0.5*Rad ) && ( mag <= 4.5 ) ||
                                        ( mag <= 1.9 ) 
            ) &&
            ( // Lunar phase
              ( i > 36.0*Rad )                   ||
              ( i > 24.0*Rad ) && ( mag <= 6.5 ) ||
              ( i > 12.0*Rad ) && ( mag <= 5.5 ) ||
                                  ( mag <= 3.5 ) 
            ) &&
            ( // Limb illumination
              ( !BrightLimb ) && (InOut==in )               ||
              ( !BrightLimb ) && (InOut==out) && (mag<=6.5) ||
              (  BrightLimb ) && (InOut==in ) && (mag<=4.5) ||
              (  BrightLimb ) && (InOut==out) && (mag<=3.5)  
            );
}


//------------------------------------------------------------------------------
//
// Examine: Checks for possibility of a lunar occultation 
//
// Input:
//
//   T1       Begin of search interval
//   T2       End of search interval
//   RA1      Right ascension of the Moon at T1 in [rad]
//   RA2      Right ascension of the Moon at T2 in [rad]
//   ET_UT    Difference ephemeris time - universal time in [s]
//   R_Obs    Geocentric position of the observer in [km]
//   aStar    Star to examine
//
// Output:    
// 
//   In       Disappearance event
//   Out      Reappearance event
//
//------------------------------------------------------------------------------
void Examine ( double T1, double T2, double RA1, double RA2, double ET_UT, 
               const Vec3D& R_Obs, Star& aStar, Event& In, Event& Out )
{
  //
  // Constants
  //
  const double cy = 876600.0;    // Hours per century
  const double dT = 0.25/cy;     // Step size
  const double DT = 2.25/cy;     // Search interval [-DT-dT,+DT+dT]
  
  
  //
  // Variables
  //
  Vec3D   e;
  bool    Conj, Visible;
  int     n_found, n_root, i;
  double  T_Conj, T, T_Cont[2], MjdUT;
  double  s_minus, s_0, s_plus, xe, ye, root[2];
  double  f, g, PosAng, a, b;


  In = Out = Event(); // Default values


  // Apparent coordinates of star
  e = aStar.Apparent();


  // Search for time of conjunction
  Conjunct ( T1, T2, RA1, RA2, e, Conj, T_Conj );
  if (!Conj) return;

  
  // Search for times of contact
  n_found = 0;
  T = T_Conj - DT;
  FG (T-dT, ET_UT, e, R_Obs, f, g, s_minus);
  
  while (true) {
  
    // Interpolate fundamental plane distance and find times of contact
    FG (T   , ET_UT, e, R_Obs, f, g, s_0   );
    FG (T+dT, ET_UT, e, R_Obs, f, g, s_plus);
    Quad (s_minus, s_0, s_plus, xe, ye, root[0], root[1], n_root);
    for (i=0;i<n_root;i++)
      T_Cont[n_found+i] = T + root[i]*dT;
    n_found += n_root;
    Conj = (n_found==2);
    
    if ( (Conj) || (T_Conj+DT<T) ) break; // Exit loop
    
    T+=2.0*dT; // Increment time

    s_minus = s_plus;
  };
    
  // Contacts (disappearance, reappearance)
  if (Conj) {
    for (i=in;i<=out;i++)  {
      MjdUT = 36525.0*T_Cont[i] + MJD_J2000 - ET_UT/86400.0;
      Contact ( T_Cont[i], ET_UT, e, R_Obs, enEvent(i), aStar.Mag(),
                PosAng, a, b, Visible );
      if (i==in)  In  = Event(aStar, MjdUT, in , Visible, PosAng, a, b);
      if (i==out) Out = Event(aStar, MjdUT, out, Visible, PosAng, a, b);
    }
  }
}



//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------
void main(int argc, char* argv[])
{
  //
  // Constants
  //
  const double Margin    = 3.0 / (36525.0*24.0);  // 3h [cy]


  //
  // Variables
  //
  Star*     Stars = NULL;     // Star catalogue
  int       NStars, iStar;
  double    MjdStart, MjdEnd, ET_UT, T, TEnd;
  double    RA_min, RA_max;
  Vec3D     R_Obs;
  Event     In, Out;
  char      InputFile[APC_MaxFilename] = "";
  char      OutputFile[APC_MaxFilename] = "";
  bool      FoundInputfile  = false;
  bool      FoundOutputfile = false;
  ofstream  OutFile;


  // Title
  cout << endl
       << "           OCCULT: occultations of stars by the Moon  " << endl
       << "          (c) 1999 Oliver Montenbruck, Thomas Pfleger " << endl 
       << endl;

  
  // Find catalogue file and (optional) name of output file
  GetFilenames( argc, argv, "Occult.dat", InputFile, FoundInputfile, 
                OutputFile, FoundOutputfile );

  // Terminate program if input file could not be found
  if (!FoundInputfile) {
    cerr << "Terminating program." << endl;
    exit(-1);
  }

  
  ReadCatalogue (InputFile, Stars, NStars);

  // Terminate program, if star catalogue could not be read
  if (NStars == 0) {
    cerr << "Terminating program." << endl;
    exit(-1);
  }
  
  // Get user input (search interval, observer location, etc.)
  GetInput (MjdStart, MjdEnd, ET_UT, R_Obs);


  // Redirect output if output file shall be created
  if (FoundOutputfile) {
    OutFile.open(OutputFile);
    if (OutFile.is_open())
      cout = OutFile;
  }


  // Print header
  cout 
   << endl
   << "    Date        Name      m_v   D/R      UT       Pos      a      b  "
   << endl
   << "                          mag                     deg    m/deg  m/deg"
   << endl;
  
  
  // Search occultations in subsequent time intervals (ephemeris time)
  T    = (MjdStart-MJD_J2000)/36525.0;
  TEnd = (MjdEnd  -MJD_J2000)/36525.0;

  do {

    // Chebyshev approximation of lunar coordinates 
    ChebMoonEqu.Fit(T-Margin, T+Interval+Margin);

    RA_min = ChebMoonEqu.Value(T)[phi];
    RA_max = ChebMoonEqu.Value(T+Interval)[phi];


    // Set epoch for computation of apparent coordinates
    Star::SetEpoch(T+Interval/2.0);

    // Star catalogue loop
    for (iStar=0; iStar<NStars; iStar++) {

      // Check for occultation
      Examine ( T, T+Interval, RA_min, RA_max,
                ET_UT, R_Obs, Stars[iStar], In, Out );

      // Print disappearance and reappearance, if visible
      if (In.IsVisible() ) cout << In;
      if (Out.IsVisible()) cout << Out;
    }
     
    T += Interval;  // Next interval
  }
  while (T<TEnd);

  if (OutFile.is_open()) OutFile.close();
  if (Stars!=NULL) delete[] Stars;
}
