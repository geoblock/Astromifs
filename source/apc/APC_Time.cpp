//------------------------------------------------------------------------------
//
// APC_Time.cpp
//
//    Implements functions and classes handling time
//
//------------------------------------------------------------------------------

#include <cmath>
#include <iomanip>

#include "APC_Const.h"
#include "APC_Math.h"
#include "APC_Time.h"

#ifdef __GNUC__  // GNU C++ adaptation
#include "GNU_iomanip.h"
#endif

using namespace std;


//------------------------------------------------------------------------------
//
// ETminUT: Difference ET-UT of ephemeris time and universal time
//
// Input:
//
//   T         Time in Julian centuries since J2000
//
// Output:
//
//   DTsec     ET-UT in [s]
//   valid     Flag indicating T in domain of approximation
//
// Notes: The approximation spans the years from 1825 to 2005
//
//------------------------------------------------------------------------------
void ETminUT (double T, double& DTsec, bool& valid)
{
  //
  // Variables
  //
  int    i = (int) floor(T/0.25);
  double t = T-i*0.25;


  if ( (T<-1.75) || (0.05<T) ) {
    valid = false;
    DTsec = 0.0;
  }
  else {
    valid = true;
    switch (i) {
      case -7: DTsec=10.4+t*(-80.8+t*( 413.9+t*( -572.3))); break; // 1825-
      case -6: DTsec= 6.6+t*( 46.3+t*(-358.4+t*(   18.8))); break; // 1850-
      case -5: DTsec=-3.9+t*(-10.8+t*(-166.2+t*(  867.4))); break; // 1875-
      case -4: DTsec=-2.6+t*(114.1+t*( 327.5+t*(-1467.4))); break; // 1900-
      case -3: DTsec=24.2+t*( -6.3+t*(  -8.2+t*(  483.4))); break; // 1925-
      case -2: DTsec=29.3+t*( 32.5+t*(  -3.8+t*(  550.7))); break; // 1950-
      case -1: DTsec=45.3+t*(130.5+t*(-570.5+t*( 1516.7))); break; // 1975-
      case  0: t+=0.25;
               DTsec=45.3+t*(130.5+t*(-570.5+t*( 1516.7)));        // 2000-
    }                                                              // 2005
  }
}


//------------------------------------------------------------------------------
//
// GMST: Greenwich mean sidereal time
//
// Input:
//
//   MJD       Time as Modified Julian Date
//
// <return>:   GMST in [rad]
//
//------------------------------------------------------------------------------
double GMST (double MJD)
{
  //
  // Constants
  //
  const double Secs = 86400.0;        // Seconds per day


  //
  // Variables
  //
  double MJD_0, UT, T_0, T, gmst;


  MJD_0 = floor ( MJD );
  UT    = Secs*(MJD-MJD_0);     // [s]
  T_0   = (MJD_0-51544.5)/36525.0; 
  T     = (MJD  -51544.5)/36525.0; 

  gmst  = 24110.54841 + 8640184.812866*T_0 + 1.0027379093*UT 
          + (0.093104-6.2e-6*T)*T*T;      // [sec]

  return  (pi2/Secs)*Modulo(gmst,Secs);   // [Rad]
}


//------------------------------------------------------------------------------
//
// Mjd: Modified Julian Date from calendar date and time
//
// Input:
//
//   Year      Calendar date components
//   Month
//   Day
//   Hour      Time components (optional)
//   Min
//   Sec
//
// <return>:   Modified Julian Date
//
//------------------------------------------------------------------------------
double Mjd ( int Year, int Month, int Day, 
             int Hour, int Min, double Sec )
{
  //
  // Variables
  //
  long    MjdMidnight;
  double  FracOfDay;
  int     b;


  if (Month<=2) { Month+=12; --Year;}
  
  if ( (10000L*Year+100L*Month+Day) <= 15821004L )
    b = -2 + ((Year+4716)/4) - 1179;     // Julian calendar 
  else
    b = (Year/400)-(Year/100)+(Year/4);  // Gregorian calendar 
    
  MjdMidnight = 365L*Year - 679004L + b + int(30.6001*(Month+1)) + Day;
  FracOfDay   = Ddd(Hour,Min,Sec) / 24.0; 

  return MjdMidnight + FracOfDay;
}


//------------------------------------------------------------------------------
//
// CalDat: Calendar date and time from Modified Julian Date
//
// Input:
//
//   Mjd       Modified Julian Date
//
// Output:
//
//   Year      Calendar date components
//   Month
//   Day
//   Hour      Decimal hours
//
//------------------------------------------------------------------------------
void CalDat ( double Mjd, 
              int& Year, int& Month, int& Day, double & Hour )
{
  //
  // Variables
  //
  long    a,b,c,d,e,f;
  double  FracOfDay;


  // Convert Julian day number to calendar date
  a = long(Mjd+2400001.0);

  if ( a < 2299161 ) {  // Julian calendar
    b = 0;
    c = a + 1524;
  }
  else {                // Gregorian calendar
    b = long((a-1867216.25)/36524.25);
    c = a +  b - (b/4) + 1525;
  }

  d     = long ( (c-122.1)/365.25 );
  e     = 365*d + d/4;
  f     = long ( (c-e)/30.6001 );

  Day   = c - e - int(30.6001*f);
  Month = f - 1 - 12*(f/14);
  Year  = d - 4715 - ((7+Month)/10);

  FracOfDay = Mjd - floor(Mjd);

  Hour = 24.0*FracOfDay;
}


//------------------------------------------------------------------------------
//
// CalDat: Calendar date and time from Modified Julian Date
//
// Input:
//
//   Mjd       Modified Julian Date
//
// Output:
//
//   Year      Calendar date components
//   Month
//   Day
//   Hour      Time components
//   Min
//   Sec
//
//------------------------------------------------------------------------------
void CalDat ( double Mjd, 
              int& Year, int& Month, int& Day,  
              int& Hour, int& Min, double& Sec )
{
  //
  // Variables
  //
  double Hours;


  CalDat (Mjd, Year, Month, Day, Hours); 
  
  DMS (Hours, Hour, Min, Sec);
}


//
// Implementation of class Time constructors and member functions
//

Time::Time ()
 : m_Hour(0.0),
   m_Format(HHMMSS)
{
}


Time::Time (double Hour, TimeFormat Format)
 : m_Hour(Hour),
   m_Format(Format)
{
}


void Time::Set (TimeFormat Format)
{
  m_Format = Format;
}


ostream& operator << (ostream& os, const Time& time)
{
  //
  // Variables
  //
  double Hour, S;
  int    H,M;


  os << right;

  switch ( time.m_Format ) {

    case HHh:

      os << fixed << setprecision(1);
      os << setw(4) << time.m_Hour;
      
      break;

    case HHMM:

      // round to 1 min
      Hour = floor(60.0*time.m_Hour+0.5)/60.0+0.00001;

      DMS (Hour, H, M, S);

      os.fill('0');
      os << setw(2) << H << ":" << setw(2) << M;
      os.fill(' ');
      
      break;

    case HHMMSS:

      // round to 1 sec
      Hour = floor(3600.0*time.m_Hour+0.5)/3600.0+0.0000001;

      DMS (Hour, H, M, S);

      os.fill('0');
      os << setw(2) << H << ":" << setw(2) << M << ":" << setw(2) << int(S);
      os.fill(' ');
  }
   
  return os;
}


//
// Implementation of class DateTime constructors and member functions
//

DateTime::DateTime ()
 : m_Mjd(0.0),
   m_Format(None)
{
}


DateTime::DateTime (double Mjd, TimeFormat Format)
 : m_Mjd(Mjd),
   m_Format(Format)
{
}


void DateTime::Set (TimeFormat Format)
{
  m_Format = Format;
}


ostream& operator << (ostream& os, const DateTime& DT)
{
  //
  // Variables
  //
  double MjdRound, Hours, S;
  int    Year, Month, Day;
  int    H, M;


  os << right;

  switch ( DT.m_Format ) {

    case None:

      CalDat (DT.m_Mjd, Year, Month, Day, Hours);
      
      os.fill('0');
      os << setw(4) << Year  << "/" 
         << setw(2) << Month << "/" 
         << setw(2) << Day; 
      os.fill(' ');
      
      break;

    case DDd:

      CalDat (DT.m_Mjd, Year, Month, Day, Hours);
      
      os.fill('0');
      os << setw(4) << Year  << "/" 
         << setw(2) << Month << "/" 
         << fixed << setw(3+os.precision()) << Day+Hours/24.0; 
      os.fill(' ');
      
      break;

    case HHh:

      // round to 0.1h
      MjdRound = floor(240.0*DT.m_Mjd+0.5)/240.0+0.0001;

      CalDat (MjdRound, Year, Month, Day, Hours);
      
      os.fill('0');
      os << setw(4) << Year  << "/" 
         << setw(2) << Month << "/" 
         << setw(2) << Day   << " "; 

      os << fixed << setprecision(1);
      os << setw(4) << Hours;
      os.fill(' ');
      
      break;

    case HHMM:

      // round to 1 min
      MjdRound = floor(1440.0*DT.m_Mjd+0.5)/1440.0+0.00001;

      CalDat (MjdRound, Year, Month, Day, H, M, S);
      
      os.fill('0');
      os << setw(4) << Year  << "/" 
         << setw(2) << Month << "/" 
         << setw(2) << Day   << " "; 
      os << setw(2) << H << ":"
         << setw(2) << M;
      os.fill(' ');
      
      break;

    case HHMMSS:

      // round to 1 sec
      MjdRound = floor(86400.0*DT.m_Mjd+0.5)/86400.0+0.000001;

      CalDat (MjdRound, Year, Month, Day, H, M, S);
      
      os.fill('0');
      os << setw(4) << Year  << "/" 
         << setw(2) << Month << "/" 
         << setw(2) << Day   << " "; 
      os << setw(2) << H << ":" 
         << setw(2) << M << ":" 
         << setw(2) << int(S);
      os.fill(' ');
  }
   
  return os;
}
