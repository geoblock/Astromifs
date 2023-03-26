//------------------------------------------------------------------------------
//
// APC_Time.h
// 
//    Functions and classes handling time
//
//------------------------------------------------------------------------------

#ifndef INC_APC_TIME_H
#define INC_APC_TIME_H

#include <iostream>


//
// Types
//

// Format tag for date and time output (used by DateTime class output operator)
enum TimeFormat { 
  None,   // don't output time (date only)
  DDd,    // output time as fractional part of a day
  HHh,    // output time as hours with one decimal place
  HHMM,   // output time as hours and minutes (rounded to next minute)
  HHMMSS  // output time as hours, minutes and seconds (rounded to next s)
};


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
// Notes: Approximation spans the years from 1825 to 2005
//
//------------------------------------------------------------------------------
void ETminUT (double T, double& DTsec, bool& valid);


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
double GMST (double MJD);


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
double Mjd ( int Year,   int Month, int Day, 
             int Hour=0, int Min=0, double Sec=0.0 );


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
              int& Year, int& Month, int& Day, double & Hour );


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
              int& Hour, int& Min, double& Sec );


//
// Auxiliary class for time output
//
class Time 
{
  public:
  
    // Constructors
    Time ();
    Time (double Hour, TimeFormat Format=HHMMSS); 
    
    // Modifiers
    void Set (TimeFormat Format=HHMMSS);

    // Time output
    friend std::ostream& operator << (std::ostream& os, const Time& time);

  private:

    double      m_Hour;
    TimeFormat  m_Format;
};


//
// Auxiliary class for date and time output
//
class DateTime 
{
  public:
  
    // Constructors
    DateTime ();
    DateTime (double Mjd, TimeFormat Format=None); 
    
    // Modifiers
    void Set (TimeFormat Format=HHMMSS);

    // Time output
    friend std::ostream& operator << (std::ostream& os, const DateTime& DT);

  private:

    double      m_Mjd;
    TimeFormat  m_Format;
};


#endif  // include blocker
