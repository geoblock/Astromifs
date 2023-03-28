//------------------------------------------------------------------------------
//
// APC_Planets.cpp
// 
//    Calculation of planetary positions
//
//------------------------------------------------------------------------------

#include <cmath>
#include <math.h>
#include <vcl.h>
#include <windows.h>    // Header file for windows
#include <stdio.h>      // Header file for standard Input/Output


#include "APC_Const.h"
#include "APC_Kepler.h"
#include "APC_Math.h"
#include "APC_Planets.h"
#include "APC_Sun.h"
#include "APC_PrecNut.h"
#include "APC_VecMat3D.h"


namespace // Unnamed namespace
{
  //
  // Constants
  //
  const int o   = 16;         // Index offset
  const int dim = 2*o+1;      // Work array dimension


  //
  // Definition of Pert class for summing up trigonometric perturbation 
  // series
  //
  class Pert
  {
    public:

      // Set time, mean anomalies and index range
      void Init ( double T, 
                  double M, int I_min, int I_max,
                  double m, int i_min, int i_max,
                  double phi=0.0 );
  
      // Sum-up perturbations in longitude, radius and latitude
      void Term ( int I, int i, int iT,
                  double dlc, double dls, 
                  double drc, double drs, 
                  double dbc, double dbs );

      // Retrieve perturbations in longitude, radius and latitude
      double dl();
      double dr();
      double db();

    private:
        
      double m_T;
      double m_cosM, m_sinM;
      double m_C[dim], m_S[dim], m_c[dim], m_s[dim];
      double m_dl, m_db, m_dr;
      double m_u, m_v;
  };


  // Set time, mean anomalies and index range
  void Pert::Init ( double T, 
                    double M, int I_min, int I_max,
                    double m, int i_min, int i_max,
                    double phi )
  {
    //
    // Variables
    //
    int i;


    m_dl=0.0; m_dr=0.0; m_db=0.0; // reset perturbations
    m_T=T;                        // set time

    
    // cosine and sine of multiples of M
    m_cosM=cos(M); m_sinM=sin(M);
  
    m_C[o]=cos(phi);
    m_S[o]=sin(phi);
  
    for (i=0; i<I_max; i++) 
      AddThe ( m_C[o+i],m_S[o+i], +m_cosM,+m_sinM, m_C[o+i+1],m_S[o+i+1] ); 
    for (i=0; i>I_min; i--) 
      AddThe ( m_C[o+i],m_S[o+i], +m_cosM,-m_sinM, m_C[o+i-1],m_S[o+i-1] ); 
  
    
    // cosine and sine of multiples of m
    m_c[o]=1.0; m_c[o+1]=cos(m); m_c[o-1]=+m_c[o+1];
    m_s[o]=0.0; m_s[o+1]=sin(m); m_s[o-1]=-m_s[o+1];
    
    for (i=1; i<i_max; i++) 
      AddThe ( m_c[o+i],m_s[o+i], m_c[o+1],m_s[o+1], m_c[o+i+1],m_s[o+i+1] );
    for (i=-1; i>i_min; i--) 
      AddThe ( m_c[o+i],m_s[o+i], m_c[o-1],m_s[o-1], m_c[o+i-1],m_s[o+i-1] );
  }
  

  // Sum-up perturbations in longitude, radius and latitude
  void Pert::Term ( int I, int i, int iT,
                    double dlc, double dls, 
                    double drc, double drs, 
                    double dbc, double dbs )
  {
    if (iT == 0) 
      AddThe ( m_C[o+I],m_S[o+I], m_c[o+i],m_s[o+i], m_u,m_v );
    else 
      { m_u *= m_T; m_v *= m_T; };

    m_dl += ( dlc*m_u + dls*m_v );
    m_dr += ( drc*m_u + drs*m_v );
    m_db += ( dbc*m_u + dbs*m_v );
  
  }


  // Retrieve perturbations in longitude, radius and latitude
  double Pert::dl() 
    { return m_dl; }

  double Pert::dr() 
    { return m_dr; }

  double Pert::db() 
    { return m_db; }

} // End of unnamed namespace


//------------------------------------------------------------------------------
//
// State: calculates position and velocity based on Keplerian orbits
//        (used by KepPosition() and KepVelocity())
//
// Input:
//
//   Planet  Identifies the planet
//   T       Time in Julian centuries since J2000
//
// Output:
//
//   r       Heliocentric position [AU], true-of-date ecliptic and equinox
//   v       Heliocentric velocity [AU/d], true-of-date ecliptic and equinox
//
//------------------------------------------------------------------------------
void State (PlanetType Planet, double T, Vec3D& r, Vec3D& v)
{
  //
  // Constants
  //
  const double p = 1.3970;        // Centential precession [deg]

  //
  // Variables
  //
  double a, e, M0, O, i, w, n, T0;
  Mat3D  PQR;


  // Orbital elements; ecliptic and equinox of J2000
  switch ( Planet )
  {
    case Sun:   
      r = v = Vec3D();  // Null vector
      return;
    case Mercury: 
      a =  0.387099; e = 0.205634; M0 = 174.7947; n = 149472.6738;
      O =  48.331;   i = 7.0048;   w  =  77.4552; T0 = 0.0;
      break;
    case Venus: 
      a =  0.723332; e = 0.006773; M0 =  50.4071; n  = 58517.8149;
      O =  76.680;   i = 3.3946;   w  = 131.5718; T0 = 0.0;
      break;
    case Earth:   
      a =  1.000000; e = 0.016709; M0 = 357.5256; n  = 35999.3720;
      O = 174.876;   i = 0.0000;   w  = 102.9400; T0 = 0.0;
      break;
    case Mars:  
      a =  1.523692; e = 0.093405; M0 =  19.3879; n  = 19140.3023;
      O =  49.557;   i = 1.8496;   w  = 336.0590; T0 = 0.0;
      break;
    case Jupiter: 
      a =  5.204267; e = 0.048775; M0 =  18.8185; n  =  3033.6272;
      O = 100.4908;  i = 1.3046;   w  =  15.5576; T0 = 0.0;
      break;
    case Saturn:  
      a =  9.582018; e = 0.055723; M0 = 320.3477; n  =  1213.8664;
      O = 113.6427;  i = 2.4852;   w  =  89.6567; T0 = 0.0;
      break;
    case Uranus: 
      a = 19.229412; e = 0.044406; M0 = 142.9559; n  =   426.9282;
      O =  73.9893;  i = 0.7726;   w  = 170.5310; T0 = 0.0;
      break;
    case Neptune:
      a = 30.103658; e = 0.011214; M0 = 267.7649; n  =   217.9599;
      O = 131.7942;  i = 1.7680;   w  =  37.4435; T0 = 0.0;
      break;
    case Pluto:
      a = 39.264230; e = 0.244672; M0 =  15.0233; n  =   146.3183;
      O = 110.2867;  i = 17.1514;  w =  224.0499; T0 = 0.0;
  }

  // State vector w.r.t. true-of-date ecliptic and equinox
  Ellip ( GM_Sun, Rad*(M0+n*(T-T0)), a, e, r, v );  // w.r.t. Orbital plane

  PQR = GaussVec ( Rad*(O+p*T), Rad*i, Rad*(w-O) ); // Transform. to ecliptic

  r = PQR*r;                                        // Ecliptic
  v = PQR*v;                                        // coordinates

}


//------------------------------------------------------------------------------
//
// KepPosition: Planetary position from Keplerian elements
//
// Input:
//
//   Planet    Identifies the planet
//   T         Time in Julian centuries since J2000
//
// <return>:   Heliocentric position [AU], true-of-date ecliptic and equinox
//
//------------------------------------------------------------------------------
Vec3D KepPosition (PlanetType Planet, double T)
{
  //
  // Variables
  //
  Vec3D r, v;

  
  State ( Planet, T, r, v );
  return r;
}


//------------------------------------------------------------------------------
//
// KepVelocity: Planetary velocity from Keplerian elements
//
// Input:
//
//   Planet    Identifies the planet
//   T         Time in Julian centuries since J2000
//
// <return>:   Heliocentric velocity [AU/d], true-of-date ecliptic and equinox
//
//------------------------------------------------------------------------------
Vec3D KepVelocity (PlanetType Planet, double T)
{
  //
  // Variables
  //
  Vec3D r, v;


  State (Planet, T, r, v);
  return v;
}


//------------------------------------------------------------------------------
//
// MercuryPos: Position of Mercury from analytical series
//
// Input:
//
//   T         Time in Julian centuries since J2000
//
// <return>:   Heliocentric position [AU], true-of-date ecliptic and equinox
//
//------------------------------------------------------------------------------
Vec3D MercuryPos (double T)
{
  //
  // Variables
  //
  double  M1,M2,M3,M5,M6;          // Mean anomalies
  Pert    Ven, Ear, Jup, Sat;      // Perturbations
  double  dl, dr, db;              // Corrections in longitude ["],
                                   // radius [AU] and latitude ["]
  double  l,b,r;                   // Ecliptic coordinates

  
  // Mean anomalies of planets in [rad]
  M1 = pi2 * Frac ( 0.4855407 + 415.2014314*T );
  M2 = pi2 * Frac ( 0.1394222 + 162.5490444*T );
  M3 = pi2 * Frac ( 0.9937861 +  99.9978139*T );
  M5 = pi2 * Frac ( 0.0558417 +   8.4298417*T );
  M6 = pi2 * Frac ( 0.8823333 +   3.3943333*T );

  
  // Keplerian motion and perturbations by Venus 
  Ven.Init ( T, M1,-1,9, M2,-5,0 );

  Ven.Term ( 1, 0,0, 259.74,84547.39,-78342.34, 0.01,11683.22,21203.79);
  Ven.Term ( 1, 0,1,   2.30,    5.04,    -7.52, 0.02,  138.55,  -71.01);
  Ven.Term ( 1, 0,2,   0.01,   -0.01,     0.01, 0.01,   -0.19,   -0.54);
  Ven.Term ( 2, 0,0,-549.71,10394.44, -7955.45, 0.00, 2390.29, 4306.79);
  Ven.Term ( 2, 0,1,  -4.77,    8.97,    -1.53, 0.00,   28.49,  -14.18);
  Ven.Term ( 2, 0,2,   0.00,    0.00,     0.00, 0.00,   -0.04,   -0.11);
  Ven.Term ( 3, 0,0,-234.04, 1748.74, -1212.86, 0.00,  535.41,  984.33);
  Ven.Term ( 3, 0,1,  -2.03,    3.48,    -0.35, 0.00,    6.56,   -2.91);
  Ven.Term ( 4, 0,0, -77.64,  332.63,  -219.23, 0.00,  124.40,  237.03);
  Ven.Term ( 4, 0,1,  -0.70,    1.10,    -0.08, 0.00,    1.59,   -0.59);
  Ven.Term ( 5, 0,0, -23.59,   67.28,   -43.54, 0.00,   29.44,   58.77);
  Ven.Term ( 5, 0,1,  -0.23,    0.32,    -0.02, 0.00,    0.39,   -0.11);
  Ven.Term ( 6, 0,0,  -6.86,   14.06,    -9.18, 0.00,    7.03,   14.84);
  Ven.Term ( 6, 0,1,  -0.07,    0.09,    -0.01, 0.00,    0.10,   -0.02);
  Ven.Term ( 7, 0,0,  -1.94,    2.98,    -2.02, 0.00,    1.69,    3.80);
  Ven.Term ( 8, 0,0,  -0.54,    0.63,    -0.46, 0.00,    0.41,    0.98);
  Ven.Term ( 9, 0,0,  -0.15,    0.13,    -0.11, 0.00,    0.10,    0.25);
  Ven.Term (-1,-2,0,  -0.17,   -0.06,    -0.05, 0.14,   -0.06,   -0.07);
  Ven.Term ( 0,-1,0,   0.24,   -0.16,    -0.11,-0.16,    0.04,   -0.01);
  Ven.Term ( 0,-2,0,  -0.68,   -0.25,    -0.26, 0.73,   -0.16,   -0.18);
  Ven.Term ( 0,-5,0,   0.37,    0.08,     0.06,-0.28,    0.13,    0.12);
  Ven.Term ( 1,-1,0,   0.58,   -0.41,     0.26, 0.36,    0.01,   -0.01);
  Ven.Term ( 1,-2,0,  -3.51,   -1.23,     0.23,-0.63,   -0.05,   -0.06);
  Ven.Term ( 1,-3,0,   0.08,    0.53,    -0.11, 0.04,    0.02,   -0.09);
  Ven.Term ( 1,-5,0,   1.44,    0.31,     0.30,-1.39,    0.34,    0.29);
  Ven.Term ( 2,-1,0,   0.15,   -0.11,     0.09, 0.12,    0.02,   -0.04);
  Ven.Term ( 2,-2,0,  -1.99,   -0.68,     0.65,-1.91,   -0.20,    0.03);
  Ven.Term ( 2,-3,0,  -0.34,   -1.28,     0.97,-0.26,    0.03,    0.03);
  Ven.Term ( 2,-4,0,  -0.33,    0.35,    -0.13,-0.13,   -0.01,    0.00);
  Ven.Term ( 2,-5,0,   7.19,    1.56,    -0.05, 0.12,    0.06,    0.05);
  Ven.Term ( 3,-2,0,  -0.52,   -0.18,     0.13,-0.39,   -0.16,    0.03);
  Ven.Term ( 3,-3,0,  -0.11,   -0.42,     0.36,-0.10,   -0.05,   -0.05);
  Ven.Term ( 3,-4,0,  -0.19,    0.22,    -0.23,-0.20,   -0.01,    0.02);
  Ven.Term ( 3,-5,0,   2.77,    0.49,    -0.45, 2.56,    0.40,   -0.12);
  Ven.Term ( 4,-5,0,   0.67,    0.12,    -0.09, 0.47,    0.24,   -0.08);
  Ven.Term ( 5,-5,0,   0.18,    0.03,    -0.02, 0.12,    0.09,   -0.03);

  dl = Ven.dl();  dr = Ven.dr();  db = Ven.db();


  // Perturbations by the Earth
  Ear.Init ( T, M1,0,2, M3,-4,-1 );

  Ear.Term ( 0,-4,0,  -0.11,   -0.07,    -0.08, 0.11,   -0.02,   -0.04);
  Ear.Term ( 1,-1,0,   0.10,   -0.20,     0.15, 0.07,    0.00,    0.00);
  Ear.Term ( 1,-2,0,  -0.35,    0.28,    -0.13,-0.17,   -0.01,    0.00);
  Ear.Term ( 1,-4,0,  -0.67,   -0.45,     0.00, 0.01,   -0.01,   -0.01);
  Ear.Term ( 2,-2,0,  -0.20,    0.16,    -0.16,-0.20,   -0.01,    0.02);
  Ear.Term ( 2,-3,0,   0.13,   -0.02,     0.02, 0.14,    0.01,    0.00);
  Ear.Term ( 2,-4,0,  -0.33,   -0.18,     0.17,-0.31,   -0.04,    0.00);

  dl += Ear.dl();  dr += Ear.dr();  db += Ear.db();

  
  // Perturbations by Jupiter 
  Jup.Init ( T, M1,-1,3, M5,-3,-1 );

  Jup.Term (-1,-1,0,  -0.08,    0.16,     0.15, 0.08,   -0.04,    0.01);
  Jup.Term (-1,-2,0,   0.10,   -0.06,    -0.07,-0.12,    0.07,   -0.01);
  Jup.Term ( 0,-1,0,  -0.31,    0.48,    -0.02, 0.13,   -0.03,   -0.02);
  Jup.Term ( 0,-2,0,   0.42,   -0.26,    -0.38,-0.50,    0.20,   -0.03);
  Jup.Term ( 1,-1,0,  -0.70,    0.01,    -0.02,-0.63,    0.00,    0.03);
  Jup.Term ( 1,-2,0,   2.61,   -1.97,     1.74, 2.32,    0.01,    0.01);
  Jup.Term ( 1,-3,0,   0.32,   -0.15,     0.13, 0.28,    0.00,    0.00);
  Jup.Term ( 2,-1,0,  -0.18,    0.01,     0.00,-0.13,   -0.03,    0.03);
  Jup.Term ( 2,-2,0,   0.75,   -0.56,     0.45, 0.60,    0.08,   -0.17);
  Jup.Term ( 3,-2,0,   0.20,   -0.15,     0.10, 0.14,    0.04,   -0.08);

  dl += Jup.dl();  dr += Jup.dr();  db += Jup.db();


  // Perturbations by Saturn 
  Sat.Init ( T, M1,1,1, M6,-2,-2 );

  Sat.Term ( 1,-2,0,  -0.19,    0.33,     0.00, 0.00,    0.00,    0.00);

  dl += Sat.dl();  dr += Sat.dr();  db += Sat.db();

  
  // Ecliptic coordinates ([rad],[AU])
  dl += +2.8 + 3.2*T;

  l = pi2 * Frac (0.2151379 + M1/pi2 + ( (5601.7+1.1*T)*T + dl ) / 1296.0e3);
  r = 0.3952829 + 0.0000016 * T  +  dr * 1.0e-6;
  b = ( -2522.15 + (-30.18 + 0.04*T) * T  +  db ) / Arcs;

  return Vec3D ( Polar(l, b, r) );  // Position vector

}


//------------------------------------------------------------------------------
//
// VenusPos: Position of Venus from analytical series
//
// Input:
//
//   T         Time in Julian centuries since J2000
//
// <return>:   Heliocentric position [AU], true-of-date ecliptic and equinox
//
//------------------------------------------------------------------------------
Vec3D VenusPos (double T)
{
  //
  // Variables
  //
  double  M1,M2,M3,M4,M5,M6;       // Mean anomalies
  Pert    Mer,Ear,Mar,Jup,Sat;     // Perturbations
  double  dl, dr, db;              // Corrections in longitude ["],
                                   // radius [AU] and latitude ["]
  double  l,b,r;                   // Ecliptic coordinates


  // Mean anomalies of planets in [rad]
  M1 = pi2 * Frac ( 0.4861431 + 415.2018375*T );
  M2 = pi2 * Frac ( 0.1400197 + 162.5494552*T );
  M3 = pi2 * Frac ( 0.9944153 +  99.9982208*T );
  M4 = pi2 * Frac ( 0.0556297 +  53.1674631*T );
  M5 = pi2 * Frac ( 0.0567028 +   8.4305083*T );
  M6 = pi2 * Frac ( 0.8830539 +   3.3947206*T );


  // Perturbations by Mercury
  Mer.Init ( T, M2,1,5, M1,-2,-1 );

  Mer.Term (1,-1,0,   0.00,   0.00,    0.06, -0.09,   0.01,   0.00);
  Mer.Term (2,-1,0,   0.25,  -0.09,   -0.09, -0.27,   0.00,   0.00);
  Mer.Term (4,-2,0,  -0.07,  -0.08,   -0.14,  0.14,  -0.01,  -0.01);
  Mer.Term (5,-2,0,  -0.35,   0.08,    0.02,  0.09,   0.00,   0.00);
  
  dl = Mer.dl();  dr = Mer.dr();  db = Mer.db();


  // Keplerian motion and perturbations by the Earth 
  Ear.Init ( T, M2,1,8, M3,-8,0 );

  Ear.Term (1, 0,0,   2.37,2793.23,-4899.07,  0.11,9995.27,7027.22);
  Ear.Term (1, 0,1,   0.10, -19.65,   34.40,  0.22,  64.95, -86.10);
  Ear.Term (1, 0,2,   0.06,   0.04,   -0.07,  0.11,  -0.55,  -0.07);
  Ear.Term (2, 0,0,-170.42,  73.13,  -16.59,  0.00,  67.71,  47.56);
  Ear.Term (2, 0,1,   0.93,   2.91,    0.23,  0.00,  -0.03,  -0.92);
  Ear.Term (3, 0,0,  -2.31,   0.90,   -0.08,  0.00,   0.04,   2.09);
  Ear.Term (1,-1,0,  -2.38,  -4.27,    3.27, -1.82,   0.00,   0.00);
  Ear.Term (1,-2,0,   0.09,   0.00,   -0.08,  0.05,  -0.02,  -0.25);
  Ear.Term (2,-2,0,  -9.57,  -5.93,    8.57,-13.83,  -0.01,  -0.01);
  Ear.Term (2,-3,0,  -2.47,  -2.40,    0.83, -0.95,   0.16,   0.24);
  Ear.Term (3,-2,0,  -0.09,  -0.05,    0.08, -0.13,  -0.28,   0.12);
  Ear.Term (3,-3,0,   7.12,   0.32,   -0.62, 13.76,  -0.07,   0.01);
  Ear.Term (3,-4,0,  -0.65,  -0.17,    0.18, -0.73,   0.10,   0.05);
  Ear.Term (3,-5,0,  -1.08,  -0.95,   -0.17,  0.22,  -0.03,  -0.03);
  Ear.Term (4,-3,0,   0.06,   0.00,   -0.01,  0.08,   0.14,  -0.18);
  Ear.Term (4,-4,0,   0.93,  -0.46,    1.06,  2.13,  -0.01,   0.01);
  Ear.Term (4,-5,0,  -1.53,   0.38,   -0.64, -2.54,   0.27,   0.00);
  Ear.Term (4,-6,0,  -0.17,  -0.05,    0.03, -0.11,   0.02,   0.00);
  Ear.Term (5,-5,0,   0.18,  -0.28,    0.71,  0.47,  -0.02,   0.04);
  Ear.Term (5,-6,0,   0.15,  -0.14,    0.30,  0.31,  -0.04,   0.03);
  Ear.Term (5,-7,0,  -0.08,   0.02,   -0.03, -0.11,   0.01,   0.00);
  Ear.Term (5,-8,0,  -0.23,   0.00,    0.01, -0.04,   0.00,   0.00);
  Ear.Term (6,-6,0,   0.01,  -0.14,    0.39,  0.04,   0.00,  -0.01);
  Ear.Term (6,-7,0,   0.02,  -0.05,    0.12,  0.04,  -0.01,   0.01);
  Ear.Term (6,-8,0,   0.10,  -0.10,    0.19,  0.19,  -0.02,   0.02);
  Ear.Term (7,-7,0,  -0.03,  -0.06,    0.18, -0.08,   0.00,   0.00);
  Ear.Term (8,-8,0,  -0.03,  -0.02,    0.06, -0.08,   0.00,   0.00);

  dl += Ear.dl();  dr += Ear.dr();  db += Ear.db();


  // Perturbations by Mars
  Mar.Init ( T, M2,1,2, M4,-3,-2 );

  Mar.Term (1,-3,0,  -0.65,   1.02,   -0.04, -0.02,  -0.02,   0.00);
  Mar.Term (2,-2,0,  -0.05,   0.04,   -0.09, -0.10,   0.00,   0.00);
  Mar.Term (2,-3,0,  -0.50,   0.45,   -0.79, -0.89,   0.01,   0.03);
  
  dl += Mar.dl();  dr += Mar.dr();  db += Mar.db();


  // Perturbations by Jupiter 
  Jup.Init ( T, M2,0,3, M5,-3,-1 );

  Jup.Term (0,-1,0,  -0.05,   1.56,    0.16,  0.04,  -0.08,  -0.04);
  Jup.Term (1,-1,0,  -2.62,   1.40,   -2.35, -4.40,   0.02,   0.03);
  Jup.Term (1,-2,0,  -0.47,  -0.08,    0.12, -0.76,   0.04,  -0.18);
  Jup.Term (2,-2,0,  -0.73,  -0.51,    1.27, -1.82,  -0.01,   0.01);
  Jup.Term (2,-3,0,  -0.14,  -0.10,    0.25, -0.34,   0.00,   0.00);
  Jup.Term (3,-3,0,  -0.01,   0.04,   -0.11, -0.02,   0.00,   0.00);

  dl += Jup.dl(); dr += Jup.dr();  db += Jup.db();


  // Perturbations by Saturn 
  Sat.Init ( T, M2,0,1, M6,-1,-1 );

  Sat.Term (0,-1,0,   0.00,   0.21,    0.00,  0.00,   0.00,  -0.01);
  Sat.Term (1,-1,0,  -0.11,  -0.14,    0.24, -0.20,   0.01,   0.00);

  dl += Sat.dl(); dr += Sat.dr();  db += Sat.db();


  // Ecliptic coordinates ([rad],[AU])
  dl +=  + 2.74*sin(pi2*(0.0764+0.4174*T))
         + 0.27*sin(pi2*(0.9201+0.3307*T));
  dl +=  +1.9 + 1.8*T;
  
  l = pi2 * Frac (0.3654783 + M2/pi2 + ( (5071.2+1.1*T)*T + dl ) / 1296.0e3);
  r = 0.7233482 - 0.0000002 * T  +  dr * 1.0e-6;
  b = ( -67.70 + ( 0.04 + 0.01*T) * T  +  db ) / Arcs;

  return Vec3D ( Polar(l,b,r) );  // Position vector

}


//------------------------------------------------------------------------------
//
// MarsPos: Position of Mars from analytical series
//
// Input:
//
//   T         Time in Julian centuries since J2000
//
// <return>:   Heliocentric position [AU], true-of-date ecliptic and equinox
//
//------------------------------------------------------------------------------
Vec3D MarsPos (double T)
{
  //
  // Variables
  //
  double  M2,M3,M4,M5,M6;          // Mean anomalies
  Pert    Ven,Ear,Jup,Sat;         // Perturbations
  double  dl, dr, db;              // Corrections in longitude ["],
                                   // radius [AU] and latitude ["]
  double  l,b,r;                   // Ecliptic coordinates


  // Mean anomalies of planets in [rad]
  M2 = pi2 * Frac ( 0.1382208 + 162.5482542*T );
  M3 = pi2 * Frac ( 0.9926208 +  99.9970236*T );
  M4 = pi2 * Frac ( 0.0538553 +  53.1662736*T );
  M5 = pi2 * Frac ( 0.0548944 +   8.4290611*T );
  M6 = pi2 * Frac ( 0.8811167 +   3.3935250*T );


  // Perturbations by Venus
  Ven.Init ( T, M4,0,7, M2,-2,-1 );

  Ven.Term ( 0,-1,0, -0.01,   -0.03,      0.10, -0.04,    0.00,   0.00);
  Ven.Term ( 1,-1,0,  0.05,    0.10,     -2.08,  0.75,    0.00,   0.00);
  Ven.Term ( 2,-1,0, -0.25,   -0.57,     -2.58,  1.18,    0.05,  -0.04);
  Ven.Term ( 2,-2,0,  0.02,    0.02,      0.13, -0.14,    0.00,   0.00);
  Ven.Term ( 3,-1,0,  3.41,    5.38,      1.87, -1.15,    0.01,  -0.01);
  Ven.Term ( 3,-2,0,  0.02,    0.02,      0.11, -0.13,    0.00,   0.00);
  Ven.Term ( 4,-1,0,  0.32,    0.49,     -1.88,  1.21,   -0.07,   0.07);
  Ven.Term ( 4,-2,0,  0.03,    0.03,      0.12, -0.14,    0.00,   0.00);
  Ven.Term ( 5,-1,0,  0.04,    0.06,     -0.17,  0.11,   -0.01,   0.01);
  Ven.Term ( 5,-2,0,  0.11,    0.09,      0.35, -0.43,   -0.01,   0.01);
  Ven.Term ( 6,-2,0, -0.36,   -0.28,     -0.20,  0.25,    0.00,   0.00);
  Ven.Term ( 7,-2,0, -0.03,   -0.03,      0.11, -0.13,    0.00,  -0.01);
  
  dl = Ven.dl();  dr = Ven.dr();  db = Ven.db();


  // Keplerian motion and perturbations by the Earth 
  Ear.Init ( T, M4,-1,16, M3,-9,0 );

  Ear.Term ( 1, 0,0, -5.32,38481.97,-141856.04,  0.40,-6321.67,1876.89);
  Ear.Term ( 1, 0,1, -1.12,   37.98,   -138.67, -2.93,   37.28, 117.48);
  Ear.Term ( 1, 0,2, -0.32,   -0.03,      0.12, -1.19,    1.04,  -0.40);
  Ear.Term ( 2, 0,0, 28.28, 2285.80,  -6608.37,  0.00, -589.35, 174.81);
  Ear.Term ( 2, 0,1,  1.64,    3.37,    -12.93,  0.00,    2.89,  11.10);
  Ear.Term ( 2, 0,2,  0.00,    0.00,      0.00,  0.00,    0.10,  -0.03);
  Ear.Term ( 3, 0,0,  5.31,  189.29,   -461.81,  0.00,  -61.98,  18.53);
  Ear.Term ( 3, 0,1,  0.31,    0.35,     -1.36,  0.00,    0.25,   1.19);
  Ear.Term ( 4, 0,0,  0.81,   17.96,    -38.26,  0.00,   -6.88,   2.08);
  Ear.Term ( 4, 0,1,  0.05,    0.04,     -0.15,  0.00,    0.02,   0.14);
  Ear.Term ( 5, 0,0,  0.11,    1.83,     -3.48,  0.00,   -0.79,   0.24);
  Ear.Term ( 6, 0,0,  0.02,    0.20,     -0.34,  0.00,   -0.09,   0.03);
  Ear.Term (-1,-1,0,  0.09,    0.06,      0.14, -0.22,    0.02,  -0.02);
  Ear.Term ( 0,-1,0,  0.72,    0.49,      1.55, -2.31,    0.12,  -0.10);
  Ear.Term ( 1,-1,0,  7.00,    4.92,     13.93,-20.48,    0.08,  -0.13);
  Ear.Term ( 2,-1,0, 13.08,    4.89,     -4.53, 10.01,   -0.05,   0.13);
  Ear.Term ( 2,-2,0,  0.14,    0.05,     -0.48, -2.66,    0.01,   0.14);
  Ear.Term ( 3,-1,0,  1.38,    0.56,     -2.00,  4.85,   -0.01,   0.19);
  Ear.Term ( 3,-2,0, -6.85,    2.68,      8.38, 21.42,    0.00,   0.03);
  Ear.Term ( 3,-3,0, -0.08,    0.20,      1.20,  0.46,    0.00,   0.00);
  Ear.Term ( 4,-1,0,  0.16,    0.07,     -0.19,  0.47,   -0.01,   0.05);
  Ear.Term ( 4,-2,0, -4.41,    2.14,     -3.33, -7.21,   -0.07,  -0.09);
  Ear.Term ( 4,-3,0, -0.12,    0.33,      2.22,  0.72,   -0.03,  -0.02);
  Ear.Term ( 4,-4,0, -0.04,   -0.06,     -0.36,  0.23,    0.00,   0.00);
  Ear.Term ( 5,-2,0, -0.44,    0.21,     -0.70, -1.46,   -0.06,  -0.07);
  Ear.Term ( 5,-3,0,  0.48,   -2.60,     -7.25, -1.37,    0.00,   0.00);
  Ear.Term ( 5,-4,0, -0.09,   -0.12,     -0.66,  0.50,    0.00,   0.00);
  Ear.Term ( 5,-5,0,  0.03,    0.00,      0.01, -0.17,    0.00,   0.00);
  Ear.Term ( 6,-2,0, -0.05,    0.03,     -0.07, -0.15,   -0.01,  -0.01);
  Ear.Term ( 6,-3,0,  0.10,   -0.96,      2.36,  0.30,    0.04,   0.00);
  Ear.Term ( 6,-4,0, -0.17,   -0.20,     -1.09,  0.94,    0.02,  -0.02);
  Ear.Term ( 6,-5,0,  0.05,    0.00,      0.00, -0.30,    0.00,   0.00);
  Ear.Term ( 7,-3,0,  0.01,   -0.10,      0.32,  0.04,    0.02,   0.00);
  Ear.Term ( 7,-4,0,  0.86,    0.77,      1.86, -2.01,    0.01,  -0.01);
  Ear.Term ( 7,-5,0,  0.09,   -0.01,     -0.05, -0.44,    0.00,   0.00);
  Ear.Term ( 7,-6,0, -0.01,    0.02,      0.10,  0.08,    0.00,   0.00);
  Ear.Term ( 8,-4,0,  0.20,    0.16,     -0.53,  0.64,   -0.01,   0.02);
  Ear.Term ( 8,-5,0,  0.17,   -0.03,     -0.14, -0.84,    0.00,   0.01);
  Ear.Term ( 8,-6,0, -0.02,    0.03,      0.16,  0.09,    0.00,   0.00);
  Ear.Term ( 9,-5,0, -0.55,    0.15,      0.30,  1.10,    0.00,   0.00);
  Ear.Term ( 9,-6,0, -0.02,    0.04,      0.20,  0.10,    0.00,   0.00);
  Ear.Term (10,-5,0, -0.09,    0.03,     -0.10, -0.33,    0.00,  -0.01);
  Ear.Term (10,-6,0, -0.05,    0.11,      0.48,  0.21,   -0.01,   0.00);
  Ear.Term (11,-6,0,  0.10,   -0.35,     -0.52, -0.15,    0.00,   0.00);
  Ear.Term (11,-7,0, -0.01,   -0.02,     -0.10,  0.07,    0.00,   0.00);
  Ear.Term (12,-6,0,  0.01,   -0.04,      0.18,  0.04,    0.01,   0.00);
  Ear.Term (12,-7,0, -0.05,   -0.07,     -0.29,  0.20,    0.01,   0.00);
  Ear.Term (13,-7,0,  0.23,    0.27,      0.25, -0.21,    0.00,   0.00);
  Ear.Term (14,-7,0,  0.02,    0.03,     -0.10,  0.09,    0.00,   0.00);
  Ear.Term (14,-8,0,  0.05,    0.01,      0.03, -0.23,    0.00,   0.03);
  Ear.Term (15,-8,0, -1.53,    0.27,      0.06,  0.42,    0.00,   0.00);
  Ear.Term (16,-8,0, -0.14,    0.02,     -0.10, -0.55,   -0.01,  -0.02);
  Ear.Term (16,-9,0,  0.03,   -0.06,     -0.25, -0.11,    0.00,   0.00);

  dl += Ear.dl();  dr += Ear.dr();  db += Ear.db();


  // Perturbations by Jupiter 
  Jup.Init ( T, M4,-2,5, M5,-5,-1 );

  Jup.Term (-2,-1,0,  0.05,    0.03,      0.08, -0.14,    0.01,  -0.01);
  Jup.Term (-1,-1,0,  0.39,    0.27,      0.92, -1.50,   -0.03,  -0.06);
  Jup.Term (-1,-2,0, -0.16,    0.03,      0.13,  0.67,   -0.01,   0.06);
  Jup.Term (-1,-3,0, -0.02,    0.01,      0.05,  0.09,    0.00,   0.01);
  Jup.Term ( 0,-1,0,  3.56,    1.13,     -5.41, -7.18,   -0.25,  -0.24);
  Jup.Term ( 0,-2,0, -1.44,    0.25,      1.24,  7.96,    0.02,   0.31);
  Jup.Term ( 0,-3,0, -0.21,    0.11,      0.55,  1.04,    0.01,   0.05);
  Jup.Term ( 0,-4,0, -0.02,    0.02,      0.11,  0.11,    0.00,   0.01);
  Jup.Term ( 1,-1,0, 16.67,  -19.15,     61.00, 53.36,   -0.06,  -0.07);
  Jup.Term ( 1,-2,0,-21.64,    3.18,     -7.77,-54.64,   -0.31,   0.50);
  Jup.Term ( 1,-3,0, -2.82,    1.45,     -2.53, -5.73,    0.01,   0.07);
  Jup.Term ( 1,-4,0, -0.31,    0.28,     -0.34, -0.51,    0.00,   0.00);
  Jup.Term ( 2,-1,0,  2.15,   -2.29,      7.04,  6.94,    0.33,   0.19);
  Jup.Term ( 2,-2,0,-15.69,    3.31,    -15.70,-73.17,   -0.17,  -0.25);
  Jup.Term ( 2,-3,0, -1.73,    1.95,     -9.19, -7.20,    0.02,  -0.03);
  Jup.Term ( 2,-4,0, -0.01,    0.33,     -1.42,  0.08,    0.01,  -0.01);
  Jup.Term ( 2,-5,0,  0.03,    0.03,     -0.13,  0.12,    0.00,   0.00);
  Jup.Term ( 3,-1,0,  0.26,   -0.28,      0.73,  0.71,    0.08,   0.04);
  Jup.Term ( 3,-2,0, -2.06,    0.46,     -1.61, -6.72,   -0.13,  -0.25);
  Jup.Term ( 3,-3,0, -1.28,   -0.27,      2.21, -6.90,   -0.04,  -0.02);
  Jup.Term ( 3,-4,0, -0.22,    0.08,     -0.44, -1.25,    0.00,   0.01);
  Jup.Term ( 3,-5,0, -0.02,    0.03,     -0.15, -0.08,    0.00,   0.00);
  Jup.Term ( 4,-1,0,  0.03,   -0.03,      0.08,  0.08,    0.01,   0.01);
  Jup.Term ( 4,-2,0, -0.26,    0.06,     -0.17, -0.70,   -0.03,  -0.05);
  Jup.Term ( 4,-3,0, -0.20,   -0.05,      0.22, -0.79,   -0.01,  -0.02);
  Jup.Term ( 4,-4,0, -0.11,   -0.14,      0.93, -0.60,    0.00,   0.00);
  Jup.Term ( 4,-5,0, -0.04,   -0.02,      0.09, -0.23,    0.00,   0.00);
  Jup.Term ( 5,-4,0, -0.02,   -0.03,      0.13, -0.09,    0.00,   0.00);
  Jup.Term ( 5,-5,0,  0.00,   -0.03,      0.21,  0.01,    0.00,   0.00);

  dl += Jup.dl();  dr += Jup.dr();  db += Jup.db();


  // Perturbations by Saturn 
  Sat.Init ( T, M4,-1,3, M6,-4,-1 );

  Sat.Term (-1,-1,0,  0.03,    0.13,      0.48, -0.13,    0.02,   0.00);
  Sat.Term ( 0,-1,0,  0.27,    0.84,      0.40, -0.43,    0.01,  -0.01);
  Sat.Term ( 0,-2,0,  0.12,   -0.04,     -0.33, -0.55,   -0.01,  -0.02);
  Sat.Term ( 0,-3,0,  0.02,   -0.01,     -0.07, -0.08,    0.00,   0.00);
  Sat.Term ( 1,-1,0,  1.12,    0.76,     -2.66,  3.91,   -0.01,   0.01);
  Sat.Term ( 1,-2,0,  1.49,   -0.95,      3.07,  4.83,    0.04,  -0.05);
  Sat.Term ( 1,-3,0,  0.21,   -0.18,      0.55,  0.64,    0.00,   0.00);
  Sat.Term ( 2,-1,0,  0.12,    0.10,     -0.29,  0.34,   -0.01,   0.02);
  Sat.Term ( 2,-2,0,  0.51,   -0.36,      1.61,  2.25,    0.03,   0.01);
  Sat.Term ( 2,-3,0,  0.10,   -0.10,      0.50,  0.43,    0.00,   0.00);
  Sat.Term ( 2,-4,0,  0.01,   -0.02,      0.11,  0.05,    0.00,   0.00);
  Sat.Term ( 3,-2,0,  0.07,   -0.05,      0.16,  0.22,    0.01,   0.01);

  dl += Sat.dl();  dr += Sat.dr();  db += Sat.db();


  // Ecliptic coordinates ([rad],[AU])
  dl += + 52.49*sin(pi2*(0.1868+0.0549*T))
        +  0.61*sin(pi2*(0.9220+0.3307*T))
        +  0.32*sin(pi2*(0.4731+2.1485*T))
        +  0.28*sin(pi2*(0.9467+0.1133*T));
  dl += + 0.14 + 0.87*T - 0.11*T*T;

  l = pi2 * Frac (0.9334591 + M4/pi2 + ( (6615.5+1.1*T)*T + dl ) / 1296.0e3);
  r = 1.5303352 + 0.0000131*T  +  dr * 1.0E-6;
  b = ( 596.32 + (-2.92 - 0.10*T) * T  +  db ) / Arcs;

  return Vec3D ( Polar(l, b, r) );  // Position vector

}


//------------------------------------------------------------------------------
//
// JupiterPos: Position of Jupiter from analytical series
//
// Input:
//
//   T         Time in Julian centuries since J2000
//
// <return>:   Heliocentric position [AU], true-of-date ecliptic and equinox
//
//------------------------------------------------------------------------------
Vec3D JupiterPos (double T)
{
  //
  // Variables
  //
  double  M5,M6,M7;                // Mean anomalies
  Pert    Sat,Ura;                 // Perturbations
  double  phi,c,s;
  double  dl, dr, db;              // Corrections in longitude ["],
                                   // radius [AU] and latitude ["]
  double  l,b,r;                   // Ecliptic coordinates


  // Mean anomalies of planets in [rad]
  M5 = pi2 * Frac ( 0.0565314 + 8.4302963*T );
  M6 = pi2 * Frac ( 0.8829867 + 3.3947688*T );
  M7 = pi2 * Frac ( 0.3969537 + 1.1902586*T );


  // Keplerian motion and perturbations by Saturn
  Sat.Init ( T, M5,-1,5, M6,-10,0 );

  Sat.Term (-1, -1,0,  -0.2,    1.4,     2.0,   0.6,    0.1, -0.2);
  Sat.Term ( 0, -1,0,   9.4,    8.9,     3.9,  -8.3,   -0.4, -1.4);
  Sat.Term ( 0, -2,0,   5.6,   -3.0,    -5.4,  -5.7,   -2.0,  0.0);
  Sat.Term ( 0, -3,0,  -4.0,   -0.1,     0.0,   5.5,    0.0,  0.0);
  Sat.Term ( 0, -5,0,   3.3,   -1.6,    -1.6,  -3.1,   -0.5, -1.2);
  Sat.Term ( 1,  0,0,-113.1,19998.6,-25208.2,-142.2,-4670.7,288.9);
  Sat.Term ( 1,  0,1, -76.1,   66.9,   -84.2, -95.8,   21.6, 29.4);
  Sat.Term ( 1,  0,2,  -0.5,   -0.3,     0.4,  -0.7,    0.1, -0.1);
  Sat.Term ( 1, -1,0,  78.8,  -14.5,    11.5,  64.4,   -0.2,  0.2);
  Sat.Term ( 1, -2,0,  -2.0, -132.4,    28.8,   4.3,   -1.7,  0.4);
  Sat.Term ( 1, -2,1,  -1.1,   -0.7,     0.2,  -0.3,    0.0,  0.0);
  Sat.Term ( 1, -3,0,  -7.5,   -6.8,    -0.4,  -1.1,    0.6, -0.9);
  Sat.Term ( 1, -4,0,   0.7,    0.7,     0.6,  -1.1,    0.0, -0.2);
  Sat.Term ( 1, -5,0,  51.5,  -26.0,   -32.5, -64.4,   -4.9,-12.4);
  Sat.Term ( 1, -5,1,  -1.2,   -2.2,    -2.7,   1.5,   -0.4,  0.3);
  Sat.Term ( 2,  0,0,  -3.4,  632.0,  -610.6,  -6.5, -226.8, 12.7);
  Sat.Term ( 2,  0,1,  -4.2,    3.8,    -4.1,  -4.5,    0.2,  0.6);
  Sat.Term ( 2, -1,0,   5.3,   -0.7,     0.7,   6.1,    0.2,  1.1);
  Sat.Term ( 2, -2,0, -76.4, -185.1,   260.2,-108.0,    1.6,  0.0);
  Sat.Term ( 2, -3,0,  66.7,   47.8,   -51.4,  69.8,    0.9,  0.3);
  Sat.Term ( 2, -3,1,   0.6,   -1.0,     1.0,   0.6,    0.0,  0.0);
  Sat.Term ( 2, -4,0,  17.0,    1.4,    -1.8,   9.6,    0.0, -0.1);
  Sat.Term ( 2, -5,0,1066.2, -518.3,    -1.3, -23.9,    1.8, -0.3);
  Sat.Term ( 2, -5,1, -25.4,  -40.3,    -0.9,   0.3,    0.0,  0.0);
  Sat.Term ( 2, -5,2,  -0.7,    0.5,     0.0,   0.0,    0.0,  0.0);
  Sat.Term ( 3,  0,0,  -0.1,   28.0,   -22.1,  -0.2,  -12.5,  0.7);
  Sat.Term ( 3, -2,0,  -5.0,  -11.5,    11.7,  -5.4,    2.1, -1.0);
  Sat.Term ( 3, -3,0,  16.9,   -6.4,    13.4,  26.9,   -0.5,  0.8);
  Sat.Term ( 3, -4,0,   7.2,  -13.3,    20.9,  10.5,    0.1, -0.1);
  Sat.Term ( 3, -5,0,  68.5,  134.3,  -166.9,  86.5,    7.1, 15.2);
  Sat.Term ( 3, -5,1,   3.5,   -2.7,     3.4,   4.3,    0.5, -0.4);
  Sat.Term ( 3, -6,0,   0.6,    1.0,    -0.9,   0.5,    0.0,  0.0);
  Sat.Term ( 3, -7,0,  -1.1,    1.7,    -0.4,  -0.2,    0.0,  0.0);
  Sat.Term ( 4,  0,0,   0.0,    1.4,    -1.0,   0.0,   -0.6,  0.0);
  Sat.Term ( 4, -2,0,  -0.3,   -0.7,     0.4,  -0.2,    0.2, -0.1);
  Sat.Term ( 4, -3,0,   1.1,   -0.6,     0.9,   1.2,    0.1,  0.2);
  Sat.Term ( 4, -4,0,   3.2,    1.7,    -4.1,   5.8,    0.2,  0.1);
  Sat.Term ( 4, -5,0,   6.7,    8.7,    -9.3,   8.7,   -1.1,  1.6);
  Sat.Term ( 4, -6,0,   1.5,   -0.3,     0.6,   2.4,    0.0,  0.0);
  Sat.Term ( 4, -7,0,  -1.9,    2.3,    -3.2,  -2.7,    0.0, -0.1);
  Sat.Term ( 4, -8,0,   0.4,   -1.8,     1.9,   0.5,    0.0,  0.0);
  Sat.Term ( 4, -9,0,  -0.2,   -0.5,     0.3,  -0.1,    0.0,  0.0);
  Sat.Term ( 4,-10,0,  -8.6,   -6.8,    -0.4,   0.1,    0.0,  0.0);
  Sat.Term ( 4,-10,1,  -0.5,    0.6,     0.0,   0.0,    0.0,  0.0);
  Sat.Term ( 5, -5,0,  -0.1,    1.5,    -2.5,  -0.8,   -0.1,  0.1);
  Sat.Term ( 5, -6,0,   0.1,    0.8,    -1.6,   0.1,    0.0,  0.0);
  Sat.Term ( 5, -9,0,  -0.5,   -0.1,     0.1,  -0.8,    0.0,  0.0);
  Sat.Term ( 5,-10,0,   2.5,   -2.2,     2.8,   3.1,    0.1, -0.2);

  dl = Sat.dl();  dr = Sat.dr();  db = Sat.db();


  // Perturbations by Uranus 
  Ura.Init ( T, M5,1,1, M7,-2,-1 );

  Ura.Term ( 1, -1,0,   0.4,    0.9,     0.0,   0.0,    0.0,  0.0);
  Ura.Term ( 1, -2,0,   0.4,    0.4,    -0.4,   0.3,    0.0,  0.0);

  dl += Ura.dl();  dr += Ura.dr();  db += Ura.db();


  // Perturbations by Saturn and Uranus
  phi = (2*M5-6*M6+3*M7); c=cos(phi); s=sin(phi);
  dl += -0.8*c+8.5*s;
  dr += -0.1*c;

  phi = (3*M5-6*M6+3*M7); c=cos(phi); s=sin(phi);
  dl += +0.4*c+0.5*s;
  dr += -0.7*c+0.5*s;
  db += -0.1*c;


  // Ecliptic coordinates ([rad],[AU])
  l = pi2 * Frac (0.0388910 + M5/pi2 + ( (5025.2+0.8*T)*T + dl ) / 1296.0e3);
  r = 5.208873 + 0.000041*T  +  dr * 1.0E-5;
  b = ( 227.3 - 0.3*T + db ) / Arcs;

  return Vec3D ( Polar(l,b,r) );  // Position vector

}


//------------------------------------------------------------------------------
//
// SaturnPos: Position of Saturn from analytical series
//
// Input:
//
//   T         Time in Julian centuries since J2000
//
// <return>:   Heliocentric position [AU], true-of-date ecliptic and equinox
//
//------------------------------------------------------------------------------
Vec3D SaturnPos (double T)
{
  //
  // Variables
  //
  double  M5,M6,M7,M8;             // Mean anomalies
  Pert    Jup,Ura,Nep;             // Perturbations
  double  phi,c,s;
  double  dl, dr, db;              // Corrections in longitude ["],
                                   // radius [AU] and latitude ["]
  double  l,b,r;                   // Ecliptic coordinates


  // Mean anomalies of planets in [rad]
  M5 = pi2 * Frac ( 0.0565314 + 8.4302963*T );
  M6 = pi2 * Frac ( 0.8829867 + 3.3947688*T );
  M7 = pi2 * Frac ( 0.3969537 + 1.1902586*T );
  M8 = pi2 * Frac ( 0.7208473 + 0.6068623*T );


  // Keplerian motion and perturbations by Jupiter
  Jup.Init ( T, M6,0,11, M5,-6,0 );

  Jup.Term ( 0,-1,0,   12.0,   -1.4,   -13.9,    6.4,    1.2,  -1.8);
  Jup.Term ( 0,-2,0,    0.0,   -0.2,    -0.9,    1.0,    0.0,  -0.1);
  Jup.Term ( 1, 1,0,    0.9,    0.4,    -1.8,    1.9,    0.2,   0.2);
  Jup.Term ( 1, 0,0, -348.3,22907.7,-52915.5, -752.2,-3266.5,8314.4);
  Jup.Term ( 1, 0,1, -225.2, -146.2,   337.7, -521.3,   79.6,  17.4);
  Jup.Term ( 1, 0,2,    1.3,   -1.4,     3.2,    2.9,    0.1,  -0.4);
  Jup.Term ( 1,-1,0,   -1.0,  -30.7,   108.6, -815.0,   -3.6,  -9.3);
  Jup.Term ( 1,-2,0,   -2.0,   -2.7,    -2.1,  -11.9,   -0.1,  -0.4);
  Jup.Term ( 2, 1,0,    0.1,    0.2,    -1.0,    0.3,    0.0,   0.0);
  Jup.Term ( 2, 0,0,   44.2,  724.0, -1464.3,  -34.7, -188.7, 459.1);
  Jup.Term ( 2, 0,1,  -17.0,  -11.3,    18.9,  -28.6,    1.0,  -3.7);
  Jup.Term ( 2,-1,0,   -3.5, -426.6,  -546.5,  -26.5,   -1.6,  -2.7);
  Jup.Term ( 2,-1,1,    3.5,   -2.2,    -2.6,   -4.3,    0.0,   0.0);
  Jup.Term ( 2,-2,0,   10.5,  -30.9,  -130.5,  -52.3,   -1.9,   0.2);
  Jup.Term ( 2,-3,0,   -0.2,   -0.4,    -1.2,   -0.1,   -0.1,   0.0);
  Jup.Term ( 3, 0,0,    6.5,   30.5,   -61.1,    0.4,  -11.6,  28.1);
  Jup.Term ( 3, 0,1,   -1.2,   -0.7,     1.1,   -1.8,   -0.2,  -0.6);
  Jup.Term ( 3,-1,0,   29.0,  -40.2,    98.2,   45.3,    3.2,  -9.4);
  Jup.Term ( 3,-1,1,    0.6,    0.6,    -1.0,    1.3,    0.0,   0.0);
  Jup.Term ( 3,-2,0,  -27.0,  -21.1,   -68.5,    8.1,  -19.8,   5.4);
  Jup.Term ( 3,-2,1,    0.9,   -0.5,    -0.4,   -2.0,   -0.1,  -0.8);
  Jup.Term ( 3,-3,0,   -5.4,   -4.1,   -19.1,   26.2,   -0.1,  -0.1);
  Jup.Term ( 4, 0,0,    0.6,    1.4,    -3.0,   -0.2,   -0.6,   1.6);
  Jup.Term ( 4,-1,0,    1.5,   -2.5,    12.4,    4.7,    1.0,  -1.1);
  Jup.Term ( 4,-2,0, -821.9,   -9.6,   -26.0, 1873.6,  -70.5,  -4.4);
  Jup.Term ( 4,-2,1,    4.1,  -21.9,   -50.3,   -9.9,    0.7,  -3.0);
  Jup.Term ( 4,-3,0,   -2.0,   -4.7,   -19.3,    8.2,   -0.1,  -0.3);
  Jup.Term ( 4,-4,0,   -1.5,    1.3,     6.5,    7.3,    0.0,   0.0);
  Jup.Term ( 5,-2,0,-2627.6,-1277.3,   117.4, -344.1,  -13.8,  -4.3);
  Jup.Term ( 5,-2,1,   63.0,  -98.6,    12.7,    6.7,    0.1,  -0.2);
  Jup.Term ( 5,-2,2,    1.7,    1.2,    -0.2,    0.3,    0.0,   0.0);
  Jup.Term ( 5,-3,0,    0.4,   -3.6,   -11.3,   -1.6,    0.0,  -0.3);
  Jup.Term ( 5,-4,0,   -1.4,    0.3,     1.5,    6.3,   -0.1,   0.0);
  Jup.Term ( 5,-5,0,    0.3,    0.6,     3.0,   -1.7,    0.0,   0.0);
  Jup.Term ( 6,-2,0, -146.7,  -73.7,   166.4, -334.3,  -43.6, -46.7);
  Jup.Term ( 6,-2,1,    5.2,   -6.8,    15.1,   11.4,    1.7,  -1.0);
  Jup.Term ( 6,-3,0,    1.5,   -2.9,    -2.2,   -1.3,    0.1,  -0.1);
  Jup.Term ( 6,-4,0,   -0.7,   -0.2,    -0.7,    2.8,    0.0,   0.0);
  Jup.Term ( 6,-5,0,    0.0,    0.5,     2.5,   -0.1,    0.0,   0.0);
  Jup.Term ( 6,-6,0,    0.3,   -0.1,    -0.3,   -1.2,    0.0,   0.0);
  Jup.Term ( 7,-2,0,   -9.6,   -3.9,     9.6,  -18.6,   -4.7,  -5.3);
  Jup.Term ( 7,-2,1,    0.4,   -0.5,     1.0,    0.9,    0.3,  -0.1);
  Jup.Term ( 7,-3,0,    3.0,    5.3,     7.5,   -3.5,    0.0,   0.0);
  Jup.Term ( 7,-4,0,    0.2,    0.4,     1.6,   -1.3,    0.0,   0.0);
  Jup.Term ( 7,-5,0,   -0.1,    0.2,     1.0,    0.5,    0.0,   0.0);
  Jup.Term ( 7,-6,0,    0.2,    0.0,     0.2,   -1.0,    0.0,   0.0);
  Jup.Term ( 8,-2,0,   -0.7,   -0.2,     0.6,   -1.2,   -0.4,  -0.4);
  Jup.Term ( 8,-3,0,    0.5,    1.0,    -2.0,    1.5,    0.1,   0.2);
  Jup.Term ( 8,-4,0,    0.4,    1.3,     3.6,   -0.9,    0.0,  -0.1);
  Jup.Term ( 9,-4,0,    4.0,   -8.7,   -19.9,   -9.9,    0.2,  -0.4);
  Jup.Term ( 9,-4,1,    0.5,    0.3,     0.8,   -1.8,    0.0,   0.0);
  Jup.Term (10,-4,0,   21.3,  -16.8,     3.3,    3.3,    0.2,  -0.2);
  Jup.Term (10,-4,1,    1.0,    1.7,    -0.4,    0.4,    0.0,   0.0);
  Jup.Term (11,-4,0,    1.6,   -1.3,     3.0,    3.7,    0.8,  -0.2);

  dl = Jup.dl();  dr = Jup.dr();  db = Jup.db();


  // Perturbations by Uranus 
  Ura.Init ( T, M6,0,3, M7,-5,-1 );

  Ura.Term ( 0,-1,0,    1.0,    0.7,     0.4,   -1.5,    0.1,   0.0);
  Ura.Term ( 0,-2,0,    0.0,   -0.4,    -1.1,    0.1,   -0.1,  -0.1);
  Ura.Term ( 0,-3,0,   -0.9,   -1.2,    -2.7,    2.1,   -0.5,  -0.3);
  Ura.Term ( 1,-1,0,    7.8,   -1.5,     2.3,   12.7,    0.0,   0.0);
  Ura.Term ( 1,-2,0,   -1.1,   -8.1,     5.2,   -0.3,   -0.3,  -0.3);
  Ura.Term ( 1,-3,0,  -16.4,  -21.0,    -2.1,    0.0,    0.4,   0.0);
  Ura.Term ( 2,-1,0,    0.6,   -0.1,     0.1,    1.2,    0.1,   0.0);
  Ura.Term ( 2,-2,0,   -4.9,  -11.7,    31.5,  -13.3,    0.0,  -0.2);
  Ura.Term ( 2,-3,0,   19.1,   10.0,   -22.1,   42.1,    0.1,  -1.1);
  Ura.Term ( 2,-4,0,    0.9,   -0.1,     0.1,    1.4,    0.0,   0.0);
  Ura.Term ( 3,-2,0,   -0.4,   -0.9,     1.7,   -0.8,    0.0,  -0.3);
  Ura.Term ( 3,-3,0,    2.3,    0.0,     1.0,    5.7,    0.3,   0.3);
  Ura.Term ( 3,-4,0,    0.3,   -0.7,     2.0,    0.7,    0.0,   0.0);
  Ura.Term ( 3,-5,0,   -0.1,   -0.4,     1.1,   -0.3,    0.0,   0.0);

  dl += Ura.dl();  dr += Ura.dr();  db += Ura.db();


  // Perturbations by Neptune
  Nep.Init ( T, M6,1,2, M8,-2,-1 );

  Nep.Term ( 1,-1,0,   -1.3,   -1.2,     2.3,   -2.5,    0.0,   0.0);
  Nep.Term ( 1,-2,0,    1.0,   -0.1,     0.1,    1.4,    0.0,   0.0);
  Nep.Term ( 2,-2,0,    1.1,   -0.1,     0.2,    3.3,    0.0,   0.0);

  dl += Nep.dl();  dr += Nep.dr();  db += Nep.db();


  // Perturbations by Jupiter and Uranus
  phi = (-2*M5+5*M6-3*M7); c=cos(phi); s=sin(phi);
 
  dl += -0.8*c-0.1*s;
  dr += -0.2*c+1.8*s;
  db += +0.3*c+0.5*s;

  phi = (-2*M5+6*M6-3*M7); c=cos(phi); s=sin(phi);
  dl += (+2.4-0.7*T)*c+(27.8-0.4*T)*s;
  dr += +2.1*c-0.2*s;
  
  phi = (-2*M5+7*M6-3*M7); c=cos(phi); s=sin(phi);
  dl += +0.1*c+1.6*s;
  dr += -3.6*c+0.3*s;
  db += -0.2*c+0.6*s;


  // Ecliptic coordinates ([rad],[AU])
  l = pi2 * Frac (0.2561136 + M6/pi2 + ( (5018.6+T*1.9)*T + dl ) / 1296.0e3);
  r = 9.557584 - 0.000186*T  +  dr*1.0e-5;
  b = ( 175.1 - 10.2*T + db ) / Arcs;

  return Vec3D ( Polar(l, b, r) );  // Position vector

}


//------------------------------------------------------------------------------
//
// UranusPos: Position of Uranus from analytical series
//
// Input:
//
//   T         Time in Julian centuries since J2000
//
// <return>:   Heliocentric position [AU], true-of-date ecliptic and equinox
//
//------------------------------------------------------------------------------
Vec3D UranusPos (double T)
{
  //
  // Variables
  //
  double  M5,M6,M7,M8;             // Mean anomalies
  Pert    Jup,Sat,Nep,JSU;         // Perturbations
  double  dl, dr, db;              // Corrections in longitude ["],
                                   // radius [AU] and latitude ["]
  double  l,b,r;                   // Ecliptic coordinates


  // Mean anomalies of planets in [rad]
  M5 = pi2 * Frac ( 0.0564472 + 8.4302889*T );
  M6 = pi2 * Frac ( 0.8829611 + 3.3947583*T );
  M7 = pi2 * Frac ( 0.3967117 + 1.1902849*T );
  M8 = pi2 * Frac ( 0.7216833 + 0.6068528*T );


  // Perturbations by Jupiter
  Jup.Init ( T, M7,-1,3, M5,-2,-1 );

  Jup.Term (-1,-1,0,  0.0,    0.0,    -0.1,   1.7,  -0.1,   0.0);
  Jup.Term ( 0,-1,0,  0.5,   -1.2,    18.9,   9.1,  -0.9,   0.1);
  Jup.Term ( 1,-1,0,-21.2,   48.7,  -455.5,-198.8,   0.0,   0.0);
  Jup.Term ( 1,-2,0, -0.5,    1.2,   -10.9,  -4.8,   0.0,   0.0);
  Jup.Term ( 2,-1,0, -1.3,    3.2,   -23.2, -11.1,   0.3,   0.1);
  Jup.Term ( 2,-2,0, -0.2,    0.2,     1.1,   1.5,   0.0,   0.0);
  Jup.Term ( 3,-1,0,  0.0,    0.2,    -1.8,   0.4,   0.0,   0.0);

  dl = Jup.dl();  dr = Jup.dr();  db = Jup.db();


  // Perturbations by Saturn
  Sat.Init ( T, M7,0,11, M6,-6,0 );

  Sat.Term ( 0,-1,0,  1.4,   -0.5,    -6.4,   9.0,  -0.4,  -0.8);
  Sat.Term ( 1,-1,0,-18.6,  -12.6,    36.7,-336.8,   1.0,   0.3);
  Sat.Term ( 1,-2,0, -0.7,   -0.3,     0.5,  -7.5,   0.1,   0.0);
  Sat.Term ( 2,-1,0, 20.0, -141.6,  -587.1,-107.0,   3.1,  -0.8);
  Sat.Term ( 2,-1,1,  1.0,    1.4,     5.8,  -4.0,   0.0,   0.0);
  Sat.Term ( 2,-2,0,  1.6,   -3.8,   -35.6, -16.0,   0.0,   0.0);
  Sat.Term ( 3,-1,0, 75.3, -100.9,   128.9,  77.5,  -0.8,   0.1);
  Sat.Term ( 3,-1,1,  0.2,    1.8,    -1.9,   0.3,   0.0,   0.0);
  Sat.Term ( 3,-2,0,  2.3,   -1.3,    -9.5, -17.9,   0.0,   0.1);
  Sat.Term ( 3,-3,0, -0.7,   -0.5,    -4.9,   6.8,   0.0,   0.0);
  Sat.Term ( 4,-1,0,  3.4,   -5.0,    21.6,  14.3,  -0.8,  -0.5);
  Sat.Term ( 4,-2,0,  1.9,    0.1,     1.2, -12.1,   0.0,   0.0);
  Sat.Term ( 4,-3,0, -0.1,   -0.4,    -3.9,   1.2,   0.0,   0.0);
  Sat.Term ( 4,-4,0, -0.2,    0.1,     1.6,   1.8,   0.0,   0.0);
  Sat.Term ( 5,-1,0,  0.2,   -0.3,     1.0,   0.6,  -0.1,   0.0);
  Sat.Term ( 5,-2,0, -2.2,   -2.2,    -7.7,   8.5,   0.0,   0.0);
  Sat.Term ( 5,-3,0,  0.1,   -0.2,    -1.4,  -0.4,   0.0,   0.0);
  Sat.Term ( 5,-4,0, -0.1,    0.0,     0.1,   1.2,   0.0,   0.0);
  Sat.Term ( 6,-2,0, -0.2,   -0.6,     1.4,  -0.7,   0.0,   0.0);

  dl += Sat.dl();  dr += Sat.dr();  db += Sat.db();


  // Keplerian motion and perturbations by Neptune
  Nep.Init ( T, M7,-1,7, M8,-8,0 );

  Nep.Term ( 1, 0,0,-78.1,19518.1,-90718.2,-334.7,2759.5,-311.9);
  Nep.Term ( 1, 0,1,-81.6,  107.7,  -497.4,-379.5,  -2.8, -43.7);
  Nep.Term ( 1, 0,2, -6.6,   -3.1,    14.4, -30.6,  -0.4,  -0.5);
  Nep.Term ( 1, 0,3,  0.0,   -0.5,     2.4,   0.0,   0.0,   0.0);
  Nep.Term ( 2, 0,0, -2.4,  586.1, -2145.2, -15.3, 130.6, -14.3);
  Nep.Term ( 2, 0,1, -4.5,    6.6,   -24.2, -17.8,   0.7,  -1.6);
  Nep.Term ( 2, 0,2, -0.4,    0.0,     0.1,  -1.4,   0.0,   0.0);
  Nep.Term ( 3, 0,0,  0.0,   24.5,   -76.2,  -0.6,   7.0,  -0.7);
  Nep.Term ( 3, 0,1, -0.2,    0.4,    -1.4,  -0.8,   0.1,  -0.1);
  Nep.Term ( 4, 0,0,  0.0,    1.1,    -3.0,   0.1,   0.4,   0.0);
  Nep.Term (-1,-1,0, -0.2,    0.2,     0.7,   0.7,  -0.1,   0.0);
  Nep.Term ( 0,-1,0, -2.8,    2.5,     8.7,  10.5,  -0.4,  -0.1);
  Nep.Term ( 1,-1,0,-28.4,   20.3,   -51.4, -72.0,   0.0,   0.0);
  Nep.Term ( 1,-2,0, -0.6,   -0.1,     4.2, -14.6,   0.2,   0.4);
  Nep.Term ( 1,-3,0,  0.2,    0.5,     3.4,  -1.6,  -0.1,   0.1);
  Nep.Term ( 2,-1,0, -1.8,    1.3,    -5.5,  -7.7,   0.0,   0.3);
  Nep.Term ( 2,-2,0, 29.4,   10.2,   -29.0,  83.2,   0.0,   0.0);
  Nep.Term ( 2,-3,0,  8.8,   17.8,   -41.9,  21.5,  -0.1,  -0.3);
  Nep.Term ( 2,-4,0,  0.0,    0.1,    -2.1,  -0.9,   0.1,   0.0);
  Nep.Term ( 3,-2,0,  1.5,    0.5,    -1.7,   5.1,   0.1,  -0.2);
  Nep.Term ( 3,-3,0,  4.4,   14.6,   -84.3,  25.2,   0.1,  -0.1);
  Nep.Term ( 3,-4,0,  2.4,   -4.5,    12.0,   6.2,   0.0,   0.0);
  Nep.Term ( 3,-5,0,  2.9,   -0.9,     2.1,   6.2,   0.0,   0.0);
  Nep.Term ( 4,-3,0,  0.3,    1.0,    -4.0,   1.1,   0.1,  -0.1);
  Nep.Term ( 4,-4,0,  2.1,   -2.7,    17.9,  14.0,   0.0,   0.0);
  Nep.Term ( 4,-5,0,  3.0,   -0.4,     2.3,  17.6,  -0.1,  -0.1);
  Nep.Term ( 4,-6,0, -0.6,   -0.5,     1.1,  -1.6,   0.0,   0.0);
  Nep.Term ( 5,-4,0,  0.2,   -0.2,     1.0,   0.8,   0.0,   0.0);
  Nep.Term ( 5,-5,0, -0.9,   -0.1,     0.6,  -7.1,   0.0,   0.0);
  Nep.Term ( 5,-6,0, -0.5,   -0.6,     3.8,  -3.6,   0.0,   0.0);
  Nep.Term ( 5,-7,0,  0.0,   -0.5,     3.0,   0.1,   0.0,   0.0);
  Nep.Term ( 6,-6,0,  0.2,    0.3,    -2.7,   1.6,   0.0,   0.0);
  Nep.Term ( 6,-7,0, -0.1,    0.2,    -2.0,  -0.4,   0.0,   0.0);
  Nep.Term ( 7,-7,0,  0.1,   -0.2,     1.3,   0.5,   0.0,   0.0);
  Nep.Term ( 7,-8,0,  0.1,    0.0,     0.4,   0.9,   0.0,   0.0);

  dl += Nep.dl();  dr += Nep.dr();  db += Nep.db();


  // Perturbations by Jupiter and Uranus
  JSU.Init ( T, M7,-2,4, M6,-6,-4, 2*M5 );

  JSU.Term (-2,-4,0, -0.7,    0.4,    -1.5,  -2.5,   0.0,   0.0);
  JSU.Term (-1,-4,0, -0.1,   -0.1,    -2.2,   1.0,   0.0,   0.0);
  JSU.Term ( 1,-5,0,  0.1,   -0.4,     1.4,   0.2,   0.0,   0.0);
  JSU.Term ( 1,-6,0,  0.4,    0.5,    -0.8,  -0.8,   0.0,   0.0);
  JSU.Term ( 2,-6,0,  5.7,    6.3,    28.5, -25.5,   0.0,   0.0);
  JSU.Term ( 2,-6,1,  0.1,   -0.2,    -1.1,  -0.6,   0.0,   0.0);
  JSU.Term ( 3,-6,0, -1.4,   29.2,   -11.4,   1.1,   0.0,   0.0);
  JSU.Term ( 3,-6,1,  0.8,   -0.4,     0.2,   0.3,   0.0,   0.0);
  JSU.Term ( 4,-6,0,  0.0,    1.3,    -6.0,  -0.1,   0.0,   0.0);

  dl += JSU.dl();  dr += JSU.dr();  db += JSU.db();


  // Ecliptic coordinates ([rad],[AU])

  l = pi2 * Frac (0.4734843 + M7/pi2 + ( (5082.3+34.2*T)*T + dl ) / 1296.0e3);
  r = 19.211991 + (-0.000333-0.000005*T)*T  +  dr*1.0e-5;
  b = ( -130.61 + (-0.54+0.04*T)*T + db ) / Arcs;

  return Vec3D ( Polar(l, b, r) );  // Position vector

}


//------------------------------------------------------------------------------
//
// NeptunePos: Position of Neptune from analytical series
//
// Input:
//
//   T         Time in Julian centuries since J2000
//
// <return>:   Heliocentric position [AU], true-of-date ecliptic and equinox
//
//------------------------------------------------------------------------------
Vec3D NeptunePos (double T)
{
  //
  // Variables
  //
  double  M5,M6,M7,M8;             // Mean anomalies
  Pert    Jup,Sat,Ura;             // Perturbations
  double  dl, dr, db;              // Corrections in longitude ["],
                                   // radius [AU] and latitude ["]
  double  l,b,r;                   // Ecliptic coordinates


  // Mean anomalies of planets in [rad]
  M5 = pi2 * Frac ( 0.0563867 + 8.4298907*T );
  M6 = pi2 * Frac ( 0.8825086 + 3.3957748*T );
  M7 = pi2 * Frac ( 0.3965358 + 1.1902851*T );
  M8 = pi2 * Frac ( 0.7214906 + 0.6068526*T );


  // Perturbations by Jupiter
  Jup.Init ( T, M8,0,2, M5,-2,0 );

  Jup.Term ( 0,-1,0,  0.1,   0.1,    -3.0,   1.8,   -0.3, -0.3);
  Jup.Term ( 1, 0,0,  0.0,   0.0,   -15.9,   9.0,    0.0,  0.0);
  Jup.Term ( 1,-1,0,-17.6, -29.3,   416.1,-250.0,    0.0,  0.0);
  Jup.Term ( 1,-2,0, -0.4,  -0.7,    10.4,  -6.2,    0.0,  0.0);
  Jup.Term ( 2,-1,0, -0.2,  -0.4,     2.4,  -1.4,    0.4, -0.3);

  dl = Jup.dl();  dr = Jup.dr();  db = Jup.db();


  // Perturbations by Saturn
  Sat.Init ( T, M8,0,2, M6,-2,0 );

  Sat.Term ( 0,-1,0, -0.1,   0.0,     0.2,  -1.8,   -0.1, -0.5);
  Sat.Term ( 1, 0,0,  0.0,   0.0,    -8.3, -10.4,    0.0,  0.0);
  Sat.Term ( 1,-1,0, 13.6, -12.7,   187.5, 201.1,    0.0,  0.0);
  Sat.Term ( 1,-2,0,  0.4,  -0.4,     4.5,   4.5,    0.0,  0.0);
  Sat.Term ( 2,-1,0,  0.4,  -0.1,     1.7,  -3.2,    0.2,  0.2);
  Sat.Term ( 2,-2,0, -0.1,   0.0,    -0.2,   2.7,    0.0,  0.0);

  dl += Sat.dl();  dr += Sat.dr();  db += Sat.db();


  // Keplerian motion and perturbations by Uranus
  Ura.Init ( T, M8,1,6, M7,-6,0 );

  Ura.Term ( 1, 0,0, 32.3,3549.5,-25880.2, 235.8,-6360.5,374.0);
  Ura.Term ( 1, 0,1, 31.2,  34.4,  -251.4, 227.4,   34.9, 29.3);
  Ura.Term ( 1, 0,2, -1.4,   3.9,   -28.6, -10.1,    0.0, -0.9);
  Ura.Term ( 2, 0,0,  6.1,  68.0,  -111.4,   2.0,  -54.7,  3.7);
  Ura.Term ( 2, 0,1,  0.8,  -0.2,    -2.1,   2.0,   -0.2,  0.8);
  Ura.Term ( 3, 0,0,  0.1,   1.0,    -0.7,   0.0,   -0.8,  0.1);
  Ura.Term ( 0,-1,0, -0.1,  -0.3,    -3.6,   0.0,    0.0,  0.0);
  Ura.Term ( 1, 0,0,  0.0,   0.0,     5.5,  -6.9,    0.1,  0.0);
  Ura.Term ( 1,-1,0, -2.2,  -1.6,  -116.3, 163.6,    0.0, -0.1);
  Ura.Term ( 1,-2,0,  0.2,   0.1,    -1.2,   0.4,    0.0, -0.1);
  Ura.Term ( 2,-1,0,  4.2,  -1.1,    -4.4, -34.6,   -0.2,  0.1);
  Ura.Term ( 2,-2,0,  8.6,  -2.9,   -33.4, -97.0,    0.2,  0.1);
  Ura.Term ( 3,-1,0,  0.1,  -0.2,     2.1,  -1.2,    0.0,  0.1);
  Ura.Term ( 3,-2,0, -4.6,   9.3,    38.2,  19.8,    0.1,  0.1);
  Ura.Term ( 3,-3,0, -0.5,   1.7,    23.5,   7.0,    0.0,  0.0);
  Ura.Term ( 4,-2,0,  0.2,   0.8,     3.3,  -1.5,   -0.2, -0.1);
  Ura.Term ( 4,-3,0,  0.9,   1.7,    17.9,  -9.1,   -0.1,  0.0);
  Ura.Term ( 4,-4,0, -0.4,  -0.4,    -6.2,   4.8,    0.0,  0.0);
  Ura.Term ( 5,-3,0, -1.6,  -0.5,    -2.2,   7.0,    0.0,  0.0);
  Ura.Term ( 5,-4,0, -0.4,  -0.1,    -0.7,   5.5,    0.0,  0.0);
  Ura.Term ( 5,-5,0,  0.2,   0.0,     0.0,  -3.5,    0.0,  0.0);
  Ura.Term ( 6,-4,0, -0.3,   0.2,     2.1,   2.7,    0.0,  0.0);
  Ura.Term ( 6,-5,0,  0.1,  -0.1,    -1.4,  -1.4,    0.0,  0.0);
  Ura.Term ( 6,-6,0, -0.1,   0.1,     1.4,   0.7,    0.0,  0.0);

  dl += Ura.dl();  dr += Ura.dr();  db += Ura.db();


  // Ecliptic coordinates ([rad],[AU])
  l = pi2 * Frac (0.1254046 + M8/pi2 + ( (4982.8-21.3*T)*T + dl ) / 1296.0e3);
  r = 30.072984 + (0.001234+0.000003*T) * T  +  dr*1.0e-5;
  b = (  54.77 + ( 0.26 + 0.06*T) * T  +  db ) / Arcs;

  return Vec3D ( Polar(l, b, r) );  // Position vector

}


//------------------------------------------------------------------------------
//
// PlutoPos: Position of Pluto from analytical series
//
// Input:
//
//   T         Time in Julian centuries since J2000
//
// <return>:   Heliocentric position [AU], true-of-date ecliptic and equinox
//
// Notes:
//
//   The series are only valid from 1890 to 2100
//
//------------------------------------------------------------------------------
Vec3D PlutoPos (double T)
{
  //
  // Variables
  //
  double  M5,M6,M9;                // Mean anomalies
  Pert    Jup,Sat;                 // Perturbations
  double  phi,c,s;
  double  dl, dr, db;              // Corrections in longitude ["],
                                   // radius [AU] and latitude ["]
  double  l,b,r;                   // Ecliptic coordinates
  Vec3D   R;                       // Position vector


  // Mean anomalies of planets in [rad]
  M5 = pi2 * Frac ( 0.0565314 + 8.4302963*T );
  M6 = pi2 * Frac ( 0.8829867 + 3.3947688*T );
  M9 = pi2 * Frac ( 0.0385795 + 0.4026667*T );


  // Perturbations by Jupiter
  Jup.Init ( T, M9,0,6, M5,-2,1 );

  Jup.Term ( 1, 0,0,   0.06,100924.08,-960396.0,15965.1,51987.68,-24288.76);
  Jup.Term ( 2, 0,0,3274.74, 17835.12,-118252.2, 3632.4,12687.49, -6049.72);
  Jup.Term ( 3, 0,0,1543.52,  4631.99, -21446.6, 1167.0, 3504.00, -1853.10);
  Jup.Term ( 4, 0,0, 688.99,  1227.08,  -4823.4,  213.5, 1048.19,  -648.26);
  Jup.Term ( 5, 0,0, 242.27,   415.93,  -1075.4,  140.6,  302.33,  -209.76);
  Jup.Term ( 6, 0,0, 138.41,   110.91,   -308.8,  -55.3,  109.52,   -93.82);
  Jup.Term ( 3,-1,0,  -0.99,     5.06,    -25.6,   19.8,    1.26,    -1.96);
  Jup.Term ( 2,-1,0,   7.15,     5.61,    -96.7,   57.2,    1.64,    -2.16);
  Jup.Term ( 1,-1,0,  10.79,    23.13,   -390.4,  236.4,   -0.33,     0.86);
  Jup.Term ( 0, 1,0,  -0.23,     4.43,    102.8,   63.2,    3.15,     0.34);
  Jup.Term ( 1, 1,0,  -1.10,    -0.92,     11.8,   -2.3,    0.43,     0.14);
  Jup.Term ( 2, 1,0,   0.62,     0.84,      2.3,    0.7,    0.05,    -0.04);
  Jup.Term ( 3, 1,0,  -0.38,    -0.45,      1.2,   -0.8,    0.04,     0.05);
  Jup.Term ( 4, 1,0,   0.17,     0.25,      0.0,    0.2,   -0.01,    -0.01);
  Jup.Term ( 3,-2,0,   0.06,     0.07,     -0.6,    0.3,    0.03,    -0.03);
  Jup.Term ( 2,-2,0,   0.13,     0.20,     -2.2,    1.5,    0.03,    -0.07);
  Jup.Term ( 1,-2,0,   0.32,     0.49,     -9.4,    5.7,   -0.01,     0.03);
  Jup.Term ( 0,-2,0,  -0.04,    -0.07,      2.6,   -1.5,    0.07,    -0.02);

  dl = Jup.dl();  dr = Jup.dr();  db = Jup.db();


  // Perturbations by Saturn
  Sat.Init ( T, M9,0,3, M6,-2,1 );

  Sat.Term ( 1,-1,0, -29.47,    75.97,   -106.4, -204.9,  -40.71,   -17.55);
  Sat.Term ( 0, 1,0, -13.88,    18.20,     42.6,  -46.1,    1.13,     0.43);
  Sat.Term ( 1, 1,0,   5.81,   -23.48,     15.0,   -6.8,   -7.48,     3.07);
  Sat.Term ( 2, 1,0, -10.27,    14.16,     -7.9,    0.4,    2.43,    -0.09);
  Sat.Term ( 3, 1,0,   6.86,   -10.66,      7.3,   -0.3,   -2.25,     0.69);
  Sat.Term ( 2,-2,0,   4.32,     2.00,      0.0,   -2.2,   -0.24,     0.12);
  Sat.Term ( 1,-2,0,  -5.04,    -0.83,     -9.2,   -3.1,    0.79,    -0.24);
  Sat.Term ( 0,-2,0,   4.25,     2.48,     -5.9,   -3.3,    0.58,     0.02);

  dl += Sat.dl();  dr += Sat.dr();  db += Sat.db();


  // Perturbations by Jupiter and Saturn
  phi=(M5-M6); c=cos(phi); s=sin(phi);
  dl += -9.11*c+0.12*s;
  dr += -3.4*c-3.3*s;
  db += +0.81*c+0.78*s;

  phi=(M5-M6+M9); c=cos(phi); s=sin(phi);
  dl += +5.92*c+0.25*s;
  dr += +2.3*c-3.8*s;
  db += -0.67*c-0.51*s;


  // Ecliptic coordinates ([rad],[AU])
  l = pi2 * Frac ( 0.6232469 + M9/pi2 + dl/1296.0E3 );
  r = 40.7247248  +  dr * 1.0E-5;
  b = -3.909434*Rad  +  db / Arcs;

  // Position vector; ecliptic and equinox of B1950.0
  R = Vec3D ( Polar(l,b,r) );

  // Position vector; ecliptic and equinox of date
  return PrecMatrix_Ecl(T_B1950,T) * R;
}


//------------------------------------------------------------------------------
//
// PertPosition: Planetary position from analytical perturbation series
//
// Input:
//
//   Planet    Identifies the planet
//   T         Time in Julian centuries since J2000
//
// <return>:   Heliocentric position [AU], true-of-date ecliptic and equinox
//
//------------------------------------------------------------------------------
Vec3D PertPosition (PlanetType Planet, double T)
{
  Vec3D r; // Null vector;

  switch ( Planet )
  {
    case Sun:      break;
    case Mercury:  r = MercuryPos(T);  break;
    case Venus:    r = VenusPos(T);    break;
    case Earth:    r = (-SunPos(T));   break;
    case Mars:     r = MarsPos(T);     break;
    case Jupiter:  r = JupiterPos(T);  break;
    case Saturn:   r = SaturnPos(T);   break;
    case Uranus:   r = UranusPos(T);   break;
    case Neptune:  r = NeptunePos(T);  break;
    case Pluto:    r = PlutoPos(T);
  }
  
  return r;
};
