//------------------------------------------------------------------------------
//
// APC_Moon.cpp
// 
//    Computation of the Moon's position
//
//------------------------------------------------------------------------------

#include <cmath>

#include <math.h>
#include <vcl.h>
#include <windows.h>    // Header file for windows
#include <stdio.h>      // Header file for standard Input/Output

#include "APC_Const.h"
#include "APC_Math.h"
#include "APC_Moon.h"
#include "APC_PrecNut.h"
#include "APC_Spheric.h"
#include "APC_VecMat3D.h"


namespace // Unnamed namespace
{
  //
  // Constants
  //
  const int o   = 6;          // Index offset
  const int dim = 2*o+1;      // Work array dimension


  //
  // Sine
  //
  double Sine (double x) { return sin(pi2*Frac(x)); }


  //
  // Definition of ILE_Pert class for summing up Brown's lunar perturbation 
  // series
  //
  class ILE_Pert
  {

    public:

      // Initialization (mean arguments, long-periodic corrections, etc.)
      void Init (double T);

      // Perturbation term 
      void Term (int p, int q, int r, int s, double& x, double& y);

      // Summation of solar perturbations
      void AddSol ( double coeffl, double coeffS, double coeffg, 
                    double coeffP, int p, int q, int r, int s );

      // Summation of perturbation in latitude
      void AddN (double coeffN, int p, int q, int r, int s);

      // Planetary perturbations
      void Planetary (double T);

      // Coordinates
      double lambda();  // ecliptic longitude in [rad]
      double beta();    // ecliptic latitude in [rad]
      double dist();    // geocentric distance in [km]

    private:

      double Dgam;                       // Longperiodic perturbation
      double Dlam, DS, gam1C, sinPi, N;  // Periodic perturbations
      double L0, l,ls,F,D;               // Mean arguments of lunar orbit
      double Cos[dim][4], Sin[dim][4];   // Cosine and sine of mean arguments
  };


  //
  // Initialization
  //
  void ILE_Pert::Init (double T)
  {
    //
    // Variables
    //
    double dL0, dl, dls, dF, dD;             // Longperiodic perturbations
    double T2, arg, fac;                     // Auxiliary variables
    double S1, S2, S3, S4, S5, S6, S7;
    int    max;

  
    T2=T*T; // Time

  
    // Reset perturbations
    Dlam=0.0; DS=0.0; gam1C=0.0; sinPi=3422.7000; N=0.0;

  
    // Longperiodic perturbations
    S1 = Sine (0.19833+0.05611*T);  S2 = Sine (0.27869+0.04508*T);
    S3 = Sine (0.16827-0.36903*T);  S4 = Sine (0.34734-5.37261*T);
    S5 = Sine (0.10498-5.37899*T);  S6 = Sine (0.42681-0.41855*T);
    S7 = Sine (0.14943-5.37511*T); 

    dL0 = 0.84*S1+0.31*S2+14.27*S3+ 7.26*S4+ 0.28*S5+0.24*S6;
    dl  = 2.94*S1+0.31*S2+14.27*S3+ 9.34*S4+ 1.12*S5+0.83*S6;
    dls =-6.40*S1                                   -1.89*S6;
    dF  = 0.21*S1+0.31*S2+14.27*S3-88.70*S4-15.30*S5+0.24*S6-1.86*S7;
    dD  = dL0-dls;

    Dgam   = -3332e-9 * Sine (0.59734-5.37261*T)
              -539e-9 * Sine (0.35498-5.37899*T)
               -64e-9 * Sine (0.39943-5.37511*T);


    // Mean arguments of the lunar orbit (incl. longperiodic corrections)
    // L0 mean longitude of the Moon
    // l  mean anomaly of the Moon     l' mean anomaly of the Sun      
    // F  mean distance from the node  D  mean elongation from the Sun 

    L0 = pi2*Frac(0.60643382+1336.85522467*T-0.00000313*T2) + dL0/Arcs;
    l  = pi2*Frac(0.37489701+1325.55240982*T+0.00002565*T2) + dl /Arcs;
    ls = pi2*Frac(0.99312619+  99.99735956*T-0.00000044*T2) + dls/Arcs;
    F  = pi2*Frac(0.25909118+1342.22782980*T-0.00000892*T2) + dF /Arcs;
    D  = pi2*Frac(0.82736186+1236.85308708*T-0.00000397*T2) + dD /Arcs;

  
    // Cosine and sine of multiples of mean arguments 
    // incl. secular correction
    for (int i=0; i<=3; i++) {
      switch(i) {      
        case 0: arg=l;  max=4; fac=1.000002208;               break;
        case 1: arg=ls; max=3; fac=0.997504612-0.002495388*T; break;
        case 2: arg=F;  max=4; fac=1.000002708+139.978*Dgam;  break;
        case 3: arg=D;  max=6; fac=1.0;                       break;
      };
      
      Cos[o][i]=1.0;  Cos[o+1][i]=cos(arg)*fac;  Cos[o-1][i]=+Cos[o+1][i];
      Sin[o][i]=0.0;  Sin[o+1][i]=sin(arg)*fac;  Sin[o-1][i]=-Sin[o+1][i];
      
      for (int j=2;j<=max;j++) {
        AddThe ( Cos[o+j-1][i],Sin[o+j-1][i], Cos[o+1][i],Sin[o+1][i],
                 Cos[o+j][i],Sin[o+j][i] ); 
        Cos[o-j][i]=+Cos[o+j][i];
        Sin[o-j][i]=-Sin[o+j][i];
      };
    };
  }


  //
  // Term: calculates x=cos(p*l+q*ls+r*F+s*D) and y=sin(p*l+q*ls+r*F+s*D)
  //
  //
  void ILE_Pert::Term (int p, int q, int r, int s, double& x, double& y)
  {
    int i[4];

    i[0]=p; i[1]=q; i[2]=r; i[3]=s;  x=1.0; y=0.0;

    for (int k=0; k<=3; k++) 
      if (i[k]!=0) AddThe(x,y,Cos[o+i[k]][k],Sin[o+i[k]][k],x,y);
  }


  //
  // AddSol: Summation of solar perturbations
  //
  void ILE_Pert::AddSol ( 
    double coeffl, double coeffS, double coeffg, double coeffP,
    int p, int q, int r, int s )
  {
    //
    // Variables
    //
    double x,y;


    Term (p,q,r,s,x,y);
    Dlam  += coeffl*y; DS    += coeffS*y;
    gam1C += coeffg*x; sinPi += coeffP*x;
  }


  //
  // AddN: Summation of perturbation in latitude
  //
  void ILE_Pert::AddN (double coeffN, int p, int q, int r, int s)
  {
    //
    // Variables
    //
    double x,y;


    Term(p,q,r,s,x,y); 
    N += coeffN*y;
  }


  //
  // Planetary: Perturbations of ecliptic latitude by Venus and Jupiter
  //
  //
  void ILE_Pert::Planetary (double T)
  {
    Dlam +=
          +0.82*Sine(0.7736  -62.5512*T)+0.31*Sine(0.0466 -125.1025*T)
          +0.35*Sine(0.5785  -25.1042*T)+0.66*Sine(0.4591+1335.8075*T)
          +0.64*Sine(0.3130  -91.5680*T)+1.14*Sine(0.1480+1331.2898*T)
          +0.21*Sine(0.5918+1056.5859*T)+0.44*Sine(0.5784+1322.8595*T)
          +0.24*Sine(0.2275   -5.7374*T)+0.28*Sine(0.2965   +2.6929*T)
          +0.33*Sine(0.3132   +6.3368*T);
  }  


  //
  // lambda, beta, dist
  //
  double ILE_Pert::lambda()
  { 
    return Modulo ( L0+Dlam/Arcs, pi2 ); 
  }


  double ILE_Pert::beta()
  {
    //
    // Variables
    //
    double S   = F + DS/Arcs;
    double fac = 1.000002708+139.978*Dgam;
    
    
    return (fac*(18518.511+1.189+gam1C)*sin(S)-6.24*sin(3*S)+N) / Arcs;
  }


  double ILE_Pert::dist()
  {
    return R_Earth * Arcs / (sinPi * 0.999953253);
  }

} // End of unnamed namespace


//------------------------------------------------------------------------------
//
// MoonPos: Computes the Moon's ecliptic position using Brown's theory
//          (Improved Lunar Ephemeris)
//
// Input:
//
//   T         Time in Julian centuries since J2000
//
// <return>:   Geocentric position of the Moon (in [km]) referred to the
//             ecliptic and equinox of date
//
// Notes: Light-time is already taken into account
//
//------------------------------------------------------------------------------
Vec3D MoonPos (double T)
{
  //
  // Variables
  //
  ILE_Pert Pert;

  
  Pert.Init(T);  // Initialization


  // Solar perturbations
  Pert.AddSol (   13.902,   14.06,-0.001,   0.2607,0, 0, 0, 4);
  Pert.AddSol (    0.403,   -4.01,+0.394,   0.0023,0, 0, 0, 3);
  Pert.AddSol ( 2369.912, 2373.36,+0.601,  28.2333,0, 0, 0, 2);
  Pert.AddSol ( -125.154, -112.79,-0.725,  -0.9781,0, 0, 0, 1);
  Pert.AddSol (    1.979,    6.98,-0.445,   0.0433,1, 0, 0, 4);
  Pert.AddSol (  191.953,  192.72,+0.029,   3.0861,1, 0, 0, 2);
  Pert.AddSol (   -8.466,  -13.51,+0.455,  -0.1093,1, 0, 0, 1);
  Pert.AddSol (22639.500,22609.07,+0.079, 186.5398,1, 0, 0, 0);
  Pert.AddSol (   18.609,    3.59,-0.094,   0.0118,1, 0, 0,-1);
  Pert.AddSol (-4586.465,-4578.13,-0.077,  34.3117,1, 0, 0,-2);
  Pert.AddSol (   +3.215,    5.44,+0.192,  -0.0386,1, 0, 0,-3);
  Pert.AddSol (  -38.428,  -38.64,+0.001,   0.6008,1, 0, 0,-4);
  Pert.AddSol (   -0.393,   -1.43,-0.092,   0.0086,1, 0, 0,-6);
  Pert.AddSol (   -0.289,   -1.59,+0.123,  -0.0053,0, 1, 0, 4);
  Pert.AddSol (  -24.420,  -25.10,+0.040,  -0.3000,0, 1, 0, 2);
  Pert.AddSol (   18.023,   17.93,+0.007,   0.1494,0, 1, 0, 1);
  Pert.AddSol ( -668.146, -126.98,-1.302,  -0.3997,0, 1, 0, 0);
  Pert.AddSol (    0.560,    0.32,-0.001,  -0.0037,0, 1, 0,-1);
  Pert.AddSol ( -165.145, -165.06,+0.054,   1.9178,0, 1, 0,-2);
  Pert.AddSol (   -1.877,   -6.46,-0.416,   0.0339,0, 1, 0,-4);
  Pert.AddSol (    0.213,    1.02,-0.074,   0.0054,2, 0, 0, 4);
  Pert.AddSol (   14.387,   14.78,-0.017,   0.2833,2, 0, 0, 2);
  Pert.AddSol (   -0.586,   -1.20,+0.054,  -0.0100,2, 0, 0, 1);
  Pert.AddSol (  769.016,  767.96,+0.107,  10.1657,2, 0, 0, 0);
  Pert.AddSol (   +1.750,    2.01,-0.018,   0.0155,2, 0, 0,-1);
  Pert.AddSol ( -211.656, -152.53,+5.679,  -0.3039,2, 0, 0,-2);
  Pert.AddSol (   +1.225,    0.91,-0.030,  -0.0088,2, 0, 0,-3);
  Pert.AddSol (  -30.773,  -34.07,-0.308,   0.3722,2, 0, 0,-4);
  Pert.AddSol (   -0.570,   -1.40,-0.074,   0.0109,2, 0, 0,-6);
  Pert.AddSol (   -2.921,  -11.75,+0.787,  -0.0484,1, 1, 0, 2);
  Pert.AddSol (   +1.267,    1.52,-0.022,   0.0164,1, 1, 0, 1);
  Pert.AddSol ( -109.673, -115.18,+0.461,  -0.9490,1, 1, 0, 0);
  Pert.AddSol ( -205.962, -182.36,+2.056,  +1.4437,1, 1, 0,-2);
  Pert.AddSol (    0.233,    0.36, 0.012,  -0.0025,1, 1, 0,-3);
  Pert.AddSol (   -4.391,   -9.66,-0.471,   0.0673,1, 1, 0,-4);
  Pert.AddSol (    0.283,    1.53,-0.111,  +0.0060,1,-1, 0,+4);
  Pert.AddSol (   14.577,   31.70,-1.540,  +0.2302,1,-1, 0, 2);
  Pert.AddSol (  147.687,  138.76,+0.679,  +1.1528,1,-1, 0, 0);
  Pert.AddSol (   -1.089,    0.55,+0.021,   0.0   ,1,-1, 0,-1);
  Pert.AddSol (   28.475,   23.59,-0.443,  -0.2257,1,-1, 0,-2);
  Pert.AddSol (   -0.276,   -0.38,-0.006,  -0.0036,1,-1, 0,-3);
  Pert.AddSol (    0.636,    2.27,+0.146,  -0.0102,1,-1, 0,-4);
  Pert.AddSol (   -0.189,   -1.68,+0.131,  -0.0028,0, 2, 0, 2);
  Pert.AddSol (   -7.486,   -0.66,-0.037,  -0.0086,0, 2, 0, 0);
  Pert.AddSol (   -8.096,  -16.35,-0.740,   0.0918,0, 2, 0,-2);
  Pert.AddSol (   -5.741,   -0.04, 0.0  ,  -0.0009,0, 0, 2, 2);
  Pert.AddSol (    0.255,    0.0 , 0.0  ,   0.0   ,0, 0, 2, 1);
  Pert.AddSol ( -411.608,   -0.20, 0.0  ,  -0.0124,0, 0, 2, 0);
  Pert.AddSol (    0.584,    0.84, 0.0  ,  +0.0071,0, 0, 2,-1);
  Pert.AddSol (  -55.173,  -52.14, 0.0  ,  -0.1052,0, 0, 2,-2);
  Pert.AddSol (    0.254,    0.25, 0.0  ,  -0.0017,0, 0, 2,-3);
  Pert.AddSol (   +0.025,   -1.67, 0.0  ,  +0.0031,0, 0, 2,-4);
  Pert.AddSol (    1.060,    2.96,-0.166,   0.0243,3, 0, 0,+2);
  Pert.AddSol (   36.124,   50.64,-1.300,   0.6215,3, 0, 0, 0);
  Pert.AddSol (  -13.193,  -16.40,+0.258,  -0.1187,3, 0, 0,-2);
  Pert.AddSol (   -1.187,   -0.74,+0.042,   0.0074,3, 0, 0,-4);
  Pert.AddSol (   -0.293,   -0.31,-0.002,   0.0046,3, 0, 0,-6);
  Pert.AddSol (   -0.290,   -1.45,+0.116,  -0.0051,2, 1, 0, 2);
  Pert.AddSol (   -7.649,  -10.56,+0.259,  -0.1038,2, 1, 0, 0);
  Pert.AddSol (   -8.627,   -7.59,+0.078,  -0.0192,2, 1, 0,-2);
  Pert.AddSol (   -2.740,   -2.54,+0.022,   0.0324,2, 1, 0,-4);
  Pert.AddSol (    1.181,    3.32,-0.212,   0.0213,2,-1, 0,+2);
  Pert.AddSol (    9.703,   11.67,-0.151,   0.1268,2,-1, 0, 0);
  Pert.AddSol (   -0.352,   -0.37,+0.001,  -0.0028,2,-1, 0,-1);
  Pert.AddSol (   -2.494,   -1.17,-0.003,  -0.0017,2,-1, 0,-2);
  Pert.AddSol (    0.360,    0.20,-0.012,  -0.0043,2,-1, 0,-4);
  Pert.AddSol (   -1.167,   -1.25,+0.008,  -0.0106,1, 2, 0, 0);
  Pert.AddSol (   -7.412,   -6.12,+0.117,   0.0484,1, 2, 0,-2);
  Pert.AddSol (   -0.311,   -0.65,-0.032,   0.0044,1, 2, 0,-4);
  Pert.AddSol (   +0.757,    1.82,-0.105,   0.0112,1,-2, 0, 2);
  Pert.AddSol (   +2.580,    2.32,+0.027,   0.0196,1,-2, 0, 0);
  Pert.AddSol (   +2.533,    2.40,-0.014,  -0.0212,1,-2, 0,-2);
  Pert.AddSol (   -0.344,   -0.57,-0.025,  +0.0036,0, 3, 0,-2);
  Pert.AddSol (   -0.992,   -0.02, 0.0  ,   0.0   ,1, 0, 2, 2);
  Pert.AddSol (  -45.099,   -0.02, 0.0  ,  -0.0010,1, 0, 2, 0);
  Pert.AddSol (   -0.179,   -9.52, 0.0  ,  -0.0833,1, 0, 2,-2);
  Pert.AddSol (   -0.301,   -0.33, 0.0  ,   0.0014,1, 0, 2,-4);
  Pert.AddSol (   -6.382,   -3.37, 0.0  ,  -0.0481,1, 0,-2, 2);
  Pert.AddSol (   39.528,   85.13, 0.0  ,  -0.7136,1, 0,-2, 0);
  Pert.AddSol (    9.366,    0.71, 0.0  ,  -0.0112,1, 0,-2,-2);
  Pert.AddSol (    0.202,    0.02, 0.0  ,   0.0   ,1, 0,-2,-4);
  Pert.AddSol (    0.415,    0.10, 0.0  ,   0.0013,0, 1, 2, 0);
  Pert.AddSol (   -2.152,   -2.26, 0.0  ,  -0.0066,0, 1, 2,-2);
  Pert.AddSol (   -1.440,   -1.30, 0.0  ,  +0.0014,0, 1,-2, 2);
  Pert.AddSol (    0.384,   -0.04, 0.0  ,   0.0   ,0, 1,-2,-2);
  Pert.AddSol (   +1.938,   +3.60,-0.145,  +0.0401,4, 0, 0, 0);
  Pert.AddSol (   -0.952,   -1.58,+0.052,  -0.0130,4, 0, 0,-2);
  Pert.AddSol (   -0.551,   -0.94,+0.032,  -0.0097,3, 1, 0, 0);
  Pert.AddSol (   -0.482,   -0.57,+0.005,  -0.0045,3, 1, 0,-2);
  Pert.AddSol (    0.681,    0.96,-0.026,   0.0115,3,-1, 0, 0);
  Pert.AddSol (   -0.297,   -0.27, 0.002,  -0.0009,2, 2, 0,-2);
  Pert.AddSol (    0.254,   +0.21,-0.003,   0.0   ,2,-2, 0,-2);
  Pert.AddSol (   -0.250,   -0.22, 0.004,   0.0014,1, 3, 0,-2);
  Pert.AddSol (   -3.996,    0.0 , 0.0  ,  +0.0004,2, 0, 2, 0);
  Pert.AddSol (    0.557,   -0.75, 0.0  ,  -0.0090,2, 0, 2,-2);
  Pert.AddSol (   -0.459,   -0.38, 0.0  ,  -0.0053,2, 0,-2, 2);
  Pert.AddSol (   -1.298,    0.74, 0.0  ,  +0.0004,2, 0,-2, 0);
  Pert.AddSol (    0.538,    1.14, 0.0  ,  -0.0141,2, 0,-2,-2);
  Pert.AddSol (    0.263,    0.02, 0.0  ,   0.0   ,1, 1, 2, 0);
  Pert.AddSol (    0.426,   +0.07, 0.0  ,  -0.0006,1, 1,-2,-2);
  Pert.AddSol (   -0.304,   +0.03, 0.0  ,  +0.0003,1,-1, 2, 0);
  Pert.AddSol (   -0.372,   -0.19, 0.0  ,  -0.0027,1,-1,-2, 2);
  Pert.AddSol (   +0.418,    0.0 , 0.0  ,   0.0   ,0, 0, 4, 0);
  Pert.AddSol (   -0.330,   -0.04, 0.0  ,   0.0   ,3, 0, 2, 0);


  // Solar perturbations in latitude
  Pert.AddN(-526.069, 0, 0,1,-2);   Pert.AddN(  -3.352, 0, 0,1,-4);
  Pert.AddN( +44.297,+1, 0,1,-2);   Pert.AddN(  -6.000,+1, 0,1,-4);
  Pert.AddN( +20.599,-1, 0,1, 0);   Pert.AddN( -30.598,-1, 0,1,-2);
  Pert.AddN( -24.649,-2, 0,1, 0);   Pert.AddN(  -2.000,-2, 0,1,-2);
  Pert.AddN( -22.571, 0,+1,1,-2);   Pert.AddN( +10.985, 0,-1,1,-2);


  // Planetary perturbations
  Pert.Planetary(T); 

  // Position vector
  return Vec3D ( Polar(Pert.lambda(),Pert.beta(),Pert.dist()) );
}


//------------------------------------------------------------------------------
//
// MoonEqu: Computes the Moon's equatorial position using Brown's theory
//          (Improved Lunar Ephemeris)
//
// Input:
//
//   T         Time in Julian centuries since J2000
//
// <return>:   Geocentric equatorial position of the Moon (in [km])
//             referred to the true equinox of date
//
// Notes: Light-time is already taken into account
//
//------------------------------------------------------------------------------
Vec3D MoonEqu (double T)
{
  return NutMatrix(T) * Ecl2EquMatrix(T) * MoonPos(T);
}


//------------------------------------------------------------------------------
//
// MiniMoon: Computes the Moon's Ra and declination using a low precision 
//           analytical series
//
// Input:
//
//   T         Time in Julian centuries since J2000
//
// Output:
//
//   Ra        Right Ascension of the Moon in [rad]
//   Dec       Declination of the Moon in [rad]
//
//------------------------------------------------------------------------------
void MiniMoon (double T, double& Ra, double& Dec)
{
  //
  // Constants
  //
  const double eps = 23.43929111*Rad; 

  
  //
  // Variables
  //
  double L_0, l,ls, F, D, dL, S, h, N;
  double l_Moon, b_Moon;
  Vec3D  e_Moon;
  
  
  // Mean elements of lunar orbit
  L_0 = Frac (0.606433 + 1336.855225*T);       // mean longitude [rev]
  
  l  = pi2*Frac ( 0.374897 + 1325.552410*T );  // Moon's mean anomaly 
  ls = pi2*Frac ( 0.993133 +   99.997361*T );  // Sun's mean anomaly 
  D  = pi2*Frac ( 0.827361 + 1236.853086*T );  // Diff. long. Moon-Sun 
  F  = pi2*Frac ( 0.259086 + 1342.227825*T );  // Dist. from ascending node 
    
  
  // Perturbations in longitude and latitude
  dL = +22640*sin(l) - 4586*sin(l-2*D) + 2370*sin(2*D) +  769*sin(2*l) 
       -668*sin(ls) - 412*sin(2*F) - 212*sin(2*l-2*D) - 206*sin(l+ls-2*D)
       +192*sin(l+2*D) - 165*sin(ls-2*D) - 125*sin(D) - 110*sin(l+ls)
       +148*sin(l-ls) - 55*sin(2*F-2*D);
  S  = F + (dL+412*sin(2*F)+541*sin(ls)) / Arcs; 
  h  = F-2*D;
  N  = -526*sin(h) + 44*sin(l+h) - 31*sin(-l+h) - 23*sin(ls+h) 
       + 11*sin(-ls+h) - 25*sin(-2*l+F) + 21*sin(-l+F);
  

  // Ecliptic longitude and latitude
  l_Moon = pi2 * Frac( L_0 + dL/1296.0e3 ); // [rad]
  b_Moon = ( 18520.0*sin(S) + N ) / Arcs;   // [rad]
    

  // Equatorial coordinates
  e_Moon = R_x(-eps) * Vec3D(Polar(l_Moon,b_Moon));
  
  Ra  = e_Moon[phi];
  Dec = e_Moon[theta];
}
