//------------------------------------------------------------------------------
//
// APC_Math.cpp
//
//    Implements mathematical functions and classes
//
//------------------------------------------------------------------------------

#include <cmath>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <vcl.h>
#include <windows.h>    // Header file for windows
#include <stdio.h>      // Header file for standard Input/Output


#include "APC_Math.h"

#ifdef __GNUC__  // GNU C++ adaptation
#include "GNU_iomanip.h"
#endif

using namespace std;



//------------------------------------------------------------------------------
//
// Frac: Gives the fractional part of a number
//
//------------------------------------------------------------------------------

double Frac (double x)
{
   return x-floor(x);
}


//------------------------------------------------------------------------------
//
// Modulo: calculates x mod y
//
//------------------------------------------------------------------------------

double Modulo (double x, double y)
{
   return y*Frac(x/y);
}


//------------------------------------------------------------------------------
//
// AddThe: Calculates cos(alpha+beta) and sin(alpha+beta) using addition 
//         theorems
//
// Input:
//
//   c1        cos(alpha)
//   s1        sin(alpha)
//   c2        cos(beta)
//   s2        sin(beta)
//
// Output:
//
//   c         cos(alpha+beta)
//   s         sin(alpha+beta)
//
//------------------------------------------------------------------------------
void AddThe ( double c1, double s1, double c2, double s2, 
              double& c, double& s )
{
  c = c1 * c2 - s1 * s2; 
  s = s1 * c2 + c1 * s2;
}


//------------------------------------------------------------------------------
//
// Ddd: Conversion of angular degrees, minutes and seconds of arc to decimal
//   representation of an angle 
//
// Input:
//
//   D        Angular degrees
//   M        Minutes of arc
//   S        Seconds of arc
//
// <return>:  Angle in decimal representation
//
//------------------------------------------------------------------------------
double Ddd(int D, int M, double S)
{
  //
  // Variables
  //
  double sign;


  if ( (D<0) || (M<0) || (S<0) ) sign = -1.0; else sign = 1.0;
    
  return  sign * ( fabs(D)+fabs(M)/60.0+fabs(S)/3600.0 );
}


//------------------------------------------------------------------------------
//
// DMS: Finds degrees, minutes and seconds of arc for a given angle 
//
// Input:
//
//   Dd        Angle in degrees in decimal representation
//
// Output:
//
//   D         Angular degrees
//   M         Minutes of arc
//   S         Seconds of arc
//
//------------------------------------------------------------------------------
void DMS (double Dd, int& D, int& M, double& S)
{
  //
  // Variables
  //
  double x;


  x = fabs(Dd);   D = int(x);
  x = (x-D)*60.0; M = int(x);  S = (x-M)*60.0;
  
  if (Dd<0.0) { if (D!=0) D*=-1; else if (M!=0) M*=-1; else S*=-1.0; }
}


//------------------------------------------------------------------------------
//
// Quad: Quadratic interpolation
//
//       Performs root finding and search for extreme values based on three 
//       equidistant values of a function.
//
// Input:
//
//   y_minus   Value of function at x = -1
//   y_0       Value of function at x =  0
//   y_plus    Value of function at x =  1
//
// Output:
//
//   xe        Abscissa of extremum (may be outside [-1, 1])
//   ye        Value of function at xe
//   root1     First root found
//   root2     Second root found
//   n_root    Number of roots found in [-1, 1]
//
// Notes:
//   
//   Roots and location of an extremum are calculated from a parabola through
//   the given functional values.
//   To investigate functional values for other abscissas than -1, 0 and 1,
//   a linear transformation of given values can be applied, if the given 
//   abscissae are spaced equally.
//
//------------------------------------------------------------------------------
void Quad ( double y_minus, double y_0, double y_plus,
            double& xe, double& ye, double& root1, 
            double& root2, int& n_root )
{
  //
  // Variables
  //
  double a,b,c, dis, dx;


  n_root = 0;

  // Coefficients of interpolating parabola y=a*x^2+b*x+c
  a  = 0.5*(y_plus+y_minus) - y_0;
  b  = 0.5*(y_plus-y_minus);
  c  = y_0;

  // Find extreme value
  xe = -b/(2.0*a); 
  ye = (a*xe+b) * xe + c;

  dis = b*b - 4.0*a*c; // Discriminant of y=a*x^2+b*x+c

  if (dis >= 0) // Parabola has roots 
  {
    dx = 0.5 * sqrt (dis) / fabs (a);

    root1 = xe - dx; 
    root2 = xe + dx;

    if (fabs(root1) <= 1.0) ++n_root;  
    if (fabs(root2) <= 1.0) ++n_root;
    if (root1       < -1.0) root1 = root2;
  }
}


//------------------------------------------------------------------------------
//
// Pegasus: Root finder using the Pegasus method
//
// Input:
//
//   PegasusFunct  Pointer to the function to be examined
//
//   LowerBound    Lower bound of search interval
//   UpperBound    Upper bound of search interval
//   Accuracy      Desired accuracy for the root
//
// Output:
//
//   Root          Root found (valid only if Success is true)
//   Success       Flag indicating success of the routine
//
// References:                                                               
//
//   Dowell M., Jarratt P., 'A modified Regula Falsi Method for Computing    
//     the root of an equation', BIT 11, p.168-174 (1971).                   
//   Dowell M., Jarratt P., 'The "PEGASUS Method for Computing the root      
//     of an equation', BIT 12, p.503-508 (1972).                            
//   G.Engeln-Muellges, F.Reutter, 'Formelsammlung zur Numerischen           
//     Mathematik mit FORTRAN77-Programmen', Bibliogr. Institut,             
//     Zuerich (1986).                                                       
//
// Notes:
//
//   Pegasus assumes that the root to be found is bracketed in the interval
//   [LowerBound, UpperBound]. Ordinates for these abscissae must therefore
//   have different signs.
//
//------------------------------------------------------------------------------
void Pegasus ( PegasusFunct f,
               double LowerBound, double UpperBound, double  Accuracy,
               double& Root, bool& Success )
{
  //
  // Constants
  //
  const int MaxIterat = 30; 

  //
  // Variables
  //
  double x1 = LowerBound; double f1 = f(x1);
  double x2 = UpperBound; double f2 = f(x2);
  double x3 = 0.0;        double f3 = 0.0;
  
  int Iterat = 0; 


  // Initialization
  Success = false;
  Root    = x1;

  
  // Iteration
  if ( f1 * f2 < 0.0 )
    do 
    {
      // Approximation of the root by interpolation
      x3 = x2 - f2/( (f2-f1)/(x2-x1) ); f3 = f(x3);

      // Replace (x1,f2) and (x2,f2) by new values, such that
      // the root is again within the interval [x1,x2]
      if ( f3 * f2 <= 0.0 ) {
        // Root in [x2,x3]
        x1 = x2; f1 = f2; // Replace (x1,f1) by (x2,f2)
        x2 = x3; f2 = f3; // Replace (x2,f2) by (x3,f3)
      }
      else {
        // Root in [x1,x3]
        f1 = f1 * f2/(f2+f3); // Replace (x1,f1) by (x1,f1')
        x2 = x3; f2 = f3;     // Replace (x2,f2) by (x3,f3)
      }

      if (fabs(f1) < fabs(f2))
        Root = x1;
      else
        Root = x2;

      Success = (fabs(x2-x1) <= Accuracy);
      Iterat++;
    }
    while ( !Success && (Iterat<MaxIterat) );
}



//
// Implementation of class Angle for angle output
//

//
// Constructor
//
Angle::Angle (double alpha, AngleFormat Format)
 : m_angle(alpha),
   m_Format(Format)
{
}


//
// Modifiers
//
void Angle::Set (AngleFormat Format)
{
  m_Format = Format;
}


//
// Angle output
//
ostream& operator << (ostream& os, const Angle& alpha)
{
  //
  // Variables
  //
  int    D,M;       // degrees, minutes
  double S;         // seconds
  long   flags;     // output stream format flags
  int    w,p,sp;    // width, precision, showpos
  double fac;       // factor for rounding
  double DegRound;  // rounded value (positive)

    
  flags = os.flags(); // output stream format flags

  w  = os.width();                 // total field width
  p  = os.precision();             // precision (digits)
  sp = flags & ios::showpos;       // flag for positive sign

  os << noshowpos << right;

  switch (alpha.m_Format) {

    case Dd:      // decimal output

      os << fixed;
      os << alpha.m_angle;
      
      break;

    case DMM:     // round to 1 arcmin

      DegRound = floor(60.0*fabs(alpha.m_angle)+0.5)/60.0+0.00001;

      DMS (DegRound, D, M, S);

      if (alpha.m_angle<0.0) 
        os << setw(1) << "-" << setw(w-4);
      else if (sp)
        os << setw(1) << "+" << setw(w-4);
      else 
        os << setw(w-3);
      os  << D << " "  << setfill('0') << setw(2) << M << setfill(' ');

      break;

    case DMMm:    // round to output precision (10^-p arcmins)
      
      fac = 60.0*pow(10.0,p);

      DegRound = floor(fac*fabs(alpha.m_angle)+0.5)/fac+0.1/fac;

      DMS (DegRound, D, M,S);

      if (alpha.m_angle<0.0) 
        os << setw(1) << "-" << setw(w-5-p);
      else if (sp)
        os << setw(1) << "+" << setw(w-5-p);
      else
        os << setw(w-4-p);
      os << D << " "  << setfill('0') << fixed 
         << setw(3+p) << M+S/60.0 << setfill(' ');   

      break;

    case DMMSS:   // round to 1 arcsec

      DegRound = floor(3600.0*fabs(alpha.m_angle)+0.5)/3600.0+0.00001;

      DMS (DegRound, D, M, S);

      if (alpha.m_angle<0.0) 
        os << setw(1) << "-" << setw(w-7);
      else if (sp)
        os << setw(1) << "+" << setw(w-7);
      else 
        os << setw(w-6);
      os  << D << " "  << setfill('0') << setw(2) 
          << M << " " << setw(2) << int(S) << setfill(' ');

      break;

    case DMMSSs:  // round to output precision (10^-p arcsecs)
      
      fac = 3600.0*pow(10.0,p);

      DegRound = floor(fac*fabs(alpha.m_angle)+0.5)/fac+0.1/fac;

      DMS (DegRound, D, M, S);

      if (alpha.m_angle<0.0) 
        os << setw(1) << "-" << setw(w-8-p);
      else if (sp)
        os << setw(1) << "+" << setw(w-8-p);
      else
        os << setw(w-7-p);
      os << D << " "  << setfill('0') << setw(2) << M << " "
         << fixed << setw(3+p) << S << setfill(' ');

      break;
  }
   
  /// os.setf(flags); // restore output stream format flags

  return os;
}


//
// Implementation of least squares estimation class SolverLSQ
//

//
// Constructor; n: Number of estimation parameters
//
SolverLSQ::SolverLSQ (int n) 
  : N(n)
{
  // Allocate storage for R and d
  d = new double[N];
  R = new double*[N];
  
  for (int i=0; i<N; i++) R[i] = new double[N];
  
  Init(); // Reset d and R
}


//
// Destructor
//
SolverLSQ::~SolverLSQ()
{
  // Deallocate d and R
  delete[] d;
  for (int i=0; i<N; i++) delete[] R[i];
  delete[] R;
};


//
// Initialize new problem
//
void SolverLSQ::Init()
{
  // Reset R and d
  for (int i=0; i<N; i++) {
    d[i]=0.0;
    for (int j=0;j<N;j++) R[i][j]=0.0;
  }
}


//
//  Add a data equation of form Ax=b to a least squares problem
//  (performs a row-wise QR transformation using Givens rotations)
//
void SolverLSQ::Accumulate (const double A[], double b)
{
  //
  // Variables
  //
  int      i,j;
  double   c,s,h;
  double*  a = new double[N];


  for (i=0;i<N;i++) a[i]=A[i];   // Copy A
  
  // Construct and apply Givens plane rotation.
  for (i=0; i<N; i++)
  {
    // Construct the rotation and apply it to
    // eliminate the i-th element of a.
    if ( R[i][i]==0.0 &&  a[i]==0.0 ) {
      c = 1.0; s = 0.0; R[i][i] = 0.0;
    }
    else {
      h = sqrt ( R[i][i]*R[i][i] + a[i]*a[i] );
      if (R[i][i]<0.0) h=-h; 
      c = R[i][i]/h;
      s = a[i]/h;
      R[i][i] = h;
    };

    a[i] = 0.0;

    // Apply the rotation to the remaining elements of a
    for (j=i+1; j<N; j++) {
      h       = +c*R[i][j]+s*a[j];
      a[j]    = -s*R[i][j]+c*a[j];
      R[i][j] = h;
    }

    // Apply the rotation to the i-th element of d
    h    = +c*d[i]+s*b;
    b    = -s*d[i]+c*b;
    d[i] = h;
  }

  delete[] a;
}


//
// Solve the LSQ problem for vector x[] by backsubstitution
//
void SolverLSQ::Solve (double x[])
{
  //
  // Variables
  //
  int i,j; i=j=0;
  double Sum=0.0;


  // Check for singular matrix 
  for (i=0;i<N;i++)
    if ( R[i][i] == 0.0 ) {
      cerr << " ERROR: Singular matrix R in SolverLSQ::Solve()" << endl;
      exit(1);
    };

  //  Solve Rx=d for x_n,...,x_1 by backsubstitution
  x[N-1] = d[N-1] / R[N-1][N-1];
  for (i=N-2;i>=0;i--) {
    Sum = 0.0;
    for (j=i+1;j<N;j++) Sum += R[i][j]*x[j];
    x[i] = ( d[i] - Sum ) / R[i][i];
  };
}
