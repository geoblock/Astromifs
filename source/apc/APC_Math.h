//------------------------------------------------------------------------------
//
// APC_Math.h
// 
// Purpose:
//
//    Mathematical functions and classes of general usefulness
//
//------------------------------------------------------------------------------

#ifndef INC_APC_MATH_H
#define INC_APC_MATH_H

#include <iostream>


//
// Types
//

// Format tag for Angle output (used by Angle class output operator)
enum AngleFormat { 
  Dd,    // decimal representation
  DMM,   // degrees and whole minutes of arc
  DMMm,  // degrees and minutes of arc in decimal representation
  DMMSS, // degrees, minutes of arc and whole seconds of arc
  DMMSSs // degrees, minutes, and seconds of arc in decimal representation
};


// Frac: Gives the fractional part of a number
double Frac (double x);


// Modulo: calculates x mod y
double Modulo (double x, double y);


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
              double& c, double& s );


//------------------------------------------------------------------------------
//
// Ddd: Conversion of angular degrees, minutes and seconds of arc to decimal
//      representation of an angle 
//
// Input:
//
//   D         Angular degrees
//   M         Minutes of arc
//   S         Seconds of arc
//
//   <return>  Angle in decimal representation
//
//------------------------------------------------------------------------------
double Ddd(int D, int M, double S);


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
void DMS (double Dd, int& D, int& M, double& S);


//------------------------------------------------------------------------------
//
// Quad: Quadratic interpolation
//
//       Performs root finding and search for extrema based on three 
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
            double& root2, int& n_root );


// Pegasus function prototype: double f(double x);
typedef double (*PegasusFunct) (double x);


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
//
// Notes:
//
//   Pegasus assumes that the root to be found is bracketed in the interval
//   [LowerBound, UpperBound]. Ordinates for these abscissae therefore must
//   have different signs.
//
//------------------------------------------------------------------------------
void Pegasus ( PegasusFunct f,
               double LowerBound, double UpperBound, double  Accuracy,
               double& Root, bool& Success  );


//
// Auxiliary class for sexagesimal angle output
//
class Angle 
{
  public:
  
    // Constructor
    Angle (double alpha, AngleFormat Format=Dd); 
    
    // Modifiers
    void Set (AngleFormat Format=Dd);

    // Angle output
    friend std::ostream& operator << (std::ostream& os, const Angle& alpha);

  private:

    double       m_angle;
    AngleFormat  m_Format;
};


//
// Least squares estimation class SolverLSQ
//
class SolverLSQ
{
  public:

    // Constructor, destructor
    SolverLSQ(int n); // n: Number of estimation parameters
    ~SolverLSQ();

    // Member functions
    
    // Reset to solve a new problem
    void Init();  
    
    // Insert a data equation of form Ax=b to a least squares estimation
    // (performs a row-wise QR transformation using Givens rotations)
    void Accumulate (const double A[], double b);
    
    // Solve the LSQ problem for vector x[] by backsubstitution
    void Solve (double x[]);

  private:

    int      N;         // Number of estimation parameters
    double   *d;        // Right hand side of transformed equations
    double   **R;       // Upper right triangular matrix

};

#endif   // include blocker
