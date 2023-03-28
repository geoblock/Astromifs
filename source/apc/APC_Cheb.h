//------------------------------------------------------------------------------
//
// APC_Cheb.h
// 
//    Declarations for three dimensional Chebyshev approximation
//
//------------------------------------------------------------------------------

#ifndef INC_APC_CHEB_H
#define INC_APC_CHEB_H


#include "APC_VecMat3D.h"


// prototype of function to fit: Vec3D f(double t)
typedef Vec3D (*C3Dfunct)(double t);


// Cheb3D: three dimensional Chebyshev approximation
class Cheb3D
{

  public:
    
    // Constructor
    Cheb3D ( 
       C3Dfunct f,   // Pointer to function to fit 
       int n,        // Degree of fit polynomial
       double dt     // Interval size
    );

    // Destructor
    ~Cheb3D();

    // Fit: approximates function m_f over the interval [ta, tb]
    void Fit (double ta, double tb);

    // Value: evaluation of approximation at argument t,
    //        performs a new fit if necessary
    Vec3D Value (double t);
    
  private:

    // Members
    
    C3Dfunct m_f;             // Function
    int      m_n;             // Degree 
    bool     m_Valid;         // Status flag for validity of Chebyshev coefficients
    double   m_dt;            // Interval size
    double   m_ta,m_tb;       // Interval
    double   *Cx,*Cy,*Cz;     // Chebyshev coefficients
   
};

#endif   // include blocker
