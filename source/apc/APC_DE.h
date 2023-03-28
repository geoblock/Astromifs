//------------------------------------------------------------------------------
//
// APC_DE.h
// 
//    ODE solver 'DE' developed by Gordon and Shampine
//
// Reference:
//
//   Shampine, Gordon: "Computer solution of Ordinary Differential Equations",
//   Freeman and Comp., San Francisco (1975)
//
//------------------------------------------------------------------------------

#ifndef INC_APC_DE_H
#define INC_APC_DE_H


//
// Types
//

// DEfunct prototype: void f (double x, double y[], double yp[])
typedef void (*DEfunct)(
  double x,    // independent variable
  double y[],  // vector with components of function
  double yp[]  // vector with components of derivative w.r.t. x
);


// State codes for the SolverDE integrator
enum DE_STATE {
  DE_INIT                  = 1,    // Starting step
  DE_DONE                  = 2,    // Successful integration step
  DE_ACCURACY_NOT_ACHIEVED = 3,    // Too stringent accuracy requirements
  DE_TOO_MANY_STEPS        = 4,    // Too many steps required
  DE_STIFF                 = 5,    // Suspect of stiff differential equation
  DE_INVALID_PARAMS        = 6     // Invalid input
};


//
// SolverDE integrator class 
//
class SolverDE
{
  public:

    // Constructor, Destructor
    SolverDE (DEfunct F, int Neqn);
    ~SolverDE ();
    
    // Integration
    void Integ ( 
      double y[],       // Solution vector; updated on output
      double& t,        // Value of independent variable; updated on output
      double tout,      // Value of independent variable to integrate to
      double& relerr,   // Desired relative accuracy of the solution; updated
      double& abserr,   // Desired absolute accuracy of the solution; updated
      DE_STATE& State,  // State code
      bool PermitTOUT = true
    );

  private:

    // Data members
    DEfunct  f;
    int      neqn;
    double   *yy,*wt,*p,*yp,*ypout;
    double   **phi;
    double   alpha[13],beta[13],v[13],w[13],psi[13];
    double   sig[14],g[14];
    double   x,h,hold,told,delsgn;
    int      ns,k,kold;
    bool     OldPermit, phase1,start,nornd;
    
    
    // Interpolation
    void Intrp ( double x, const double y[], double xout, 
                 double yout[], double ypout[] );

    // Elementary integration step
    void Step (double& x, double y[], double& eps, bool& crash);
};

#endif   // include blocker
