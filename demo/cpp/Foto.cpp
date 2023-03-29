//------------------------------------------------------------------------------
//
// Foto program (Foto.cpp)
//
// Astrometric analysis of photographic plates
//
//------------------------------------------------------------------------------

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "APC_Const.h"
#include "APC_IO.h"
#include "APC_Math.h"
#include "APC_Spheric.h"

#ifdef __GNUC__  // GNU C++ Adaptations
#include "GNU_iomanip.h"
#endif

using namespace std;


//
// Constants
//
const int MaxObj = 30;    // Maximum number of objects on the photo 
const int StrLen = 13;


//------------------------------------------------------------------------------
//
// GetInput: Reads input from file
//
// Input:
//
//   Filename  Name of input file
//
// Output:
//
//   RA0       Right ascension of plate centre
//   Dec0      Declination of plate centre
//   NObj      Number of objects read
//   Name[]    Names of the objects
//   Ra[]      Right ascensions of reference objects
//   Dec[]     Declinations of reference objects
//   x[]       Measured coordinates of reference objects...
//   y[]       ... with respect to the plate centre [mm]
//
// Note: all angles in [rad]
//
//------------------------------------------------------------------------------
void GetInput ( char* Filename,
                double& RA0, double& Dec0, int& NObj, char Name[][StrLen],
                double Ra[], double Dec[], double x[], double y[] )
{
  //
  // Variables
  //
  int       i,k,h,m,deg,min;
  double    s,sec;
  char      c,sign;
  ifstream  inp;

  
  // Open file for reading
  inp.open(Filename); 
  cout << " Input data file: " << Filename << endl << endl;

  // Read coordinates of the plate centre 
  for (k=0;k<StrLen;k++) inp.get(c);

  inp >> h >> m >> s;       inp.get(c).get(c).get(sign);
  inp >> deg >> min >> sec; inp.ignore(81,'\n'); 

  RA0  = Rad*15.0*Ddd(h,m,s);
  Dec0 = Rad*Ddd(deg,min,sec); if (sign=='-') Dec0 = -Dec0;


  // Read name, plate coordinates (and equatorial coordinates) 
  i = -1;
  do {
    ++i;
    inp.get(Name[i],StrLen);  if (inp.fail()) break;
    if ( Name[i][0] == '*' ) { 
      
      // Reference star 
      inp >> x[i] >> y[i];
      inp >> h >> m >> s;       inp.get(c).get(c).get(sign);
      inp >> deg >> min >> sec; inp.ignore(81,'\n'); 

      Ra[i]  = Rad*15.0*Ddd(h,m,s);
      Dec[i] = Rad*Ddd(deg,min,sec); if (sign=='-') Dec[i] = -Dec[i];
    }
    else { 
      
      // Unknown object
      inp >> x[i] >> y[i]; inp.ignore(81,'\n'); 
      Ra[i]=0.0; Dec[i]=0.0;
    };
    
    NObj = i+1;
    c=inp.peek(); // Required to enable EOF check
  }
  while ( !inp.eof() && i<MaxObj );

  inp.close();  // Close input file
}


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  //
  // Constants
  //
  const int NEst = 3;   // Number of estimation parameters

  
  //
  // Variables
  //
  int       i, NObj;
  char      Name[MaxObj][StrLen];
  double    RA0, Dec0, RA_Obs, Dec_Obs, D_RA, D_Dec;
  double    Ra[MaxObj], Dec[MaxObj], Delta[MaxObj];
  double    x[MaxObj], y[MaxObj];
  double    X[MaxObj], Y[MaxObj];
  double    A[NEst],C[NEst];
  double    a, b, c, d, e, f;
  double    Det, FocLen, Scale;
  SolverLSQ LSQx(NEst);
  SolverLSQ LSQy(NEst);
  char      InputFile[APC_MaxFilename] = "";
  char      OutputFile[APC_MaxFilename] = "";
  bool      FoundInputfile = false;
  bool      FoundOutputfile = false;
  ofstream  OutFile;


  // Title
  cout << endl
       << "     FOTO: astrometric analysis of photographic plates  " << endl
       << "           (c) 1999 Oliver Montenbruck, Thomas Pfleger  " << endl
       << endl;


  // Find input file
  GetFilenames( argc, argv, "Foto.dat", InputFile, FoundInputfile,
                OutputFile, FoundOutputfile );
  
  // Terminate program if input file could not be found
  if (!FoundInputfile) {
    cerr << " Terminating program." << endl;
    exit(-1);
  }

  
  // Read input file 
  GetInput (InputFile, RA0,Dec0, NObj, Name, Ra,Dec, x,y);


  // Calculate standard coordinates of reference stars and 
  // fill least squares systems
  for (i=0;i<NObj;i++)
    if (Name[i][0]=='*') {
      EquStd (RA0, Dec0, Ra[i], Dec[i], X[i], Y[i]);
      A[0]=x[i]; A[1]=y[i]; A[2]=1.0;
      LSQx.Accumulate (A, X[i]);
      LSQy.Accumulate (A, Y[i]);
    };

  
  // Calculate plate constants
  LSQx.Solve(C); a=C[0]; b=C[1]; c=C[2];
  LSQy.Solve(C); d=C[0]; e=C[1]; f=C[2];


  // Calculate equatorial coordinates (and errors for reference stars(*) )
  for (i=0;i<NObj;i++) {
    X[i] = a*x[i]+b*y[i]+c;
    Y[i] = d*x[i]+e*y[i]+f;
    StdEqu (RA0, Dec0, X[i], Y[i], RA_Obs, Dec_Obs);
    if (Name[i][0]=='*') {  
      D_RA  = (RA_Obs-Ra[i]) * cos(Dec[i]);
      D_Dec = (Dec_Obs-Dec[i]);
      Delta[i] = Arcs * sqrt(D_RA*D_RA+D_Dec*D_Dec);  // in [arcsec]
    };
    Ra[i] = RA_Obs;  Dec[i] = Dec_Obs;
  };
  
  
  // Calculate focal length and scale
  Det    = a*e-d*b;
  FocLen = 1.0/sqrt(fabs(Det)); // in [mm]
  Scale  = Arcs / FocLen;       // in ["/mm]
  
  

  // Redirect output if output file shall be created
  if (FoundOutputfile) {
    OutFile.open(OutputFile);
	if (OutFile.is_open())
      cout = OutFile;
  }
  
  
  // Output
  cout << " Plate constants:" << endl << endl
       << fixed << setprecision(8)
       << "   a =" << setw(12) << a << "   b =" << setw(12) << b 
       << "   c =" << setw(12) << c << endl
       << "   d =" << setw(12) << d << "   e =" << setw(12) << e 
       << "   f =" << setw(12) << f << endl << endl;

  cout << " Effective focal length and image scale:" << endl << endl
       << fixed << setprecision(2)
       << "   F =" << setw(9) << FocLen << " mm" << endl 
       << "   m =" << setw(7) << Scale  << " \"/mm" << endl << endl;

  cout << " Coordinates:" << endl << endl
       << "        Name        x      y       X       Y "
       << "         Ra           Dec    Error" << endl
       << "                   mm     mm                 "
       << "      h  m  s       o  '  \"     \" " << endl;

  for (i=0;i<NObj;i++) {
    cout << "   " << Name[i] << fixed 
         << setprecision(1) << setw(7) << x[i] << setw(7) << y[i]
         << setprecision(4) << setw(9) << X[i] << setw(8) << Y[i]
         << setprecision(2) << setw(14) << Angle(Deg*Ra[i]/15,DMMSSs)
         << showpos << "  "
         << setprecision(1) << setw(11) << Angle(Deg*Dec[i],DMMSSs)
         << noshowpos;
    if (Name[i][0]=='*') cout << setprecision(1) << setw(6) << Delta[i];
    cout << endl;
  }

  if (OutFile.is_open()) OutFile.close();
}
