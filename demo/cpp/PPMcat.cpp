//------------------------------------------------------------------------------
//
// PPMcat.cpp
// 
// Selection of stars from the PPM star catalogue
//
//------------------------------------------------------------------------------

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "APC_Const.h"
#include "APC_Math.h"
#include "APC_Spheric.h"

#ifdef __GNUC__  // GNU C++ Adaptation
#include "GNU_iomanip.h"
#endif

using namespace std;


//------------------------------------------------------------------------------
//
// GetInp: Input of the image center and the field of view
//
// Ausgabe:
//
//   RA_0      Right ascension of the image centre [rad]
//   Dec_0     Declination of the image center [rad]
//   Radius    Radius of the field of view [rad]
//   T_Epoch   Epoch of the image (Julian. cent. from J2000)
//
//------------------------------------------------------------------------------
void GetInp (double& RA_0, double& Dec_0, double& Radius, double& T_Epoch)
{
  int     h;
  double  m,d,y;
  
  cout << endl
       << " Image centre" << endl;
  cout << "   Right ascension [hh mm.m] ... ";
  cin >> h >> m; cin.ignore(81,'\n'); RA_0=15.0*Rad*(h+m/60.0);
  cout << "   Declination [deg]         ... ";
  cin >> d; cin.ignore(81,'\n'); Dec_0=Rad*d;
  cout << " Radius [deg]                ... ";
  cin >> d; cin.ignore(81,'\n'); Radius=Rad*d;
  cout << " Epoch [yyyy.y]              ... ";
  cin >> y; cin.ignore(81,'\n'); T_Epoch = (y-2000.0)/100.0;
}

//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{

  // Variablen

  ifstream  inp;        // Input file
  ofstream  out;        // Output file
  
  double    RA0,Dec0;   // Right ascension, declination [rad]
  double    r;          // Field of view radius [rad]
  double    T;          // Epoch

  long      PPM;        // PPM number
  float     mag;        // Brightness 
  double    Ra,Dec;     // Right ascension, declination
  float     pmRA,pmDec; // Proper motion 
  Vec3D     e,e0;       // Unit vectors

  struct {
    long   PPM;
    double Ra,Dec;
    float  pmRA,pmDec,mag;
  } Data;


  // Title
  cout << endl
       << "       PPMcat: Selection of stars from the PPM Catalogue " << endl
       << "          (c) 1999 Oliver Montenbruck, Thomas Pfleger    " << endl 
       << endl;

  // Open input file for reading 
  if (argc>1) 
    inp.open(argv[1],ios::binary|ios::in);
  else {
    inp.open("PPM.bin",ios::binary|ios::in);
  };
  if (!inp) {
    cerr << " Abort. Input file not found." << endl;
    exit(1);
  }

  //User input
  GetInp ( RA0,Dec0,r, T );
  e0 = Vec3D(Polar(RA0,Dec0));

  // Redirect output, if output file shall be written
  if (argc>2) {
    out.open(argv[2]);
    if (out.is_open()) cout = out;
  }

  // Header 
  cout 
   << endl
   << " PPM Number       Ra          Dec        PM(Ra)  PM(Dec)   m   " 
   << endl
   << "               h  m  s        o  '  \"     s/cy    \"/cy    mag  " 
   << endl;
  
 
  // Loop over all stars
 
  while (true) {

    // Read one input line
    inp.read( (char *) &Data, sizeof(Data) );
    if (inp.fail()) break;

    PPM   = Data.PPM;
    Ra    = Data.Ra  + Data.pmRA *T;
    Dec   = Data.Dec + Data.pmDec*T;
    pmRA  = Data.pmRA;
    pmDec = Data.pmDec;
    mag   = Data.mag;

    // Unit vector
    e = Vec3D(Polar(Ra,Dec));

    // Ouput of stars within the field of view
    if ( Dot(e,e0)>cos(r) ) {
      cout << fixed 
           << " PPM " << setw(6) << PPM << "  " 
           << setw(12) << setprecision(3) << Angle(Deg*Ra/15.0,DMMSSs)
           << "  " 
           << showpos
           << setw(12) << setprecision(2) << Angle(Deg*Dec,DMMSSs)
           << setw(8) << setprecision(3) << Deg*pmRA *3600.0/15.0
           << setw(8) << setprecision(2) << Deg*pmDec*3600.0
           << noshowpos
           << setw(6) << setprecision(1) << mag
           << endl;
    };

  };

}
