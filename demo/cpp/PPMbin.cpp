//------------------------------------------------------------------------------
//
// PPMbin.cpp
// 
// Conversion of the PPM star catalogue (ASCII->binary)
//
//------------------------------------------------------------------------------

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "APC_Const.h"
#include "APC_Math.h"

#ifdef __GNUC__  // GNU C++ Adaptations
#include "GNU_iomanip.h"
#endif

using namespace std;


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{

  // Variables

  ifstream  inp;        // Input file
  ofstream  out;        // Output file
  int       count;      // Counter

  long      PPM;        // PPM number
  char      DM[10];     // Bonner/Cordoba Durchmusterung number
  float     mag;        // Brightness 
  char      Sp[3];      // Spectral type
  int       h,m;        // Right ascension (hours,minutes)
  double    s;          // Right ascension (seconds)
  char      sign;       // Declination (sign)
  int       deg,min;    // Declination (degrees, minutes)
  double    sec;        // Declination (seconds)
  double    Ra,Dec;     // Right ascesnion, declination
  float     pmRA,pmDec; // Proper motion 

  struct {
    long   PPM;
    double Ra,Dec;
    float  pmRA,pmDec,mag;
  } Data;

  // Titel
  cout << endl
       << "     PPMbin: Conversion of the PPM Catalogue (ASCII->binary) "
       << endl
       << "          (c) 1999 Oliver Montenbruck, Thomas Pfleger       " 
       << endl << endl;

  // Open input file for reading
  if (argc>1) 
    inp.open(argv[1]);
  else {
    inp.open("PPM.dat");
  };
  if (!inp) {
    cerr << " Abort. Input file not found." << endl;
    exit(1);
  }

  // Open output file
  if (argc>2) 
    out.open(argv[2],ios::binary|ios::out);
  else
    out.open("PPM.bin",ios::binary|ios::out);
  if (!out) {
    cerr << " Abort. Error while opening output file." << endl;
    exit(1);
  }
  cout << endl << " ";
    
  // Read all stars

  count = 0;
  
  while (true) {

    // Read one input line
                             //  Spalten  Inhalt
    inp.ignore(1);           //   1       (empty)
    inp >> PPM;              //   2-  7   PPM number
    inp.ignore(2);           //   8-  9   (empty)
    inp.get(DM,10);          //  10- 18   Bonner/Cordoba Durchmusterung number
    inp.ignore(1);           //  19       (empty)
    inp >> mag;              //  20- 23   Brightness
    inp.ignore(1);           //  24       (empty)
    inp.get(Sp,3);           //  25- 26   Spectral type
    inp.ignore(1);           //  27       (empty)
    inp >> h;                //  28- 29   Right ascension J2000 (hours)
    inp.ignore(1);           //  30       (empty)
    inp >> m;                //  31- 32   Right ascension J2000 (minutes)
    inp.ignore(1);           //  33       (empty)
    inp >> s;                //  34- 39   Right ascension J2000 (seconds)
    inp.ignore(2);           //  40- 41   (empty)
    inp.get(sign);           //  42       Declination J2000 (sign)
    inp >> deg;              //  43- 44   Declination J2000 (hours)
    inp.ignore(1);           //  45       (empty)
    inp >> min;              //  46- 47   Declination J2000 (minutes)
    inp.ignore(1);           //  48       (empty)
    inp >> sec;              //  49- 53   Declination J2000 (seconds)
    inp.ignore(2);           //  54- 55   (empty)
    inp >> pmRA;             //  56- 62   Proper motion in Ra [s/year]
    inp.ignore(1);           //  63       (empty)
    inp >> pmDec;            //  64- 69   Proper motion in Dec ["/year]
    
    inp.ignore(133,'\n'); 

    if (inp.fail()) break;


    // Coordinates epoch J2000  [rad]
    
    Ra  = Rad*15.0*Ddd(h,m,s);
    Dec = Rad*Ddd(deg,min,sec); if (sign=='-') Dec = -Dec;

    // Proper motion [rad/cent.]

    pmRA  = Rad*15.0*(pmRA*100.0)/3600.0;
    pmDec = Rad*(pmDec*100.0)/3600.0;

    // Output (binary format)

    Data.PPM   = PPM;
    Data.Ra    = Ra;
    Data.Dec   = Dec;
    Data.pmRA  = pmRA;
    Data.pmDec = pmDec;
    Data.mag   = mag;

    out.write( (char *) &Data, sizeof Data ); 
    count++;

    if ( PPM % 10000 == 0 ) cout << "*"; cout.flush();

  };

  out.close();
  
  cout << endl << endl
       << " " << count << " stars processed" << endl;
  
}
