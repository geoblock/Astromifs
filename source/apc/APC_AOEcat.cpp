//------------------------------------------------------------------------------
//
// Module:     APC_AOEcat.cpp
//
// Selection of Minor Planets from Lowell Observatory's
//             "Asteroids Orbital Elements Database"
//
//------------------------------------------------------------------------------

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>


#ifdef __GNUC__  // GNU C++ Adaptations
#include <cstdlib>
#include <ctype.h>
#include <string.h>
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

  ifstream  inp;                  // Input file
  ofstream  out;                  // Output file
  int       k;                    // Counter
  char      Query[19],query[19];  // Search string
  bool      found,numeric;        // Flags

  char      Num[6];               // Number of the minor planet
  int       n,Number;             // Number of the minor planet
  char      Name[19],name[19];    // Name/designation
  double    H,G;                  // Brightness 
  long      ymd;                  // Epoch
  int       year,month,day;       // Epoch
  double    a,e,i,Omega,omega,M;  // Orbital elements


  // Header
  
  cout << endl
       << "  AOEcat: Selection of Data from the "
       << "Asteroid Orbital Elements Database"
       << endl
       << "               (c) 1999 Oliver Montenbruck, Thomas Pfleger" 
       << endl << endl;
  
  // Open input file for reading
  
  if (argc>1) 
    inp.open(argv[1]);
  else {
    inp.open("astorb.dat");
  };
  if (!inp) {
    cerr << " Abort. Input file not found." << endl;
    exit(1);
  }

  // Query
   
  cout << " Name/number of the minor planet: ";
  cin.get(Query,19,'\n'); 
  cout << endl;

  strcpy(query,Query);
  for (k=0;query[k]!=0;k++) query[k]=tolower(query[k]);  // Lower case copy

  numeric = true;
  for (k=0;query[k]!=0;k++) numeric = numeric & (isdigit(Query[k])!=0);
  n = atoi(Query);


  // Redirect output if output file shall be created
  if (argc>2) {
    out.open(argv[2]);
	if (out.is_open()) cout = out;
  };
  cout << endl;
    
  // Read all data sets

  while (true) {

    // Lesen einer Zeile
                                //  Columns  Contents
    inp.get(Num,6);             //   1-  5   Number of the minor planet
    inp.ignore(1);              //   6       (empty)
    inp.get(Name,19);           //   7- 24   Name/designation
    inp.ignore(1);              //  25       (empty)
    inp.ignore(15);             //  26- 40   Author  
    inp.ignore(1);              //  41       (empty)
    inp >> H;                   //  42- 46   Brightness [mag]
    inp.ignore(1);              //  47       (empty)
    inp >> G;                   //  48- 52   Slope parameter [mag]
    inp.ignore(53);             //  53-105   (unused)
    inp >> ymd;                 // 106-113   Epoch (YearMonthDay)
    inp.ignore(1);              // 114       (empty)
    inp >> M;                   // 115-124   Mean anomaly
    inp.ignore(1);              // 125       (empty)
    inp >> omega;               // 126-135   Argument of latitude
    inp.ignore(1);              // 136       (empty)
    inp >> Omega;               // 137-146   Longitude of the ascending node
    inp.ignore(1);              // 147       (empty)
    inp >> i;                   // 148-156   Inclination
    inp.ignore(1);              // 157       (empty)
    inp >> e;                   // 158-167   Eccentricity
    inp.ignore(1);              // 168       (empty)
    inp >> a;                   // 169-180   Semimajor axis

    // Extracting data

    Number = atoi(Num);
    
    year  = ymd/10000;
    month = (ymd%10000)/100;
    day   = ymd%100;

    if (inp.fail()) break;

    // Test
    
    strcpy(name,Name);
    for (k=0;name[k]!=0;k++) name[k]=tolower(name[k]); // Lower case copy
    
    found = ( strncmp(query,name,strlen(query)) == 0 );

    found = ( found | (numeric & (n==Number)) );
    
    // Ausgabe

    if (found) {

      cout << setw(5) << year  << " "  << setfill('0')
           << setw(2)  << month << " " << setw(2)  << day << ".0"
           << setfill(' ') 
           << "   ! Epoch        Minor planet";
      if (Number!=0) cout << "(" << Number << ") ";
      cout << Name << endl
           << fixed
           << setprecision(8) << setw(13) << a
           << "   ! a [AU]" << endl
           << setprecision(8) << setw(13) << e
           << "   ! e" << endl
           << setprecision(6) << setw(13) << i
           << "   ! i [deg]" << endl
           << setprecision(6) << setw(13) << Omega
           << "   ! Omega [deg]" << endl
           << setprecision(6) << setw(13) << omega
           << "   ! omega [deg]" << endl
           << setprecision(6) << setw(13) << M
           << "   ! M [deg]" << endl
           << setprecision(1) << setw(13) << 2000.0 
           << "   ! Equinox" << endl;

      break;

    };
    
    // Skip remaining characters
    
    inp.ignore(266,'\n'); 

  };

  // Error message

  if (!found) {
    cout << "Minor planet " << Query << " not found" << endl;
  }
  
  // Close input file

  inp.close();

}
