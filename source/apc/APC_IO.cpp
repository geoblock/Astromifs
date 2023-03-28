//------------------------------------------------------------------------------
//
// APC_IO.cpp
// 
//    Implements file I/O utilities
//
//------------------------------------------------------------------------------

#include <iostream>
#include <fstream>     

#ifdef __GNUC__   // GNU C++ adaptation
#include <memory>       
#endif

#include "APC_IO.h"


using namespace std;


//------------------------------------------------------------------------------
//
// FileExists: Checks for existence of file
//
// Input:
//
//   Filename   Name of file to check for
//
//  <return>    Flag indicating existence of file
//
//------------------------------------------------------------------------------
bool FileExists(char* Filename)
{
  //
  // Variables
  //
  
  bool exists = false;

  if (Filename != NULL) {
    ifstream File;                 
    File.open(Filename,ios::in);   
    exists = File.is_open();  
    if (exists) File.close();  
  }
  
  return exists;
}


//------------------------------------------------------------------------------
//
// GetFilenames: Evaluates command line arguments and checks for file names
//
// Input:
//
//   argc           Argument count as defined by main()
//   argv           Argument vector as defined by main()
//   DefInputFile   Filename to be used when no command line argument is given
//
// Output:
//
//   InputFile        Name of input file if found
//   FoundInputFile   Flag indicating whether input file could be found
//   OutputFile       Name of output file if found
//   FoundOutputFile  Flag indicating whether input file could be found
//
// Notes: file names are limited to max. APC_MaxFilename characters
//
//------------------------------------------------------------------------------
void GetFilenames( int argc, char* argv[], 
                   char* DefInputFile, 
                   char* InputFile,  bool& FoundInputFile,
                   char* OutputFile, bool& FoundOutputFile )
{
  FoundInputFile  = false;
  FoundOutputFile = false;

  // Clear buffers
  memset(InputFile,  0, APC_MaxFilename);
  memset(OutputFile, 0, APC_MaxFilename);

  
  // Check for command line arguments
  if (argc>=2) {  
    if( FileExists(argv[1]) ) {
      strncpy(InputFile, argv[1], APC_MaxFilename-1);
      FoundInputFile = true;
      if (argc>2) {
        strncpy(OutputFile, argv[2], APC_MaxFilename-1);
        FoundOutputFile = true;
      }
    }
    else cout << " File " << argv[1] << " not found." << endl;
  }
  else {  // no command line argument given; using default filename
    if ( FileExists(DefInputFile) ) {
      strncpy(InputFile, DefInputFile, APC_MaxFilename-1);
      FoundInputFile = true;
      cout << " Using default input file " << DefInputFile << endl;
    }
    else 
      cout << " Default file " << DefInputFile 
           << " not found." << endl;
  }
}
