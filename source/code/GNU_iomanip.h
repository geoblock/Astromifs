//------------------------------------------------------------------------------
//
// GNU_iomanip.h
//
//    Temporaray implementation of ostream manipulators from the 
//    C++ Standard Library, which are not contained in <iomanip> 
//    as provided with GNU C++.
//
//------------------------------------------------------------------------------

#include <iomanip>
#include <iostream>

namespace{
ostream& left (ostream& os){os.setf(ios::left ,ios::adjustfield); return os;};
ostream& right(ostream& os){os.setf(ios::right,ios::adjustfield); return os;};
ostream& fixed(ostream& os){os.setf(ios::fixed,ios::floatfield);  return os;};
ostream& showpos  (ostream& os){os.setf(ios::showpos); return os;};
ostream& noshowpos(ostream& os){os.unsetf(ios::showpos); return os;};
}
