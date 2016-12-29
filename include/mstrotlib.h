#ifndef _MSTLIB_H
#define _MSTLIB_H

#include "msttypes.h"

namespace MST {

class RotamerLibrary {
  public:
    RotamerLibrary();
    RotamerLibrary(string rotLibFile);
    ~RotamerLibrary();
  
    void readRotamerLibrary(string rotLibFile);

  private:
    
};

}
#endif
