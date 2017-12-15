#ifndef _MSTGEOMETRY_H
#define _MSTGEOMETRY_H

/* More complex geometry than in "msttypes.h", which often includes some linear
 * algebra and requires additional #includes and capabilities. */

#include "mstlinalg.h"
#include <chrono>

using namespace MST;

class MstGeometry {
  public:
    MstGeometry();

  protected:
    template <class T>
    static T& ref(T& obj) { return obj; }
    template <class T>
    static T& ref(T* obj) { return *obj; }

  private:
    mstreal rot[3][3];
    mstreal trans[3];
    vector<vector<mstreal> > residuals;
};


#endif
