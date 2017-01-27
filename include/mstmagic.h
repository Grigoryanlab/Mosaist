#ifndef _MSTMAGIC_H
#define _MSTMAGIC_H

#include "msttypes.h"
#include <algorithm>

using namespace MST;

class TERMUtils {
  public:
    static vector<AtomPointerVector> mostDesignableFragments(Structure& C, vector<Structure*>& TERMs, int n = 10, CartesianPoint* cen = NULL, CartesianPoint* ext = NULL, string outBase = "", bool verb = false);

};

#endif
