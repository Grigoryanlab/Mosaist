#include "mstgeometry.h"

MstGeometry::MstGeometry() {
  for (int i = 0; i < 3; i++) {
    trans[i] = 0;
    for (int j = 0; j < 3; j++) {
      rot[i][j] = 0;
    }
  }
}
