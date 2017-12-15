#include "msttypes.h"
#include "mstgeometry.h"

using namespace MST;

int main(int argc, char** argv) {
  CartesianGeometry::testPrimitiveGradients();
  RMSDCalculator::testQCP(true);
}
