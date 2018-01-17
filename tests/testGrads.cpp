#include "msttypes.h"

using namespace MST;

int main(int argc, char** argv) {
  CartesianGeometry::testPrimitiveGradients(true);
  CartesianGeometry::testPrimitiveGradients(false);
  RMSDCalculator::testQCP(true);
}
