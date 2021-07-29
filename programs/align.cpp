#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>
#include <typeinfo>
#include <stdio.h>
#include "msttypes.h"
#include "mstoptions.h"
#include "mstlinalg.h"
#include "msttransforms.h"

using namespace std;
using namespace MST;

int main(int argc, char** argv) {
  MstOptions op;
  op.setTitle("Centers and aligns the structure by principal components. Options:");
  op.addOption("p", "input PDB file.", true);
  op.addOption("o", "output PDB file.", true);
  op.addOption("c", "center only, do not align.");
  op.addOption("a", "align only, do not center.");
  op.addOption("ord", "order of axes in the alignment. By default, the largest principal component will align along Z, the next largest along Y, and then X. This can be changed by specifying this. E.g., --o '1 2 3' will align the primary axis along X, secondary along Y, and the last one along Z.");
  op.setOptions(argc, argv);
  Structure S(op.getString("p"));

  if (!op.isGiven("a")) AtomPointerVector(S.getAtoms()).center();
  if (!op.isGiven("c")) {
    CartesianPoint O = AtomPointerVector(S.getAtoms()).getGeometricCenter();
    Matrix axes = MstLinAlg::getPrincipalAxes(S.getAtoms());
    vector<int> axisOrder = {3, 2, 1};
    if (op.isGiven("ord")) {
      vector<int> ord = MstUtils::splitToInt(MstUtils::trim(op.getString("ord")));
      if (ord.size() != 3) MstUtils::error("--ord must specify a quoted list of three values");
      if (!(MstUtils::setdiff(axisOrder, ord).empty() && MstUtils::setdiff(ord, axisOrder).empty())) MstUtils::error("--ord must specify have 1, 2, and 3 in some order.");
      axisOrder = ord;
    }
    CartesianPoint X = CartesianPoint((vector<mstreal>) axes.column(axisOrder[0] - 1)).getUnit();
    CartesianPoint Y = CartesianPoint((vector<mstreal>) axes.column(axisOrder[1] - 1)).getUnit();
    CartesianPoint Z = CartesianPoint((vector<mstreal>) axes.column(axisOrder[2] - 1)).getUnit();
    // make sure system is right handed
    if ((X.cross(Y)).dot(Z) < 0) X = -X;
    // align along axes by switching frames
    Transform T = TransformFactory::switchFrames(Frame(O, CartesianPoint(1, 0, 0), CartesianPoint(0, 1, 0), CartesianPoint(0, 0, 1)), Frame(O, X, Y, Z));
    T.apply(S);
  }

  S.writePDB(op.getString("o"));
  return 0;
}
