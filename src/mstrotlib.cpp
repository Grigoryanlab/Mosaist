#include "mstrotlib.h"

using namespace MST;

bool RotamerLibrary::placeRotamer(Residue& res, string aa, int rotIndex, bool strict) {
  double phi = res.getPhi();
  double psi = res.getPsi();
  if (phi == Residue::badDihedral) {
    // SHOULD REALLY ASSUME THE DEFAULT BIN!!!!!!!!!!!!!!!!!!!
    if (strict) MstUtils::error("could not compute PHI for", "RotamerLibrary::placeRotamer");
    return false;
  }
  if (psi == Residue::badDihedral) {
    // SHOULD REALLY ASSUME THE DEFAULT BIN!!!!!!!!!!!!!!!!!!!
    if (strict) MstUtils::error("could not compute PSI for", "RotamerLibrary::placeRotamer");
    return false;
  }

  int binInd = getBackboneBin(aa, phi, psi);
  if (rotamers.find(aa) == rotamers.end()) MstUtils::error("rotamer library does not contain amino acid '" + aa + "'", "RotamerLibrary::placeRotamer");
  if (rotamers[aa].size() <= bindInd) {
    MstUtils::error("rotamer library has " + MstTools:toString(rotamers[aa].size()) + " backbone bins for amino acid '" + aa + "', but bin number " + + MstTools:toString(binIndex+1) + " was requested", "RotamerLibrary::placeRotamer");
  }
  Residue& rots = rotamers[aa][binInd];

  // get the transformation to go from the standard position of the backbone, which is the
  // position for which the rotamer library stores side-chain coordinates, to the actual
  // backbone position in the given residue. That's the transformation we will need to apply
  // to go from rotamer-library coordinates to the final placed coordinates.
  CartesianPoint CA = CartesianPoint(res.findAtom(&res, "CA"));
  CartesianPoint C = CartesianPoint(res.findAtom(&res, "C"));
  CartesianPoint N = CartesianPoint(res.findAtom(&res, "N"));
  CartesianPoint X = C - CA;          // X-axis of the residue frame defined to be along CA -> C
  CartesianPoint Z = X.cross(N - CA); // since N-CA is defined to be in the XY plane, (CA -> C) x (CA -> N) will be a vector along Z
  CartesianPoint Y = Z.cross(X);      // finally, Z x X is Y
  Frame R(CA, X, Y, Z);               // residue frame
  Transform3D T = switchFrames(Frame(CA, X, Y, Z), Frame()); // Frame class defaults to the laboratory frame
  
  // fish out the right rotamer and copy it over to the residue, destroying previous atoms
  vector<Atom*> newAtoms(rots.atomSize(), NULL);
  for (int i = 0; i < rots.atomSize(); i++) {
    Atom& a = rots[i];
    if (a.numAlternatives() < rotIndex) {
      MstUtils::error("rotamer library contains " + MstTools:toString(a.numAlternatives() + 1) + " rotamers for amino-acid, but rotamer number " + MstTools:toString(rotIndex+1) + "was requested", "RotamerLibrary::placeRotamer");
    }
    if (rotIndex > 0) a.swapWithAlternative(rotIndex-1);
    Atom* newAtom = new Atom(a);
    T.appy(newAtom);
    newAtoms[i] = newAtom;
    if (rotIndex > 0) a.swapWithAlternative(rotIndex-1);
  }
  vector<Atom*> oldAtoms;
  for (int i = 0; i < res.atomSize(); i++) {
    Atom& a = res[i];
    if (!isBackbone(a.getName()) || ((aa.compare("PRO") == 0) && a.isNamed("H"))) {
      oldAtoms.push_back(&a);
    }
  }
  res->replaceAtoms(newAtoms, &oldAtoms);
  res->setName(aa);
}

int RotamerLibrary::getBackboneBin(string aa, real phi, real psi) {
  MstUtils::assert(binPhiCenters.find(aa) != binPhiCenters.end(), "no PHI bin for amino acid '" + aa + "'", "RotamerLibrary::getBackboneBin");
  MstUtils::assert(binPsiCenters.find(aa) != binPsiCenters.end(), "no PSI bin for amino acid '" + aa + "'", "RotamerLibrary::getBackboneBin");

  int phiInd = findClosestAngle(binPhiCenters[aa], phi);
  int psiInd = findClosestAngle(binPsiCenters[aa], psi);

  // by convention, linear bin index is encoded in phi-major order
  return phiInd * binPsiCenters[aa].size() + psiInd;
}

int RotamerLibrary::findClosestAngle(vector<real>& array, real value) {
  MstUtils::assert(array.size() != 0, "empty vector passed", "RotamerLibrary::findClosestAngle");

  // binary search on a circle. First determine where the angle maps with
  // respect to the first (smallest) and last (largest) angles in the array,
  // and treat some special cases.
  int L = 0; int U = array.size() - 1;
  double dL = angleDiff(value, array[L]);
  double dU = angleDiff(value, array[U]);
  if ((dL < 0) && (dU > 0)) {
    // this means the angle maps "before" the smallest angle and "after" the
    // largest angle, meaning that it is sandwiched between the two, so pick
    // the closest one. E.g., the first angle is 10, the last one is -10 (aka
    // 350), and the angle is 5.
    if (abs(dL) < abs(dU)) return L;
    else return U;
  } else if ((dL > 0) && (dU > 0)) {
    // if the angle maps after both the first and the last angle in the
    // array, then it is the last angle that is the closest to it
    return U;
  }

  // and now the main component of the binary search
  while (L != U) {
    // treat consequitive upper/lower indices as a special case so that
    // index arithmetic below simple (no special cases to worry about)
    if (U == L + 1) {
      if (abs(angleDiff(value, array[L])) < abs(angleDiff(value, array[U]))) return L;
      else return U;
    }

    int C = (L + U)/2; // try the mid point next (integer division; floor)
    double d = angleDiff(value, array[C]); // value - angle[C]

    if (d < 0) {
      U = C;
    } else if (d > 0) {
      L = C;
    } else {
      return C;
    }
  }
  return L;
}

real RotamerLibrary::angleDiff(real a, real b) {

  real da = MstUtils::mod((MstUtils::mod(a, 360.0) - MstUtils::mod(b, 360.0)), 360.0);
  if (da > 180.0) da -= 360.0;

  return da;
}
