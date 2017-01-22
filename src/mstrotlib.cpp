#include "mstrotlib.h"

using namespace MST;

RotamerLibrary::RotamerLibrary(string rotLibFile) {
  readRotamerLibrary(rotLibFile);
}

RotamerLibrary::~RotamerLibrary() {
  for (map<string, vector<Residue*> >::iterator it = rotamers.begin(); it != rotamers.end(); ++it) {
    vector<Residue*>& residues = it->second;
    for (int i = 0; i < residues.size(); i++) {
      delete residues[i];
    }
  }
}


void RotamerLibrary::readRotamerLibrary(string rotLibFile) {
  int nc, nb, na, nr;
  Residue* rot;
  real x, y, z, phi, psi;
  float val; // the binary file is written with single precision values, so read as floats and cast as necessary
  fstream inp;
  MstUtils::openFile(inp, rotLibFile, ios_base::in | ios_base::binary, "RotamerLibrary::readRotamerLibrary(string rotLibFile)");

  while (true) {
    string aa = MstUtils::readNullTerminatedString(inp);
    if (inp.eof()) break;
    inp.read((char*) &nc, sizeof(nc)); // number of chi angles in this amino acid
    inp.read((char*) &na, sizeof(na)); // number of atoms in the side-chain of this amino acid
    inp.read((char*) &nb, sizeof(nb)); // number of phi/psi bins defined for this amino acid
    chidef[aa] = vector<vector<string> >(nc);
    chi[aa] = vector<vector<vector<pair<real, real> > > >(nb);
    rotamers[aa] = vector<Residue*>(nb, NULL);

    // read definitions of chi angles
    for (int i = 0; i < nc; i++) {
      chidef[aa][i] = vector<string>(4);
      for (int j = 0; j < 4; j++) {
        chidef[aa][i][j] = MstUtils::readNullTerminatedString(inp);
      }
    }

    // read rotamer atom names and initialize a residue object
    Residue* res = new Residue(aa, 1);
    for (int i = 0; i < na; i++) {
      string atomname = MstUtils::readNullTerminatedString(inp);
      res->appendAtom(new Atom(1, atomname, 0, 0, 0, 0, 0, false, ' '));
    }

    // read phi/psi angles for each bin, record unique ones.
    map<real, bool> uniquePhi, uniquePsi;
    map<pair<real, real>, bool> uniqueBins;
    vector<pair<real, real> > bins;
    for (int i = 0; i < nb; i++) {
      pair<real, real> bin;
      inp.read((char*) &val, sizeof(val));
      phi = angleToStandardRange((real) val);
      uniquePhi[phi] = true;
      bin.first = phi;
      inp.read((char*) &val, sizeof(val));
      psi = angleToStandardRange((real) val);
      uniquePsi[psi] = true;
      bin.second = psi;
      bins.push_back(bin);
      uniqueBins[bin] = true;
    }

    // make sure the phi/psi bins are defined on a grid and bins are not repeated
    if (uniqueBins.size() != uniquePhi.size()*uniquePsi.size()) {
      MstUtils::error("PHI/PSI bins not on a grid for amino acid '" + aa + "' in rotamer library '" + rotLibFile + "'", "RotamerLibrary::readRotamerLibrary(string)");
    }
    if (uniqueBins.size() != nb) {
      MstUtils::error("some PHI/PSI bins are not unique for amino acid '" + aa + "' in rotamer library '" + rotLibFile + "'", "RotamerLibrary::readRotamerLibrary(string)");
    }

    // assign bin centers along phi and psi axes and then standard bin indices for each bin
    binPhiCenters[aa] = keys(uniquePhi, true);
    binPsiCenters[aa] = keys(uniquePsi, true);
    vector<int> binIndex(nb);
    for (int i = 0; i < nb; i++) binIndex[i] = getBackboneBin(aa, bins[i].first, bins[i].second);    

    // read rotamers in each phi/psi bin
    for (int ii = 0; ii < nb; ii++) {
      int i = binIndex[ii];
      inp.read((char*) &nr, sizeof(nr));  // number of rotamers in this bin
      chi[aa][i].resize(nr);
      rot = new Residue(*res);
      for (int j = 0; j < nr; j++) {
        // read chi and chi sigma values for the rotamer
        chi[aa][i][j].resize(nc);
        for (int k = 0; k < nc; k++) {
          inp.read((char*) &val, sizeof(val));
          chi[aa][i][j][k].first = (real) val;
          inp.read((char*) &val, sizeof(val));
          chi[aa][i][j][k].second = (real) val;
        }
        // read atom coordinates for the rotamer
        for (int k = 0; k < na; k++) {
          inp.read((char*) &val, sizeof(val)); x = (real) val;
          inp.read((char*) &val, sizeof(val)); y = (real) val;
          inp.read((char*) &val, sizeof(val)); z = (real) val;
          if (j == 0) {
            (*rot)[k].setCoor(x, y, z);
          } else {
            (*rot)[k].addAlternative(x, y, z, 0, 0, ' ');
          }
        }
      }
      // insert rotamer into the library
      rotamers[aa][i] = rot;
    }
  }
}

real RotamerLibrary::angleToStandardRange(real angle) {
  if ((angle >= -180) && (angle < 180)) return angle;
  return MstUtils::mod(angle + 180, 360.0) - 180;
}

int RotamerLibrary::numberOfRotamers(string aa, real phi, real psi) {
  int bi = getBackboneBin(aa, phi, psi);
  return rotamers[aa][bi]->getAtom(0).numAlternatives();
}

//vector<string> RotamerLibrary::availableAminoAcids() {
//  vector<string> aas;
//  for (map<string, vector<Residue*> >::iterator it = rotamers.begin(); it != rotamers.end(); ++it) {
//    aas.push_back(it->first);
//  }
//  return aas;
//}

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
  if (rotamers[aa].size() <= binInd) {
    MstUtils::error("rotamer library has " + MstUtils::toString(rotamers[aa].size()) + " backbone bins for amino acid '" + aa + "', but bin number " + MstUtils::toString(binInd+1) + " was requested", "RotamerLibrary::placeRotamer");
  }
  Residue& rots = *(rotamers[aa][binInd]);

  // get the transformation to go from the standard position of the backbone, which is the
  // position for which the rotamer library stores side-chain coordinates, to the actual
  // backbone position in the given residue. That's the transformation we will need to apply
  // to go from rotamer-library coordinates to the final placed coordinates.
  CartesianPoint CA = CartesianPoint(res.findAtom("CA"));
  CartesianPoint C = CartesianPoint(res.findAtom("C"));
  CartesianPoint N = CartesianPoint(res.findAtom("N"));
  CartesianPoint X = C - CA;          // X-axis of the residue frame defined to be along CA -> C
  CartesianPoint Z = X.cross(N - CA); // since N-CA is defined to be in the XY plane, (CA -> C) x (CA -> N) will be a vector along Z
  CartesianPoint Y = Z.cross(X);      // finally, Z x X is Y
  Frame R(CA, X, Y, Z);               // residue frame
  Frame L;                            // Frame class defaults to the laboratory frame
  Transform T = TransformFactory::switchFrames(R, L);
  
  // fish out the right rotamer and copy it over to the residue, destroying previous atoms
  vector<Atom*> newAtoms(rots.atomSize(), NULL);
  for (int i = 0; i < rots.atomSize(); i++) {
    Atom& a = rots[i];
    if (a.numAlternatives() < rotIndex) {
      MstUtils::error("rotamer library contains " + MstUtils::toString(a.numAlternatives() + 1) + " rotamers for amino-acid, but rotamer number " + MstUtils::toString(rotIndex+1) + "was requested", "RotamerLibrary::placeRotamer");
    }
    if (rotIndex > 0) a.swapWithAlternative(rotIndex-1);
    Atom* newAtom = new Atom(a);
    T.apply(newAtom);
    newAtoms[i] = newAtom;
    if (rotIndex > 0) a.swapWithAlternative(rotIndex-1);
  }
  vector<int> oldAtoms;
  for (int i = 0; i < res.atomSize(); i++) {
    Atom& a = res[i];
    if (!isBackboneAtom(a.getName()) || ((aa.compare("PRO") == 0) && a.isNamed("H"))) {
      oldAtoms.push_back(i);
    }
  }
  res.replaceAtoms(newAtoms, &oldAtoms);
  res.setName(aa);
}

bool RotamerLibrary::isBackboneAtom(string atomName) {
  /* backbone atoms can be either nitrogens, carbons, oxigens, or hydrogens.
   * specifically, possible known names in each category are:
   * 'N', 'NT'
   * 'CA', 'C', 'CY', 'CAY'
   * 'OY', 'O', 'OCT*', 'OXT', 'OT1', 'OT2'
   * 'H', 'HY*', 'HA1', 'HN', 'HT*'
   */
  if (atomName.size() == 0) return false;
  switch (atomName[0]) {
    case 'N':
      if ((atomName.compare("N") == 0) || (atomName.compare("NT") == 0)) return true;
      return false;
    case 'C':
      if ((atomName.compare("C") == 0) || (atomName.compare("CA") == 0) || (atomName.compare("CY") == 0) || (atomName.compare("CAY") == 0)) return true;
      return false;
    case 'O':
      if ((atomName.compare("O") == 0) || (atomName.compare("OY") == 0) || (atomName.compare("OXT") == 0) || 
          (atomName.compare("OT1") == 0) || (atomName.compare("OT2") == 0) || 
          ((atomName.size() >= 3) && (atomName.compare(0, 3, "OCT") == 0))) return true;
      return false;
    case 'H':
      if ((atomName.compare("H") == 0) || (atomName.compare("HA1") == 0) || (atomName.compare("HN") == 0) || 
          ((atomName.size() >= 2) && ((atomName.compare(0, 2, "HT") == 0) || (atomName.compare(0, 2, "HY") == 0)))) return true;
      return false;
    default:
      return false;
  }
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
  if (dL == 0) {
    return L
  } else if (dU == 0) {
    return U;
  } else if ((dL < 0) && (dU > 0)) {
    // this means the angle maps "before" the smallest angle and "after" the
    // largest angle, meaning that it is sandwiched between the two, so pick
    // the closest one. E.g., the first angle is 10, the last one is -10 (aka
    // 350), and the angle is 5.
    if (abs(dL) < abs(dU)) return L;
    else return U;
  } else if ((dL > 0) && (dU > 0) && (dL > dU)) {
    // this means that the entire set of agles in the array comes after
    // the query angle, so return the lower limit
    return L;
  } else if ((dL < 0) && (dU < 0) && (dL > dU)) {
    // this means that the entire set of agles in the array comes before
    // the query angle, so return the upper limit
    return U;
  }

  // if we have reached here, that means the array has angles both before
  // and after the query angle, so need to find a part of angles to
  // sandwich the query, and then narrow the sandwich
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
