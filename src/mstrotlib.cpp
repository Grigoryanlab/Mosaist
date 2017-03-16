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
    prob[aa] = vector<vector<real> >(nb);
    binFreq[aa] = vector<real>(nb);

    // read definitions of chi angles
    for (int i = 0; i < nc; i++) {
      chidef[aa][i] = vector<string>(4);
      for (int j = 0; j < 4; j++) {
        chidef[aa][i][j] = MstUtils::readNullTerminatedString(inp);
      }
    }

    // read rotamer atom names and initialize a residue object
    Residue res(aa, 1);
    for (int i = 0; i < na; i++) {
      string atomname = MstUtils::readNullTerminatedString(inp);
      res.appendAtom(new Atom(1, atomname, 0, 0, 0, 0, 0, false, ' '));
    }

    // read phi/psi angles for each bin, record unique ones.
    map<real, bool> uniquePhi, uniquePsi;
    map<pair<real, real>, bool> uniqueBins;
    vector<vector<real> > bins(nb);
    for (int i = 0; i < nb; i++) {
      bins[i] = vector<real>(3);
      pair<real, real> bin;
      inp.read((char*) &val, sizeof(val));
      phi = angleToStandardRange((real) val);
      uniquePhi[phi] = true;
      bins[i][0] = phi;
      inp.read((char*) &val, sizeof(val));
      psi = angleToStandardRange((real) val);
      uniquePsi[psi] = true;
      bins[i][1] = psi;
      uniqueBins[pair<real, real>(phi, psi)] = true;
      inp.read((char*) &val, sizeof(val));
      bins[i][2] = (real) val;
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
    for (int i = 0; i < nb; i++) {
      binIndex[i] = getBackboneBin(aa, bins[i][0], bins[i][1], false);
      binFreq[aa][binIndex[i]] = bins[i][2];
    }

    // make the default phi/psi bin for each amino acid be the most frequent one
    int defaultBinIndex = 0;
    for (int i = 0; i < nb; i++) {
      if (binFreq[aa][i] > binFreq[aa][defaultBinIndex]) defaultBinIndex = i;
    }
    defaultBin[aa] = defaultBinIndex;

    // read rotamers in each phi/psi bin
    for (int ii = 0; ii < nb; ii++) {
      int i = binIndex[ii];
      inp.read((char*) &nr, sizeof(nr));  // number of rotamers in this bin
      chi[aa][i].resize(nr);
      prob[aa][i].resize(nr);
      rot = new Residue(res);
      for (int j = 0; j < nr; j++) {
        // read rotamer probability
        inp.read((char*) &val, sizeof(val));
        prob[aa][i][j] = (real) val;

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

int RotamerLibrary::numberOfRotamers(string aa, real phi, real psi, bool strict) {
  int bi = getBackboneBin(aa, phi, psi, !strict);
  if (rotamers[aa][bi]->atomSize() == 0) return 1; // for GLY
  return rotamers[aa][bi]->getAtom(0).numAlternatives() + 1;
}

real RotamerLibrary::rotamerProbability(string aa, int ri, real phi, real psi, bool strict) {
  int bi = getBackboneBin(aa, phi, psi, !strict);
  return prob[aa][bi][ri];
}

real RotamerLibrary::rotamerProbability(rotamerID* rot) {
  return prob[rot->aminoAcid()][rot->binIndex()][rot->rotIndex()];
}

rotamerID RotamerLibrary::getRotamer(Residue& res, string aa, int rotIndex, bool strict) {
  double phi = res.getPhi();
  double psi = res.getPsi();
  int binInd = getBackboneBin(aa, phi, psi, !strict);
  if (rotamers.find(aa) == rotamers.end()) MstUtils::error("rotamer library does not contain amino acid '" + aa + "'", "RotamerLibrary::placeRotamer");
  if (rotamers[aa].size() <= binInd) {
    MstUtils::error("rotamer library has " + MstUtils::toString(rotamers[aa].size()) + " backbone bins for amino acid '" + aa + "', but bin number " + MstUtils::toString(binInd+1) + " was requested", "RotamerLibrary::placeRotamer");
  }
  return rotamerID(aa, binInd, rotIndex);
}

rotamerID RotamerLibrary::placeRotamer(Residue& res, string aa, int rotIndex, Residue* dest_ptr, bool strict) {
  double phi = res.getPhi();
  double psi = res.getPsi();
  int binInd = getBackboneBin(aa, phi, psi, !strict);
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
  CartesianPoint X = CA - N;          // X-axis of the residue frame defined to be along N -> CA
  CartesianPoint Z = X.cross(C - CA); // since CA-C is defined to be in the XY plane, (N -> CA) x (CA -> C) will be a vector along Z
  CartesianPoint Y = Z.cross(X);      // finally, Z x X is Y
  Frame R(CA, X, Y, Z);               // residue frame
  Frame L;                            // Frame class defaults to the laboratory frame
  Transform T = TransformFactory::switchFrames(R, L);

  // fish out the right rotamer and transform it onto the residue
  vector<Atom*> newAtoms; transformRotamerAtoms(T, rots, rotIndex, newAtoms);
  if (dest_ptr == NULL) {
    vector<int> oldAtoms; // atoms to be deleted (mostly side-chain atoms of the previous residue)
    for (int i = 0; i < res.atomSize(); i++) {
      Atom& a = res[i];
      if (!isBackboneAtom(a.getName()) || ((aa.compare("PRO") == 0) && a.isNamed("H"))) {
        oldAtoms.push_back(i);
      }
    }

    // copy rotamer atoms over to the residue, destroying previous atoms
    res.replaceAtoms(newAtoms, &oldAtoms);
    res.setName(aa);
  } else {
    Residue& dest = *dest_ptr;
    // this means a previous residue has been specified. It should be either empty
    // or have been previously filled by the same function.
    if (dest.atomSize() == 0) {
      dest.setName(aa);
      dest.appendAtoms(newAtoms); // first the side-chain
      for (int i = 0; i < res.atomSize(); i++) {
        Atom& a = res[i];
        if (isBackboneAtom(a.getName()) && !((aa.compare("PRO") == 0) && a.isNamed("H"))) {
          dest.appendAtom(new Atom(a));
        }
      }
    } else {
      if (!dest.isNamed(aa)) MstUtils::error("specified destination residue is neither empty nor does it contain a rotamer of the same amino acid", "RotamerLibrary::placeRotamer");
      for (int i = 0; i < newAtoms.size(); i++) {
        if (!dest[i].isNamed(newAtoms[i]->getNameC())) MstUtils::error("specified destination residue is neither empty nor does it contain a rotamer of the same amino acid", "RotamerLibrary::placeRotamer");
        dest[i].setCoor(newAtoms[i]->getCoor());
        delete newAtoms[i];
      }
    }
  }
  return rotamerID(aa, binInd, rotIndex);
}

void RotamerLibrary::transformRotamerAtoms(Transform& T, Residue& rots, int rotIndex, vector<Atom*>& newAtoms) {
  newAtoms.resize(rots.atomSize(), NULL);
  for (int i = 0; i < rots.atomSize(); i++) {
    Atom& a = rots[i];
    if (a.numAlternatives() < rotIndex) {
      MstUtils::error("rotamer library contains " + MstUtils::toString(a.numAlternatives() + 1) + " rotamers for amino-acid, but rotamer number " + MstUtils::toString(rotIndex+1) + "was requested", "RotamerLibrary::placeRotamer");
    }
    if (rotIndex > 0) a.swapWithAlternative(rotIndex-1);
    Atom* newAtom = new Atom(a, false);
    T.apply(newAtom);
    newAtoms[i] = newAtom;
    if (rotIndex > 0) a.swapWithAlternative(rotIndex-1);
  }
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

bool RotamerLibrary::isHydrogen(string atomName) {
  if (atomName.length() == 0) return false;
  if (atomName.compare(0, 1, "H") == 0) return true;
  return false;
}

int RotamerLibrary::getBackboneBin(string aa, real phi, real psi, bool assumeDefault) {
  if (binPhiCenters.find(aa) == binPhiCenters.end()) MstUtils::error("no PHI bin for amino acid '" + aa + "'", "RotamerLibrary::getBackboneBin");
  if (binPsiCenters.find(aa) == binPsiCenters.end()) MstUtils::error("no PSI bin for amino acid '" + aa + "'", "RotamerLibrary::getBackboneBin");

  if ((phi == Residue::badDihedral) || (psi == Residue::badDihedral)) {
    if (assumeDefault) {
      // if no valid phi/psi pair provided, find the most frequent bin for the given amino acid
      return defaultBin[aa];
    } else {
      MstUtils::error("could not find backbone bin for amino acid '" + aa + "', because provided phi/psi pair is invalid: " + MstUtils::toString(phi) + "/" + MstUtils::toString(psi), "RotamerLibrary::getBackboneBin");
    }
  }

  int phiInd = findClosestAngle(binPhiCenters[aa], phi);
  int psiInd = findClosestAngle(binPsiCenters[aa], psi);

  // by convention, linear bin index is encoded in phi-major order
  return phiInd * binPsiCenters[aa].size() + psiInd;
}

pair<real, real> RotamerLibrary::getBinPhiPsi(string aa, int bi) {
  if (binPhiCenters.find(aa) == binPhiCenters.end()) MstUtils::error("no PHI bin for amino acid '" + aa + "'", "RotamerLibrary::getBackboneBin");
  if (binPsiCenters.find(aa) == binPsiCenters.end()) MstUtils::error("no PSI bin for amino acid '" + aa + "'", "RotamerLibrary::getBackboneBin");
  if (bi >= binPhiCenters[aa].size() * binPsiCenters[aa].size()) MstUtils::error("bin index " + MstUtils::toString(bi) + " out of range", "RotamerLibrary::getBackboneBin");

  int phiInd = bi / binPsiCenters[aa].size();
  int psiInd = bi % binPsiCenters[aa].size();
  return pair<real, real>(binPhiCenters[aa][phiInd], binPsiCenters[aa][psiInd]);
}

int RotamerLibrary::findClosestAngle(vector<real>& array, real value) {
  MstUtils::assert(array.size() != 0, "empty vector passed", "RotamerLibrary::findClosestAngle");
  int minInd = 0;
  real minDist = fabs(angleDiff(value, array[0]));
  real dist;
  for (int i = 1; i < array.size(); i++) {
    dist = fabs(angleDiff(value, array[i]));
    if (minDist > dist) {
      minDist = dist;
      minInd = i;
    }
  }
  return minInd;
}

real RotamerLibrary::angleDiff(real a, real b) {

  real da = MstUtils::mod((MstUtils::mod(a, 360.0) - MstUtils::mod(b, 360.0)), 360.0);
  if (da > 180.0) da -= 360.0;

  return da;
}
