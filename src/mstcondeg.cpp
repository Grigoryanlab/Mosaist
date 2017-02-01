#include "mstcondeg.h"

aaProp["ALA"] = 7.73; aaProp["CYS"] = 1.84; aaProp["ASP"] = 5.82; aaProp["GLU"] = 6.61; aaProp["PHE"] = 4.05;
aaProp["GLY"] = 7.11; aaProp["HIS"] = 2.35; aaProp["HSD"] = 2.35; aaProp["ILE"] = 5.66; aaProp["LYS"] = 6.27;
aaProp["LEU"] = 8.83; aaProp["MET"] = 2.08; aaProp["ASN"] = 4.50; aaProp["PRO"] = 4.52; aaProp["GLN"] = 3.94;
aaProp["ARG"] = 5.03; aaProp["SER"] = 6.13; aaProp["THR"] = 5.53; aaProp["VAL"] = 6.91; aaProp["TRP"] = 1.51; aaProp["TYR"] = 3.54;

ConFind::ConFind(string rotLibFile, Structure& _S) {
  rotLib = new RotamerLibrary(rotLibFile);
  init(_S);
}
ConFind::ConFind(RotamerLibrary* _rotLib, Structure& S) {
  rotLib = _rotLib;
  init(_S);
}

ConFind::~ConFind() {
  delete rotLib;
}

ConFind::init(Structure& _S) {
  // TODO: initilialize these data structures:
  AtomPointerVector backbone, nonHyd, ca;
  ProximitySearch bbNN, nonHydNN, caNN;
  ProximitySearch rotamerHeavySC, rotamerHeavyBB; // point cloud of rotamer atoms from ALL rotamers
}

bool ConFind::cache(Residue* residue) {
  if (fullyDescribed.find(residue) != fullyDescribed.end()) return true;
  bool doNotCountCB = true; // if true, CB is not counted as a side-chain atom for counting clashes (except for ALA)

  // find all residues around this one that are within cutoff distance and can affects it
  vector<int> close = caNN.getPointsWithin(residue.findAtom("CA"), 0, dcut);
  vector<Residue*> residues(surr.size(), NULL);
  bool foundSelf = false
  for (int i = 0; i < close.size(); i++) {
    residues[i] = ca[surr[i]]->getResidue();
    if (residues[i] == residue) foundSelf = true;
  }
  if (!foundSelf)
    MstUtils::error("when looking for residues around " + MstUtils::toString(*residue) + " did not find self!", "ConFind::cache(Residue*)");

  // rotamers are cached or everybody, while the rest of the stuff only for the residue in question
  for (int i = 0; i < residues.size(); i++) {
    Residue &res = *(residues[i]);
    bool ofInterest = (&res == residue); // is this a position of interest or a supporting position?
    // has this has already been visited (as either support for another residue or as a residue of interest)?
    if (rotamers.find(residues[i]) != rotamers.end()) continue;
    double phi = res.getPhi(); double psi = res.getPsi();

    // make sure this residue has a proper backbone, otherwise adding rotamers will fail
    if (!res.atomExists("N") || !res.atomExists("CA") || !res.atomExists("C")) {
      if (ofInterest) {
        MstUtils::warning("residue lacks proper backbone for placing rotamers!", "ConFind::cache(Residue*)");
        return false;
      }
      MstUtils::warning("will not build rotamers at position " + MstUtils::toString(*res) + " as it lacks proper backbone. " +
                        "NOTE: this will affect the correctness of some values, e.g., for residue " << MstUtils::toString(*residue));
      continue;
    }

    // load rotamers of each amino acid
    int numRemRotsInPosition = 0; int totNumRotsInPosition = 0;
    double freeVolumeNum = 0; double freeVolumeDen = 0;
    for (int j = 0; j < aaNames.size(); j++) {
      string aa = aaNames[j];
      if (aaProp.find(aa) == aaProp.end()) MstUtils::error("no propensity defined for amino acid " + aa);
      double aaProp = aaProp[aa];
      int nr = rotLib.numberOfRotamers(aaNames[i], phi, psi);
      Residue resCopy;
      resCopy.copyAtoms(res); // make a copy of the original residue for the purposes of placing rotamers, so as not to destroy the structure
      rotamers[res].push_back(new aaRotamers(aa));
      for (int ri = 0; ri < nr; ri++) {
        rotLib.placeRotamer(resCopy, aaNames[i], ri);
        double rotProb = rotLib.rotamerProbability(aaNames[i], ri, phi, psi);

        // see if the rotamer needs to be pruned (clash with the backbone).
        // also measure contribution to the free volume
        bool prune = false;
        vector<int> closeOnes;
        for (int k = 0; k < resCopy.atomSize(); k++) {
          Atom& a = resCopy[k];
          if (RotamerLibrary::isBackboneAtom(a)) continue;
          // free volume contributions
          closeOnes.clear();
          nonHydNN.pointsWithin(a.getCoor(), 0.0, _opt.contDist, closeOnes);
          for (int ci = 0; ci < closeOnes.size(); ci++) {
            // self clashes do not count (the rotamer library should not allow true clashes with own backbone or own side-chain)
            if (nonHyd[closeOnes[ci]]->getResidue() != res) {
              freeVolumeNum += aaProp*rotProb;
              break;
            }
          }
          freeVolumeDen += aaProp*rotProb;

          // shuld the rotamer be pruned?
          if (doNotCountCB && !resCopy.isNamed("ALA") && a.isNamed("CB")) continue;
          closeOnes.clear();
          bbNN.pointsWithin(a.getCoor(), 0.0, _opt.clashDist, closeOnes);
          for (int ci = 0; ci < closeOnes.size(); ci++) {
            // backbone atoms of the same residue do not count as clashing (the rotamer library should not allow true clashes with own backbone)
            if (backbone[closeOnes[ci]]->getResidue() != res) {
              prune = true;
              // clashes with ALA have a special meaning (permanent "unavoidable" contacts; need to find all of them, though unlikely to have more than one)
              if (resCopy.isNamed("ALA")) permanentContacts[residue].insert(closeOnes[ci]);
              else break;
            }
          }
          if (prune) break;
        }
        if (prune) continue;

        // if not pruned, cache info about it
        rotamers[res].back()->addRotamer(resCopy);
        for (int ai = 0; ai < resCopy.atomSize(); ai++) {
          if (!RotamerLibrary::isHydrogen(resCopy[ai])) {
            if (RotamerLibrary::isBackboneAtom(resCopy[ai]) || (doNotCountCB && !resCopy.isNamed("ALA") && resCopy[ai].isNamed("CB"))) {
              rotamerHeavyBB.addPoint(resCopy[ai]);
            } else {
              rotamerHeavySC.addPoint(resCopy[ai]);
            }
          }
        }
        numRemRots++;
        numRemRotsInPosition++;
      }
      totNumRotsInPosition += nr;
    }
    fractionPruned[residue] = (totNumRotsInPosition - numRemRotsInPosition)*1.0/totNumRotsInPosition;
    origNumRots[residue] = totNumRotsInPosition;
    freeVolume[residue] = 1 - freeVolumeNum/freeVolumeDen;
  }
  fullyDescribed[residue] = true;
  return true;
}

bool ConFind::cache(vector<Residue*>& residues) {
  bool ret = true;
  for (int i = 0; i < residues.size(); i++) {
    ret = ret && cache(residues[i]);
  }
  return ret;
}

bool ConFind::cache(Structure& S) {
  vector<Residue*> residues = S.getResidues();
  return cache(residues);
}
