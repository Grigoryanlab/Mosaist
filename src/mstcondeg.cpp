#include "mstcondeg.h"

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
}

void ConFind::cache(Residue* res) {
  bool doNotCountCB = true; // if true, CB is not counted as a side-chain atom for counting clashes (except for ALA)

  // find all residues around this one that are within cutoff distance and can affects it
  vector<int> surr = caNN.getPointsWithin(res.findAtom("CA"), 0, dcut);
  vector<Residue*> surrResidues(surr.size(), NULL);
  for (int i = 0; i < surr.size(); i++) surrResidues[i] = ca[surr[i]]->getResidue();

  // place rotamers onto residue res, and cache that info (NN rotamer object)

  // TODO: I have not yet expanded the selection to include everything with 25 A of the specified residues!!!!!!!!!!
  for (int i = 0; i < residues.size(); i++) {
    Residue &res = *(residues[i]);
    if (_opt.verbose) cout << "position " << res << ", " << i+1 "/" << (int) residues.size() << endl;
    double phi = res.getPhi(); double psi = res.getPsi();
    pp[i][0] = phi;
    pp[i][1] = psi;
    pp[i][2] = res.getOmega(); // ???? CHECK THE MSL NOMENCLATURE ABOUT WHICH RESIDUE THE OMEGA BELONGS TO!

    // make sure this residue has a proper backbone, otherwise adding rotamers will fail
    if (!res.atomExists("N") || !res.atomExists("CA") || !res.atomExists("C")) {
      cout << "\nWarning: will not build rotamers at position " << res << ", because N-CA-C atoms are not defined\n";
      cout << "NOTE: this can affect the correctness of contact degree and crowdedness calculations at neighboring positions.\n" << endl;
      for (int ii = 0; ii < res.atomSize(); ii++) cout << res[ii] << endl;
      continue;
    }

    // load rotamers of each amino acid
    int numRemRotsInPosition = 0; int totNumRotsInPosition = 0;
    double freeVolumeNum = 0; double freeVolumeDen = 0;
    for (int j = 0; j < aaNames.size(); j++) {
      if (_opt.aaProp.find(aaNames[j]) == _opt.aaProp.end()) MstUtils::error("no propensity defined for amino acid " + aaNames[j]);
      double aaProp = _opt.aaProp[aaNames[j]];
      if (_opt.verbose) printf("%s %.3f: ", aaNames[j].c_str(), aaProp);

      int nr = rotLib.numberOfRotamers(aaNames[i], phi, psi);
      Residue resCopy;
      resCopy.copyAtoms(res); // make a copy of the original residue for the purposes of placing rotamers, so as not to destroy the structure
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
              if (resCopy.isNamed("ALA")) permanentContacts[i].insert(closeOnes[ci]);
              else break;
            }
          }
          if (prune) break;
        }
        if (prune) continue;

        if (!_opt.rotOutFile.empty()) {
          rof << "REM " << resCopy << ", rotamer " << r+1 << endl;
          Structure rot;
          rot.addAtoms(resCopy.getAtoms());
          rot.writePDB(rof);
        }
        AtomPointerVector heavySC, heavyBB;
        for (int ai = 0; ai < resCopy.atomSize(); ai++) {
          if (!RotamerLibrary::isHydrogen(resCopy[ai])) {
            if (RotamerLibrary::isBackboneAtom(resCopy[ai]) || (doNotCountCB && !resCopy.isNamed("ALA") && resCopy[ai].isNamed("CB"))) {
              heavyBB.push_back(&(resCopy[ai]));
            } else {
              heavySC.push_back(&(resCopy[ai]));
            }
          }
        }
        rotamers[i].push_back(rotamer(new ProximitySearch(heavySC, _opt.contDist), new ProximitySearch(heavyBB, _opt.contDist), aaProp, rotProb, aaNames[j], r));
        if (_opt.verbose) printf("%d %.3f; ", r, rotProb);
        numRemRots++;
        numRemRotsInPosition++;
      }
      totNumRotsInPosition += nr;
    }
    if (_opt.verbose) cout << numRemRots << "/" << totNumRotsInPosition << " remaining at position " << res << endl;
    fractionPruned[i] = (totNumRotsInPosition - numRemRotsInPosition)*1.0/totNumRotsInPosition;
    origNumRots[i] = totNumRotsInPosition;
    freeVolume[i] = 1 - freeVolumeNum/freeVolumeDen;
  }
  if (_opt.verbose) cout << "end of rotamer filtering..." << endl;
  if (!_opt.rotOutFile.empty()) rof.close();
}
