#include "mstcondeg.h"

ConFind::ConFind(string rotLibFile, Structure& S) {
  setParams();
  rotLib = new RotamerLibrary(rotLibFile);
  isRotLibLocal = true;
  init(S);
}
ConFind::ConFind(RotamerLibrary* _rotLib, Structure& S) {
  setParams();
  rotLib = _rotLib;
  isRotLibLocal = false;
  init(S);
}

void ConFind::setParams() {
  dcut = 25.0;
  clashDist = 2.0;
  contDist = 3.0;
  aaNames = {"ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "ILE", "LEU",
             "LYS", "MET", "PHE", "SER", "THR", "TRP", "TYR", "VAL", "ALA"}; // all but GLY and PRO
  aaProp["ALA"] = 7.73; aaProp["CYS"] = 1.84; aaProp["ASP"] = 5.82; aaProp["GLU"] = 6.61; aaProp["PHE"] = 4.05;
  aaProp["GLY"] = 7.11; aaProp["HIS"] = 2.35; aaProp["HSD"] = 2.35; aaProp["ILE"] = 5.66; aaProp["LYS"] = 6.27;
  aaProp["LEU"] = 8.83; aaProp["MET"] = 2.08; aaProp["ASN"] = 4.50; aaProp["PRO"] = 4.52; aaProp["GLN"] = 3.94;
  aaProp["ARG"] = 5.03; aaProp["SER"] = 6.13; aaProp["THR"] = 5.53; aaProp["VAL"] = 6.91; aaProp["TRP"] = 1.51; aaProp["TYR"] = 3.54;
}

ConFind::~ConFind() {
  if (isRotLibLocal) delete rotLib;
  for (map<Residue*, vector<aaRotamers*> >::iterator it = rotamers.begin(); it != rotamers.end(); ++it) {
    vector<aaRotamers*>& aaRots = it->second;
    for (int i = 0; i < aaRots.size(); i++) delete(aaRots[i]);
  }
  delete bbNN;
  delete nonHydNN;
  delete caNN;
  delete rotamerHeavySC;
}

void ConFind::init(Structure& S) {
  AtomPointerVector atoms = S.getAtoms();
  for (int i = 0; i < atoms.size(); i++) {
    Atom* a = atoms[i];
    if (a->isNamed("H")) continue;
    if (RotamerLibrary::isBackboneAtom(a)) backbone.push_back(a);
    nonHyd.push_back(a);
    if (a->isNamed("CA")) ca.push_back(a);
  }
  bbNN = new ProximitySearch(backbone, clashDist/2);
  nonHydNN = new ProximitySearch(nonHyd, clashDist/2);
  caNN = new ProximitySearch(ca, dcut/2);

  // add a padding of 50 A to the structure to cover all future rotamer atoms
  rotamerHeavySC = new DecoratedProximitySearch<rotamerAtomInfo>(atoms, contDist/2, 50.0);
}

bool ConFind::cache(Residue* res) {
cout << "working on " << MstUtils::toString(*res) << endl;
  if (rotamers.find(res) != rotamers.end()) return true;
  bool doNotCountCB = true; // if true, CB is not counted as a side-chain atom for counting clashes (except for ALA)

  // make sure this residue has a proper backbone, otherwise adding rotamers will fail
  if (!res->atomExists("N") || !res->atomExists("CA") || !res->atomExists("C")) {
    MstUtils::warn("cannot build rotamers at position " + MstUtils::toString(*res) + " as it lacks proper backbone!", "ConFind::cache(Residue*)");
    return false;
  }
cout << "trying to get phi/psi for " << MstUtils::toString(*res) << endl;
  real phi = res->getPhi(); real psi = res->getPsi();
cout << "phi = " << phi << ", psi = " << psi << endl;

  // load rotamers of each amino acid
  int numRemRotsInPosition = 0; int totNumRotsInPosition = 0;
  double freeVolumeNum = 0; double freeVolumeDen = 0;
  for (int j = 0; j < aaNames.size(); j++) {
    string aa = aaNames[j];
    if (aaProp.find(aa) == aaProp.end()) MstUtils::error("no propensity defined for amino acid " + aa);
    double aaP = aaProp[aa];
    int nr = rotLib->numberOfRotamers(aa, phi, psi);
    Residue resCopy;
    resCopy.copyAtoms(*res); // make a copy of the original residue for the purposes of placing rotamers, so as not to destroy the structure
    rotamers[res].push_back(new aaRotamers(aa, res));
    aaRotamers& aaRots = *(rotamers[res].back());
    for (int ri = 0; ri < nr; ri++) {
      rotLib->placeRotamer(resCopy, aa, ri);
      double rotP = rotLib->rotamerProbability(aa, ri, phi, psi);

      // see if the rotamer needs to be pruned (clash with the backbone).
      // also measure contribution to the free volume
      bool prune = false;
      vector<int> closeOnes;
      for (int k = 0; k < resCopy.atomSize(); k++) {
        Atom& a = resCopy[k];
        if (RotamerLibrary::isBackboneAtom(a)) continue;
        // free volume contributions
        closeOnes.clear();
        nonHydNN->pointsWithin(a.getCoor(), 0.0, contDist, &closeOnes);
        for (int ci = 0; ci < closeOnes.size(); ci++) {
          // self clashes do not count (the rotamer library should not allow true clashes with own backbone or own side-chain)
          if (nonHyd[closeOnes[ci]]->getResidue() != res) {
            freeVolumeNum += aaP*rotP;
            break;
          }
        }
        freeVolumeDen += aaP*rotP;

        // shuld the rotamer be pruned?
        if (doNotCountCB && !resCopy.isNamed("ALA") && a.isNamed("CB")) continue;
        closeOnes.clear();
        bbNN->pointsWithin(a.getCoor(), 0.0, clashDist, &closeOnes);
        for (int ci = 0; ci < closeOnes.size(); ci++) {
          // backbone atoms of the same residue do not count as clashing (the rotamer library should not allow true clashes with own backbone)
          if (backbone[closeOnes[ci]]->getResidue() != res) {
            prune = true;
            // clashes with ALA have a special meaning (permanent "unavoidable" contacts; need to find all of them, though unlikely to have more than one)
            if (resCopy.isNamed("ALA")) permanentContacts[res].insert(closeOnes[ci]);
            else break;
          }
        }
        if (prune) break;
      }
      if (prune) continue;

      // if not pruned, cache info about it
      aaRots.addRotamer(resCopy, ri, rotP);
      for (int ai = 0; ai < resCopy.atomSize(); ai++) {
        if (!RotamerLibrary::isHydrogen(resCopy[ai])) {
          Atom& a = resCopy[ai];
          if (RotamerLibrary::isBackboneAtom(a)) continue;
          if (doNotCountCB && a.isNamed("CB") && !resCopy.isNamed("ALA")) continue;
          rotamerAtomInfo rotAtom(&aaRots, aaRots.numberOfRotamers(), ai);
          rotamerHeavySC->addPoint(a, rotAtom);
        }
      }
      numRemRotsInPosition++;
    }
    totNumRotsInPosition += nr;
  }
  fractionPruned[res] = (totNumRotsInPosition - numRemRotsInPosition)*1.0/totNumRotsInPosition;
  origNumRots[res] = totNumRotsInPosition;
  freeVolume[res] = 1 - freeVolumeNum/freeVolumeDen;

  return true;
}

bool ConFind::cache(vector<Residue*>& residues) {
cout << "and then here" << endl;
  bool ret = true;
  for (int i = 0; i < residues.size(); i++) {
    ret = ret && cache(residues[i]);
  }
  return ret;
}

bool ConFind::cache(Structure& S) {
cout << "and here" << endl;
  vector<Residue*> residues = S.getResidues();
  return cache(residues);
}

contactList ConFind::getContacts(Residue* res, real cdcut) {
  vector<Residue*> neighborhood = getNeighbors(res);
  cache(neighborhood);
  map<Residue*, real> conDegNum;

  // go over all the rotamers at the current residue and see what they contact
  vector<aaRotamers*>& posRots = rotamers[res];
  for (int aai = 0; aai < posRots.size(); aai++) {
    aaRotamers& aaRots = *(posRots[aai]);
    string aani = aaRots.aaName();
    for (int ri = 0; ri < aaRots.numberOfRotamers(); ri++) {
      real pi = aaRots.rotProb(ri);
      map<Residue*, map<int, bool> > conRots; // collection of rotamers at other positions that this rotamer contacts
      for (int ai = 0; ai < aaRots.atomSize(); ai++) {
        vector<rotamerAtomInfo> conts = rotamerHeavySC->getPointsWithin(aaRots.rotamerAtom(ri, ai), 0, contDist);
        for (int ci = 0; ci < conts.size(); ci++) {
          Residue* cres = conts[ci].position();
          if (cres == res) continue; // clashes with rotamers at the same position do not count
          int rotInd = conts[ci].rotIndex();
          if (conRots.find(cres) == conRots.end()) conRots[cres] = map<int, bool>();
          if (conRots[cres].find(rotInd) == conRots[cres].end()) {
            // a new contact detected, add its contribution
            string aanj = conts[ci].aaName();
            real pj = conts[ci].rotProb();
            if (conDegNum.find(cres) == conDegNum.end()) conDegNum[cres] = 0.0;
            conDegNum[cres] += aaProp[aani] * aaProp[aanj] * pi * pj;
          }
          conRots[cres][rotInd] = true;
        }
      }
    }
  }

  // go over all contacting residues and finish contact degree calculation
  contactList L(res);
  real wi = weightOfAvailableRotamers(res);
  for (map<Residue*, real>::iterator it = conDegNum.begin(); it != conDegNum.end(); ++it) {
    Residue* cres = it->first;
    real cd = it->second / (wi * weightOfAvailableRotamers(cres));
    if (cd > cdcut) L.addContact(cres, cd);
  }
  return L;
}

real ConFind::weightOfAvailableRotamers(Residue* res) {
  real weightAA = 0; real weightRot = 0;
  if (rotamers.find(res) == rotamers.end()) MstUtils::error("residue not cached: " + MstUtils::toString(*res), "ConFind::weightOfAvailableRotamers");
  for (int i = 0; i < rotamers[res].size(); i++) {
    aaRotamers& aaRots = *(rotamers[res][i]);
    for (int ri = 0; ri < aaRots.numberOfRotamers(); ri++) {
      weightRot += aaRots.rotProb(ri);
      weightAA += aaProp[aaRots.aaName()];
    }
  }
  return weightAA * weightRot;
}

vector<Residue*> ConFind::getNeighbors(Residue* residue) {
  // find all residues around this one that are within cutoff distance and can affects it
  vector<int> close = caNN->getPointsWithin(residue->findAtom("CA"), 0, dcut);
  vector<Residue*> neighborhood(close.size(), NULL);
  bool foundSelf = false;
  for (int i = 0; i < close.size(); i++) {
    neighborhood[i] = ca[close[i]]->getResidue();
    if (neighborhood[i] == residue) foundSelf = true;
  }
  if (!foundSelf) {
    MstUtils::error("when looking for residues around " + MstUtils::toString(*residue) + " did not find self!", "ConFind::cache(Residue*)");
  }
  return neighborhood;
}

vector<Residue*> ConFind::getNeighbors(vector<Residue*>& residues) {
  // find all residues that are within cutoff distance of the given list of residues
  map<Residue*, bool> within;
  for (int k = 0; k < residues.size(); k++) {
    vector<int> close = caNN->getPointsWithin(residues[k]->findAtom("CA"), 0, dcut);
    bool foundSelf = false;
    for (int i = 0; i < close.size(); i++) {
      Residue* wres = ca[close[i]]->getResidue();
      within[wres] = true;
      if (residues[i] == wres) foundSelf = true;
    }
    if (!foundSelf) {
      MstUtils::error("when looking for residues around " + MstUtils::toString(*(residues[k])) + " did not find self!", "ConFind::cache(Residue*)");
    }
  }
  vector<Residue*> neighborhood(within.size(), NULL);
  map<Residue*, bool>::iterator it; int i;
  for (it = within.begin(), i = 0; it != within.end(); it++, i++) {
    neighborhood[i] = it->first;
  }
  return neighborhood;
}
