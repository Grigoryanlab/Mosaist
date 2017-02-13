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
  doNotCountCB = true;
  aaNames = {"ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "ILE", "LEU",
             "LYS", "MET", "PHE", "SER", "THR", "TRP", "TYR", "VAL", "ALA"}; // all but GLY and PRO
  aaProp["ALA"] = 7.73; aaProp["CYS"] = 1.84; aaProp["ASP"] = 5.82; aaProp["GLU"] = 6.61; aaProp["PHE"] = 4.05;
  aaProp["GLY"] = 7.11; aaProp["HIS"] = 2.35; aaProp["HSD"] = 2.35; aaProp["ILE"] = 5.66; aaProp["LYS"] = 6.27;
  aaProp["LEU"] = 8.83; aaProp["MET"] = 2.08; aaProp["ASN"] = 4.50; aaProp["PRO"] = 4.52; aaProp["GLN"] = 3.94;
  aaProp["ARG"] = 5.03; aaProp["SER"] = 6.13; aaProp["THR"] = 5.53; aaProp["VAL"] = 6.91; aaProp["TRP"] = 1.51; aaProp["TYR"] = 3.54;
}

ConFind::~ConFind() {
  if (isRotLibLocal) delete rotLib;
  delete bbNN;
  delete nonHydNN;
  delete caNN;
  for (map<Residue*, vector<rotamerID*> >::iterator it = survivingRotamers.begin(); it != survivingRotamers.end(); ++it) {
    vector<rotamerID*>& rots = it->second;
    for (int i = 0; i < rots.size(); i++) delete(rots[i]);
  }
  for (map<Residue*, DecoratedProximitySearch<rotamerID*>* >::iterator it = rotamerHeavySC.begin(); it != rotamerHeavySC.end(); ++it) {
    delete it->second;
  }
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
}

bool ConFind::cache(Residue* res, fstream* rotOut) {
  if (rotamerHeavySC.find(res) != rotamerHeavySC.end()) return true;
  AtomPointerVector pointCloud;      // side-chain atoms of surviving rotames
  vector<rotamerID*> pointCloudTags; // corresponding tags (i.e.,  rotamer identity)

  // make sure this residue has a proper backbone, otherwise adding rotamers will fail
  if (!res->atomExists("N") || !res->atomExists("CA") || !res->atomExists("C")) {
    MstUtils::warn("cannot build rotamers at position " + MstUtils::toString(*res) + " as it lacks proper backbone!", "ConFind::cache(Residue*)");
    return false;
  }
  real phi = res->getPhi(); real psi = res->getPsi();

  // load rotamers of each amino acid
  int numRemRotsInPosition = 0; int totNumRotsInPosition = 0;
  double freeVolumeNum = 0; double freeVolumeDen = 0;
  for (int j = 0; j < aaNames.size(); j++) {
    string aa = aaNames[j];
    if (aaProp.find(aa) == aaProp.end()) MstUtils::error("no propensity defined for amino acid " + aa);
    double aaP = aaProp[aa];
    int nr = rotLib->numberOfRotamers(aa, phi, psi);
    Residue rot;
    for (int ri = 0; ri < nr; ri++) {
      rotamerID rID = rotLib->placeRotamer(*res, aa, ri, &rot);
      double rotP = rotLib->rotamerProbability(aa, ri, phi, psi);

      // see if the rotamer needs to be pruned (clash with the backbone).
      // also measure contribution to the free volume
      bool prune = false;
      vector<int> closeOnes;
      for (int k = 0; k < rot.atomSize(); k++) {
        Atom& a = rot[k];
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
        if (!countsAsSidechain(a)) continue;
        closeOnes.clear();
        bbNN->pointsWithin(a.getCoor(), 0.0, clashDist, &closeOnes);
        for (int ci = 0; ci < closeOnes.size(); ci++) {
          // backbone atoms of the same residue do not count as clashing (the rotamer library should not allow true clashes with own backbone)
          if (backbone[closeOnes[ci]]->getResidue() != res) {
            prune = true;
            // clashes with ALA have a special meaning (permanent "unavoidable" contacts; need to find all of them, though unlikely to have more than one)
            if (rot.isNamed("ALA")) permanentContacts[res].insert(closeOnes[ci]);
            else break;
          }
        }
        if (prune) break;
      }
      if (prune) continue;
      if (rotOut != NULL) { Structure S(rot); S.writePDB(*rotOut); }

      // if not pruned, collect atoms needed later
      rotamerID* rotTag = new rotamerID(rID);
      survivingRotamers[res].push_back(rotTag);
      for (int ai = 0; ai < rot.atomSize(); ai++) {
        if (!countsAsSidechain(rot[ai])) continue;
        pointCloud.push_back(new Atom(rot[ai]));
        pointCloudTags.push_back(rotTag);
      }
      numRemRotsInPosition++;
    }
    totNumRotsInPosition += nr;
  }
  fractionPruned[res] = (totNumRotsInPosition - numRemRotsInPosition)*1.0/totNumRotsInPosition;
  origNumRots[res] = totNumRotsInPosition;
  freeVolume[res] = 1 - freeVolumeNum/freeVolumeDen;

  // cash all atoms for faster distance-based searches
  rotamerHeavySC[res] = new DecoratedProximitySearch<rotamerID*>(pointCloud, contDist/2, pointCloudTags);
  pointCloud.deletePointers();

  return true;
}

bool ConFind::countsAsSidechain(Atom& a) {
  if (RotamerLibrary::isHydrogen(a) || RotamerLibrary::isBackboneAtom(a)) return false;
  if (doNotCountCB && a.isNamed("CB") && !(a.getResidue()->isNamed("ALA"))) return false;
  return true;
}
bool ConFind::cache(vector<Residue*>& residues, fstream* rotOut) {
  bool ret = true;
  for (int i = 0; i < residues.size(); i++) {
    ret = ret && cache(residues[i], rotOut);
  }
  return ret;
}

bool ConFind::cache(Structure& S, fstream* rotOut) {
  vector<Residue*> residues = S.getResidues();
  return cache(residues, rotOut);
}

real ConFind::contactDegree(Residue* resA, Residue* resB, bool doNotCache) {
  if (!doNotCache) { cache(resA); cache(resB); }
  DecoratedProximitySearch<rotamerID*>& cloudA = *(rotamerHeavySC[resA]);
  DecoratedProximitySearch<rotamerID*>& cloudB = *(rotamerHeavySC[resB]);

  /* find rotamer pairs that clash (could account for how many times and
   * with which atoms, if we wanted) */
  map<rotamerID*, map<rotamerID*, bool> > clashing;
  for (int ai = 0; ai < cloudA.pointSize(); ai++) {
    vector<rotamerID*> p = cloudB.getPointsWithin(cloudA.getPoint(ai), 0, contDist);
    if (p.size() == 0) continue;
    rotamerID* rID = cloudA.getPointTag(ai);
    for (int i = 0; i < p.size(); i++) clashing[rID][p[i]] = true;
  }

  // compute contact degree
  real cd = 0;
  for (map<rotamerID*, map<rotamerID*, bool> >::iterator itA = clashing.begin(); itA != clashing.end(); ++itA) {
    rotamerID& rotA = *(itA->first);
    real rotProbA = rotLib->rotamerProbability(rotA);
    real aaPropA = aaProp[rotA.aminoAcid()];
    for (map<rotamerID*, bool>::iterator itB = itA->second.begin(); itB != itA->second.end(); ++itB) {
      rotamerID& rotB = *(itB->first);
      real rotProbB = rotLib->rotamerProbability(rotB);
      real aaPropB = aaProp[rotB.aminoAcid()];
      cd +=  aaPropA * aaPropB * rotProbA * rotProbB;
    }
  }
  cd /= weightOfAvailableRotamers(resA) * weightOfAvailableRotamers(resB);
  return cd;
}

contactList ConFind::getContacts(Residue* res, real cdcut) {
  vector<Residue*> neighborhood = getNeighbors(res);
  cache(neighborhood);

  // compute contact degree between this residue and every one of its neighbors
  contactList L;
  for (int i = 0; i < neighborhood.size(); i++) {
    if (res == neighborhood[i]) continue;
    real cd = contactDegree(res, neighborhood[i], true);
    if (cd > cdcut) L.addContact(res, neighborhood[i], cd);
  }
  return L;
}

contactList ConFind::getContacts(Structure& S, real cdcut) {
  cache(S);
  vector<Residue*> allRes = S.getResidues();
  map<Residue*, map<Residue*, bool> > checked;

  contactList L;
  for (int i = 0; i < allRes.size(); i++) {
    Residue* resi = allRes[i];
    vector<Residue*> neighborhood = getNeighbors(resi);
    for (int j = 0; j < neighborhood.size(); j++) {
      Residue* resj = neighborhood[j];
      if (checked[resi].find(resj) == checked[resi].end()) {
        checked[resi][resj] = true;
        real cd = contactDegree(resi, resj);
        if (cd > cdcut) {
          L.addContact(resi, resj, cd);
        }
      }
    }
  }

  return L;
}

real ConFind::weightOfAvailableRotamers(Residue* res) {
  real weight = 0;
  if (survivingRotamers.find(res) == survivingRotamers.end()) MstUtils::error("residue not cached: " + MstUtils::toString(*res), "ConFind::weightOfAvailableRotamers");
  vector<rotamerID*>& rots = survivingRotamers[res];
  for (int i = 0; i < rots.size(); i++) {
    rotamerID& rot = *(rots[i]);
    weight += aaProp[rot.aminoAcid()] * rotLib->rotamerProbability(rot);
  }
  return weight;
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

real contactList::degree(Residue* _resi, Residue* _resj) {
  if (inContact.find(_resi) == inContact.end()) return 0;
  if (inContact[_resi].find(_resj) == inContact[_resi].end()) return 0;
  return degrees[inContact[_resi][_resj]];
}
