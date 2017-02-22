#include "mstcondeg.h"

vector<pair<Residue*, Residue*> > contactList::getOrderedContacts() {
  vector<pair<Residue*, Residue*> > conts(orderedContacts.size());
  int i = 0;
  for (set<pair<Residue*, Residue*>, contComp>::iterator it = orderedContacts.begin(); it != orderedContacts.end(); ++it, i++) {
    conts[i] = *it;
  }
  return conts;
}

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
  delete caNN;
  for (map<Residue*, vector<rotamerID*> >::iterator it = survivingRotamers.begin(); it != survivingRotamers.end(); ++it) {
    vector<rotamerID*>& rots = it->second;
    for (int i = 0; i < rots.size(); i++) delete(rots[i]);
  }
  for (map<Residue*, DecoratedProximitySearch<rotamerID*>* >::iterator it = rotamerHeavySC.begin(); it != rotamerHeavySC.end(); ++it) {
    if (it->second != NULL) delete it->second;
  }
}

void ConFind::init(Structure& S) {
  AtomPointerVector atoms = S.getAtoms();
  for (int i = 0; i < atoms.size(); i++) {
    Atom* a = atoms[i];
    if (a->isNamed("H")) continue;
    if (RotamerLibrary::isBackboneAtom(a)) backbone.push_back(a);
    if (a->isNamed("CA")) ca.push_back(a);
  }
  bbNN = new ProximitySearch(backbone, clashDist/2);
  caNN = new ProximitySearch(ca, dcut/2);
}

void ConFind::cache(Residue* res) {
  if (rotamerHeavySC.find(res) != rotamerHeavySC.end()) return;
  AtomPointerVector pointCloud;      // side-chain atoms of surviving rotames
  vector<rotamerID*> pointCloudTags; // corresponding tags (i.e.,  rotamer identity)
  survivingRotamers[res].resize(0);
  rotamerHeavySC[res] = NULL;
  bool writeLog = rotOut.is_open();

  // make sure this residue has a proper backbone, otherwise adding rotamers will fail
  if (!res->atomExists("N") || !res->atomExists("CA") || !res->atomExists("C")) {
    MstUtils::error("cannot build rotamers at position " + MstUtils::toString(*res) + " as it lacks proper backbone!", "ConFind::cache(Residue*)");
  }
  real phi = res->getPhi(); real psi = res->getPsi();

  // load rotamers of each amino acid
  int numRemRotsInPosition = 0; int totNumRotsInPosition = 0;
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
      bool prune = false;
      vector<int> closeOnes;
      for (int k = 0; k < rot.atomSize(); k++) {
        Atom& a = rot[k];
        if (!countsAsSidechain(a)) continue;

        // shuld the rotamer be pruned based on this atom's clash(es)?
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
      if (writeLog) { Structure S(rot); S.writePDB(rotOut); }

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
  numLibraryRotamers[res] = totNumRotsInPosition;

  // cash all atoms for faster distance-based searches
  if (pointCloud.size() != 0) rotamerHeavySC[res] = new DecoratedProximitySearch<rotamerID*>(pointCloud, contDist/2, pointCloudTags);
  pointCloud.deletePointers();
}

bool ConFind::countsAsSidechain(Atom& a) {
  if (RotamerLibrary::isHydrogen(a) || RotamerLibrary::isBackboneAtom(a)) return false;
  if (doNotCountCB && a.isNamed("CB") && !(a.getResidue()->isNamed("ALA"))) return false;
  return true;
}

void ConFind::cache(vector<Residue*>& residues) {
  for (int i = 0; i < residues.size(); i++) cache(residues[i]);
}

void ConFind::cache(Structure& S) {
  vector<Residue*> residues = S.getResidues();
  cache(residues);
}

real ConFind::contactDegree(Residue* resA, Residue* resB, bool doNotCache, bool checkNeighbors) {
  if ((degrees.find(resA) != degrees.end()) && (degrees[resA].find(resB) != degrees[resA].end())) return degrees[resA][resB];
  if (!doNotCache) { cache(resA); cache(resB); }
  if (checkNeighbors && !areNeighbors(resA, resB)) return 0;
  DecoratedProximitySearch<rotamerID*>* cloudA = rotamerHeavySC[resA];
  DecoratedProximitySearch<rotamerID*>* cloudB = rotamerHeavySC[resB];
  if (cloudA == NULL) collProb[resA] = map<rotamerID*, real>();
  if (cloudB == NULL) collProb[resB] = map<rotamerID*, real>();
  if ((cloudA == NULL) || (cloudB == NULL)) return 0.0;

  // check if the point clouds representing to two rotamer trees even overlap
  if (!cloudA->overlaps(*cloudB)) return 0;

  // if so, find rotamer pairs that clash
  map<rotamerID*, map<rotamerID*, bool> > clashing;
  for (int ai = 0; ai < cloudA->pointSize(); ai++) {
    vector<rotamerID*> p = cloudB->getPointsWithin(cloudA->getPoint(ai), 0, contDist);
    if (p.size() == 0) continue;
    rotamerID* rID = cloudA->getPointTag(ai);
    for (int i = 0; i < p.size(); i++) clashing[rID][p[i]] = true;
  }

  // compute contact degree
  real cd = 0;
  for (map<rotamerID*, map<rotamerID*, bool> >::iterator itA = clashing.begin(); itA != clashing.end(); ++itA) {
    rotamerID* rotA = itA->first;
    real rotProbA = rotLib->rotamerProbability(rotA);
    real aaPropA = aaProp[rotA->aminoAcid()];
    for (map<rotamerID*, bool>::iterator itB = itA->second.begin(); itB != itA->second.end(); ++itB) {
      rotamerID* rotB = itB->first;
      real rotProbB = rotLib->rotamerProbability(rotB);
      real aaPropB = aaProp[rotB->aminoAcid()];
      cd +=  aaPropA * aaPropB * rotProbA * rotProbB;
      collProb[resA][rotA] += aaPropA * rotProbA;
      collProb[resB][rotB] += aaPropB * rotProbB;
    }
  }
  cd /= weightOfAvailableRotamers(resA) * weightOfAvailableRotamers(resB);
  degrees[resA][resB] = cd;
  degrees[resB][resA] = cd;

  return cd;
}

contactList ConFind::getContacts(Residue* res, real cdcut, contactList* list) {
  vector<Residue*> neighborhood = getNeighbors(res);
  cache(neighborhood);
  contactList L;

  // compute contact degree between this residue and every one of its neighbors
  if (list == NULL) list = &L;
  for (int i = 0; i < neighborhood.size(); i++) {
    if (res == neighborhood[i]) continue;
    real cd = contactDegree(res, neighborhood[i], true, false);
    if (cd > cdcut) list->addContact(res, neighborhood[i], cd);
  }
  computeFreedom(res); // since all contacts for this residue have been visited

  return *list;
}

vector<Residue*> ConFind::getContactingResidues(Residue* res, real cdcut) {
  vector<Residue*> neighborhood = getNeighbors(res);
  cache(neighborhood);

  // compute contact degree between this residue and every one of its neighbors
  vector<Residue*> partners;
  for (int i = 0; i < neighborhood.size(); i++) {
    if (res == neighborhood[i]) continue;
    real cd = contactDegree(res, neighborhood[i], true, false);
    if (cd > cdcut) partners.push_back(neighborhood[i]);
  }
  computeFreedom(res); // since all contacts for this residue have been visited

  return partners;
}

contactList ConFind::getContacts(vector<Residue*>& residues, real cdcut, contactList* list) {
  cache(residues);
  map<Residue*, map<Residue*, bool> > checked;

  contactList L;
  if (list == NULL) list = &L;
  for (int i = 0; i < residues.size(); i++) {
    Residue* resi = residues[i];
    vector<Residue*> neighborhood = getNeighbors(resi);
    for (int j = 0; j < neighborhood.size(); j++) {
      Residue* resj = neighborhood[j];
      if ((resi != resj) && (checked[resi].find(resj) == checked[resi].end())) {
        checked[resj][resi] = true;
        real cd = contactDegree(resi, resj, false, false);
        if (cd > cdcut) {
          list->addContact(resi, resj, cd);
        }
      }
    }
    computeFreedom(resi); // since all contacts for this residue have been visited
  }

  return *list;
}

contactList ConFind::getContacts(Structure& S, real cdcut, contactList* list) {
  contactList L;
  if (list == NULL) list = &L;
  vector<Residue*> allRes = S.getResidues();
  getContacts(allRes, cdcut, list);

  return *list;
}

real ConFind::getCrowdedness(Residue* res) {
  cache(res);
  return fractionPruned[res];
}

vector<real> ConFind::getCrowdedness(vector<Residue*>& residues) {
  vector<real> crowdedness(residues.size());
  for (int i = 0; i < residues.size(); i++) crowdedness[i] = getCrowdedness(residues[i]);
  return crowdedness;
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

real ConFind::getFreedom(Residue* res) {
  if (freedom.find(res) == freedom.end()) getContacts(res);
  return freedom[res];
}

vector<real> ConFind::getFreedom(vector<Residue*>& residues) {
  vector<real> freedoms(residues.size());
  for (int i = 0; i < residues.size(); i++) freedoms[i] = getFreedom(residues[i]);
  return freedoms;
}

real ConFind::computeFreedom(Residue* res) {
  if (freedom.find(res) != freedom.end()) return freedom[res];
  if (collProb.find(res) == collProb.end())
    MstUtils::error("residue not cached " + MstUtils::toString(*res), "ConFind::computeFreedom");

  int type = 2;
  map<rotamerID*, real>& cp = collProb[res];
  switch (type) {
    case 1: {
      // number of rotamers with < 0.5 collision probability mass
      real n = 0;
      for (map<rotamerID*, real>::iterator it = cp.begin(); it != cp.end(); ++it) {
        if (it->second/100 < 0.5) n += 1;
      }
      freedom[res] = n/numLibraryRotamers[res];
      break;
    }
    case 2: {
      // a combination of the number of rotamers with < 0.5 and < 2.0 collision probability masses
      real n1 = 0; real n2 = 0;
      for (map<rotamerID*, real>::iterator it = cp.begin(); it != cp.end(); ++it) {
        if (it->second/100 < 0.5) n1 += 1;
        if (it->second/100 < 2.0) n2 += 1;
      }
      freedom[res] = sqrt((n1*n1 + n2*n2)/2)/numLibraryRotamers[res];
      break;
    }
    default:
      MstUtils::error("unknown freedom type '" + MstUtils::toString(type) + "'", "ConFind::computeFreedom");
  }
  return freedom[res];
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

bool ConFind::areNeighbors(Residue* resA, Residue* resB) {
  Atom* CAA = resA->findAtom("CA");
  Atom* CAB = resB->findAtom("CA");
  return (CAA->distance(*CAB) <= dcut);
}

void ConFind::openLogFile(string fname) {
  if (rotOut.is_open()) rotOut.close();
  MstUtils::openFile(rotOut, fname, fstream::app);
}

void ConFind::closeLogFile() {
  rotOut.close();
}

real contactList::degree(Residue* _resi, Residue* _resj) {
  if (inContact.find(_resi) == inContact.end()) return 0;
  if (inContact[_resi].find(_resj) == inContact[_resi].end()) return 0;
  return degrees[inContact[_resi][_resj]];
}
