#include "mstcondeg.h"

vector<pair<Residue*, Residue*> > contactList::getOrderedContacts() {
  vector<pair<Residue*, Residue*> > conts(orderedContacts.size());
  int i = 0;
  for (set<pair<Residue*, Residue*>, contComp>::iterator it = orderedContacts.begin(); it != orderedContacts.end(); ++it, i++) {
    conts[i] = *it;
  }
  return conts;
}

mstreal contactList::degree(Residue* _resi, Residue* _resj) {
  if (inContact.find(_resi) == inContact.end()) return 0;
  if (inContact[_resi].find(_resj) == inContact[_resi].end()) return 0;
  return degrees[inContact[_resi][_resj]];
}

bool contactList::areInContact(Residue* A, Residue* B) {
  if (inContact.find(A) == inContact.end()) return false;
  if (inContact[A].find(B) == inContact[A].end()) return false;
  return true;
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
  loCollProbCut = 0.5;
  hiCollProbCut = 2.0;
  freedomType = 2;
}

ConFind::~ConFind() {
  if (isRotLibLocal) delete rotLib;
  delete bbNN;
  delete caNN;
  for (fastmap<Residue*, vector<rotamerID*> >::iterator it = survivingRotamers.begin(); it != survivingRotamers.end(); ++it) {
    vector<rotamerID*>& rots = it->second;
    for (int i = 0; i < rots.size(); i++) delete(rots[i]);
  }
  for (fastmap<Residue*, DecoratedProximitySearch<rotamerID*>* >::iterator it = rotamerHeavySC.begin(); it != rotamerHeavySC.end(); ++it) {
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
  mstreal phi = res->getPhi(false); mstreal psi = res->getPsi(false);

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
          // backbone atoms of the same residue do not count as clashing (the
          // rotamer library should not allow true clashes with own backbone)
          if (backbone[closeOnes[ci]]->getResidue() != res) {
            prune = true;
            // clashes with ALA have a special meaning (permanent "unavoidable" contacts;
            // need to find all of them, though unlikely to have more than one)
            if (rot.isNamed("ALA")) permanentContacts[res].insert(closeOnes[ci]);
            else break;
          }
        }

        // compute contributions to this residue's interference
        set<Residue*> interfering;
        for (int ci = 0; ci < closeOnes.size(); ci++) {
          Residue* resB = backbone[closeOnes[ci]]->getResidue();
          if (interfering.find(resB) != interfering.end()) continue;
          if (resB != res) {
            interfering.insert(resB);
            interference[res][resB] += aaP * rotP/100.0;
          }
        }

        if (prune) break;
      }
      if (prune) continue;
      if (writeLog) {
        rotOut << "REM " << *res << " (" << rot.getName() << "), rotamer " << ri+1 << endl;
        Structure S(rot); S.writePDB(rotOut, "RENUMBER");
      }

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

mstreal ConFind::contactDegree(Residue* resA, Residue* resB, bool cacheA, bool cacheB, bool checkNeighbors) {
  if ((degrees.find(resA) != degrees.end()) && (degrees[resA].find(resB) != degrees[resA].end())) return degrees[resA][resB];
  if (cacheA) cache(resA);
  if (cacheB) cache(resB);
  if (checkNeighbors && !areNeighbors(resA, resB)) return 0;
  DecoratedProximitySearch<rotamerID*>* cloudA = rotamerHeavySC[resA];
  DecoratedProximitySearch<rotamerID*>* cloudB = rotamerHeavySC[resB];
  if ((cloudA == NULL) || (cloudB == NULL)) return 0.0;

  // check if the point clouds representing to two rotamer trees even overlap
  if (!cloudA->overlaps(*cloudB)) return 0;

  // if so, find rotamer pairs that clash
  fastmap<rotamerID*, fastmap<rotamerID*, bool> > clashing;
  for (int ai = 0; ai < cloudA->pointSize(); ai++) {
    vector<rotamerID*> p = cloudB->getPointsWithin(cloudA->getPoint(ai), 0, contDist);
    if (p.size() == 0) continue;
    rotamerID* rID = cloudA->getPointTag(ai);
    for (int i = 0; i < p.size(); i++) clashing[rID][p[i]] = true;
  }
  bool updateA = (updateCollProb.find(resA) != updateCollProb.end()) && updateCollProb[resA];
  bool updateB = (updateCollProb.find(resB) != updateCollProb.end()) && updateCollProb[resB];

  // compute contact degree
  mstreal cd = 0;
  for (fastmap<rotamerID*, fastmap<rotamerID*, bool> >::iterator itA = clashing.begin(); itA != clashing.end(); ++itA) {
    rotamerID* rotA = itA->first;
    mstreal rotProbA = rotLib->rotamerProbability(rotA);
    mstreal aaPropA = aaProp[rotA->aminoAcid()];
    for (fastmap<rotamerID*, bool>::iterator itB = itA->second.begin(); itB != itA->second.end(); ++itB) {
      rotamerID* rotB = itB->first;
      mstreal rotProbB = rotLib->rotamerProbability(rotB);
      mstreal aaPropB = aaProp[rotB->aminoAcid()];
      cd += aaPropA * aaPropB * rotProbA * rotProbB;
      if (updateA) collProb[resA][rotA] += aaPropB * rotProbB;
      if (updateB) collProb[resB][rotB] += aaPropA * rotProbA;
    }
  }
  cd /= weightOfAvailableRotamers(resA) * weightOfAvailableRotamers(resB);
  degrees[resA][resB] = cd;
  degrees[resB][resA] = cd;

  return cd;
}

contactList ConFind::getContacts(Residue* res, mstreal cdcut, contactList* list) {
  vector<Residue*> neighborhood = getNeighbors(res);
  cache(neighborhood);
  contactList L;

  // compute contact degree between this residue and every one of its neighbors
  if (list == NULL) list = &L;
  collProbUpdateOn(res);
  for (int i = 0; i < neighborhood.size(); i++) {
    if (res == neighborhood[i]) continue;
    mstreal cd = contactDegree(res, neighborhood[i], false, false, false);
    if (cd > cdcut) list->addContact(res, neighborhood[i], cd);
  }
  collProbUpdateOff(res);
  computeFreedom(res); // since all contacts for this residue have been visited

  return *list;
}

vector<Residue*> ConFind::getContactingResidues(Residue* res, mstreal cdcut) {
  vector<Residue*> neighborhood = getNeighbors(res);
  cache(neighborhood);

  // compute contact degree between this residue and every one of its neighbors
  vector<Residue*> partners;
  collProbUpdateOn(res);
  for (int i = 0; i < neighborhood.size(); i++) {
    if (res == neighborhood[i]) continue;
    mstreal cd = contactDegree(res, neighborhood[i], false, false, false);
    if (cd > cdcut) partners.push_back(neighborhood[i]);
  }
  collProbUpdateOff(res);
  computeFreedom(res); // since all contacts for this residue have been visited

  return partners;
}

contactList ConFind::getContacts(vector<Residue*>& residues, mstreal cdcut, contactList* list) {
  cache(residues);
  fastmap<Residue*, fastmap<Residue*, bool> > checked;
  fastmap<Residue*, bool> ofInterest;
  for (int i = 0; i < residues.size(); i++) ofInterest[residues[i]] = true;

  contactList L;
  if (list == NULL) list = &L;
  for (int i = 0; i < residues.size(); i++) {
    Residue* resi = residues[i];
    collProbUpdateOn(resi);
    vector<Residue*> neighborhood = getNeighbors(resi);
    for (int j = 0; j < neighborhood.size(); j++) {
      Residue* resj = neighborhood[j];
      if ((resi != resj) && (checked[resi].find(resj) == checked[resi].end())) {
        if (ofInterest.find(resj) != ofInterest.end()) collProbUpdateOn(resj);
        checked[resj][resi] = true;
        mstreal cd = contactDegree(resi, resj, false, true, false);
        if (cd > cdcut) {
          list->addContact(resi, resj, cd);
        }
        if (ofInterest.find(resj) != ofInterest.end()) collProbUpdateOff(resj);
      }
    }
    collProbUpdateOff(resi);
    // since all contacts for this residue have now been visited, we have all the
    // information we need to compute freedom. If the rotamers surviving at this
    // position happen not to clash with any other rotamers at all (e.g., no
    // rotamers survive at the position but could be that some survive and never
    // clash), the collision probably map for this position will not exist, so
    // make it empty.
    if (collProb.find(resi) == collProb.end()) collProb[resi] = fastmap<rotamerID*, mstreal>();
    computeFreedom(resi);
  }

  return *list;
}

contactList ConFind::getInterference(vector<Residue*>& residues, mstreal incut, contactList* list) {
  cache(residues);
  contactList L;
  if (list == NULL) list = &L;
  set<Residue*> wanted;
  for (int i = 0; i < residues.size(); i++) wanted.insert(residues[i]);

  // because interference is directional, need to check if the desired residues
  // are involved in either direction
  for (auto itA = interference.begin(); itA != interference.end(); ++itA) {
    fastmap<Residue*, mstreal>& interB = itA->second;
    bool wantA = (wanted.find(itA->first) != wanted.end());
    for (auto itB = interB.begin(); itB != interB.end(); ++itB) {
      mstreal in = itB->second;
      if (wantA || (wanted.find(itB->first) != wanted.end())) {
        if (in >= incut) {
          list->addContact(itA->first, itB->first, in, "", true);
        }
      }
    }
  }
  return *list;
}

contactList ConFind::getInterference(Structure& S, mstreal incut, contactList* list) {
  contactList L;
  if (list == NULL) list = &L;
  vector<Residue*> allRes = S.getResidues();
  getInterference(allRes, incut, list);

  return *list;
}

contactList ConFind::getContacts(Structure& S, mstreal cdcut, contactList* list) {
  contactList L;
  if (list == NULL) list = &L;
  vector<Residue*> allRes = S.getResidues();
  getContacts(allRes, cdcut, list);

  return *list;
}

mstreal ConFind::getCrowdedness(Residue* res) {
  cache(res);
  return fractionPruned[res];
}

vector<mstreal> ConFind::getCrowdedness(vector<Residue*>& residues) {
  vector<mstreal> crowdedness(residues.size());
  for (int i = 0; i < residues.size(); i++) crowdedness[i] = getCrowdedness(residues[i]);
  return crowdedness;
}

mstreal ConFind::weightOfAvailableRotamers(Residue* res) {
  mstreal weight = 0;
  if (survivingRotamers.find(res) == survivingRotamers.end()) MstUtils::error("residue not cached: " + MstUtils::toString(res), "ConFind::weightOfAvailableRotamers");
  vector<rotamerID*>& rots = survivingRotamers[res];
  for (int i = 0; i < rots.size(); i++) {
    rotamerID& rot = *(rots[i]);
    weight += aaProp[rot.aminoAcid()] * rotLib->rotamerProbability(rot);
  }
  return weight;
}

mstreal ConFind::getFreedom(Residue* res) {
  if (freedom.find(res) == freedom.end()) getContacts(res);
  return freedom[res];
}

vector<mstreal> ConFind::getFreedom(vector<Residue*>& residues) {
  vector<mstreal> freedoms(residues.size());
  for (int i = 0; i < residues.size(); i++) freedoms[i] = getFreedom(residues[i]);
  return freedoms;
}

mstreal ConFind::computeFreedom(Residue* res) {
  if (freedom.find(res) != freedom.end()) return freedom[res];
  if (collProb.find(res) == collProb.end())
    MstUtils::error("residue not cached " + MstUtils::toString(res), "ConFind::computeFreedom");

  // rotamers are in this map only if they do actually collide, so those that
  // do not collide with anything thus can be assumed to have zero collision probability mass
  fastmap<rotamerID*, mstreal>& cp = collProb[res];
  mstreal n, n1, n2;
  switch (freedomType) {
    case 1:
      // number of rotamers with < 0.5 collision probability mass
      n = survivingRotamers[res].size() - cp.size();
      for (fastmap<rotamerID*, mstreal>::iterator it = cp.begin(); it != cp.end(); ++it) {
        if (it->second/100 < 0.5) n += 1;
      }
      freedom[res] = n/numLibraryRotamers[res];
      break;
    case 2:
    case 3:
      // a combination of the number of rotamers with < 0.5 and < 2.0 collision probability masses
      n1 = survivingRotamers[res].size() - cp.size(); n2 = n1;
      for (fastmap<rotamerID*, mstreal>::iterator it = cp.begin(); it != cp.end(); ++it) {
        if (it->second/100 < loCollProbCut) n1 += 1;
        if (it->second/100 < hiCollProbCut) n2 += 1;
      }
      if (freedomType == 2) {
        freedom[res] = sqrt((n1*n1 + n2*n2)/2)/numLibraryRotamers[res];
      } else {
         n1 /= numLibraryRotamers[res]; n2 /= numLibraryRotamers[res];
         freedom[res] = sqrt((n2 * n2 + n2*n1)/2);
      }
//      freedom[res] = sqrt((n1*n1 + 0.2*n2*n2)/(1 + 0.2))/numLibraryRotamers[res];
//cout << "REPORT: " << *res << " " << n1 << " " << n2 << endl;
      break;
    default:
      MstUtils::error("unknown freedom type '" + MstUtils::toString(freedomType) + "'", "ConFind::computeFreedom");
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
    MstUtils::error("when looking for residues around " + MstUtils::toString(residue) + " did not find self!", "ConFind::getNeighbors");
  }
  return neighborhood;
}

vector<Residue*> ConFind::getNeighbors(vector<Residue*>& residues) {
  // find all residues that are within cutoff distance of the given list of residues
  fastmap<Residue*, bool> within;
  for (int k = 0; k < residues.size(); k++) {
    vector<int> close = caNN->getPointsWithin(residues[k]->findAtom("CA"), 0, dcut);
    bool foundSelf = false;
    for (int i = 0; i < close.size(); i++) {
      Residue* wres = ca[close[i]]->getResidue();
      within[wres] = true;
      if (residues[i] == wres) foundSelf = true;
    }
    if (!foundSelf) {
      MstUtils::error("when looking for residues around " + MstUtils::toString(residues[k]) + " did not find self!", "ConFind::getNeighbors(vector<Residue*>&)");
    }
  }
  vector<Residue*> neighborhood(within.size(), NULL);
  fastmap<Residue*, bool>::iterator it; int i;
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

void ConFind::openLogFile(string fname, bool append) {
  if (rotOut.is_open()) rotOut.close();
  MstUtils::openFile(rotOut, fname, append ? fstream::app : fstream::out);
}

void ConFind::closeLogFile() {
  rotOut.close();
}
