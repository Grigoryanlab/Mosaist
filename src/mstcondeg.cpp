#include "mstcondeg.h"

/* contactList */

void contactList::sortByDegree() {
  vector<int> sortedIndices = MstUtils::sortIndices(degrees, true);
  vector<Residue*> resiOld = resi, resjOld = resj;
  vector<mstreal> degreesOld = degrees;
  vector<string> infosOld = infos;
  for (int i = 0; i < sortedIndices.size(); i++) {
    resi[i] = resiOld[sortedIndices[i]];
    resj[i] = resjOld[sortedIndices[i]];
    degrees[i] = degreesOld[sortedIndices[i]];
    infos[i] = infosOld[sortedIndices[i]];
  }
}

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
  return degrees[*inContact[_resi][_resj].begin()];
}

bool contactList::areInContact(Residue* _resi, Residue* _resj) {
  if (inContact.find(_resi) == inContact.end()) return false;
  if (inContact[_resi].find(_resj) == inContact[_resi].end()) return false;
  return true;
}

/* ConFind */

ConFind::ConFind(string rotLibFile, const Structure& S, bool _strict) {
  setParams();
  rotLib = new RotamerLibrary(rotLibFile);
  isRotLibLocal = true;
  strict = _strict;
  init(S);
}
ConFind::ConFind(RotamerLibrary* _rotLib, const Structure& S, bool _strict) {
  setParams();
  rotLib = _rotLib;
  isRotLibLocal = false;
  strict = _strict;
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
  for (auto it = survivingRotamers.begin(); it != survivingRotamers.end(); ++it) {
    vector<rotamerID*>& rots = it->second;
    for (int i = 0; i < rots.size(); i++) delete(rots[i]);
  }
  for (auto res_it = rotamerHeavySC.begin(); res_it != rotamerHeavySC.end(); ++res_it) {
    for (auto aa_it = rotamerHeavySC[res_it->first].begin(); aa_it != rotamerHeavySC[res_it->first].end(); ++aa_it) {
      if (aa_it->second != NULL) delete(aa_it->second);
    }
  }
}

void ConFind::init(const Structure& S) {
  AtomPointerVector atoms = S.getAtoms();
  for (int i = 0; i < atoms.size(); i++) {
    Atom* a = atoms[i];
    if (RotamerLibrary::isHydrogen(a)) continue;
    int idx = RotamerLibrary::backboneAtomType(a);
    if (idx >= 0) backbone.push_back(a);
    if (idx == RotamerLibrary::bbCA) ca.push_back(a);
  }
  bbNN = new ProximitySearch(backbone, clashDist/2);
  caNN = new ProximitySearch(ca, dcut/2);
}

void ConFind::cache(Residue* res) {
  string res_name = res->getName();
  if (rotamerHeavySC.find(res) != rotamerHeavySC.end()) return;
  AtomPointerVector pointCloud;      // side-chain atoms of surviving rotames
  vector<rotamerID*> pointCloudTags; // corresponding tags (i.e.,  rotamer identity)
  survivingRotamers[res].resize(0);
  bool writeLog = rotOut.is_open();

  // make sure this residue has a proper backbone, otherwise adding rotamers will fail
  vector<Atom*> bb = RotamerLibrary::getBackbone(res);
  if ((bb[RotamerLibrary::bbN] == NULL) || (bb[RotamerLibrary::bbCA] == NULL) || (bb[RotamerLibrary::bbC] == NULL)) {
    MstUtils::error("cannot build rotamers at position " + MstUtils::toString(*res) + " as it lacks proper backbone!", "ConFind::cache(Residue*)");
  }
  mstreal phi = res->getPhi(false); mstreal psi = res->getPsi(false);

  // load rotamers of each amino acid
  int numRemRotsInPosition = 0; int totNumRotsInPosition = 0;
  for (string aa : aaNames) {
    if (aaProp.find(aa) == aaProp.end()) MstUtils::error("no propensity defined for amino acid " + aa);
    rotamerHeavySC[res][aa] = NULL;
    double aaP = aaProp[aa];
    if (strict && res_name != aa && res_name != "UNK") continue;
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

        // should the rotamer be pruned based on this atom's clash(es)?
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
            if (interference[res][resB].count(aa) == 0) interference[res][resB][aa] = 0.0;
            interference[res][resB][aa] += aaP * rotP/100.0;
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
    
    // cash all the rotamer heavy atoms from rotamers of this amino acid for faster distance-based searches
    if (pointCloud.size() != 0) rotamerHeavySC[res][aa] = new DecoratedProximitySearch<rotamerID*>(pointCloud, contDist/2, pointCloudTags);
    pointCloud.deletePointers();
    pointCloudTags.clear();
    
    totNumRotsInPosition += nr;
  }
  fractionPruned[res] = (totNumRotsInPosition - numRemRotsInPosition)*1.0/totNumRotsInPosition;
  numLibraryRotamers[res] = totNumRotsInPosition;
}

bool ConFind::countsAsSidechain(Atom& a) {
  if (RotamerLibrary::isHydrogen(a) || RotamerLibrary::isBackboneAtom(a)) return false;
  if (doNotCountCB && a.isNamed("CB") && !(a.getResidue()->isNamed("ALA"))) return false;
  return true;
}

void ConFind::cache(const vector<Residue*>& residues) {
  for (int i = 0; i < residues.size(); i++) cache(residues[i]);
}

void ConFind::cache(const Structure& S) {
  vector<Residue*> residues = S.getResidues();
  cache(residues);
}

mstreal ConFind::contactDegree(Residue* resA, Residue* resB, bool cacheA, bool cacheB, bool checkNeighbors, set<string> aaAllowedA, set<string> aaAllowedB) {
  // only cache CD value/collision probabities if amino acids are not restricted at either position
  bool no_aa_restriction = (aaAllowedA.empty() && aaAllowedA.empty());
  if (aaAllowedA.empty()) aaAllowedA = aaNames;
  if (aaAllowedB.empty()) aaAllowedB = aaNames;
  bool updateA = (no_aa_restriction && updateCollProb.find(resA) != updateCollProb.end()) && updateCollProb[resA];
  bool updateB = (no_aa_restriction && updateCollProb.find(resB) != updateCollProb.end()) && updateCollProb[resB];
  // if collision probabilities need to be updated, then go into CD calculation even if value already available
  if (no_aa_restriction && (degrees.find(resA) != degrees.end()) && (degrees[resA].find(resB) != degrees[resA].end()) && !updateA && !updateB) {
    return degrees[resA][resB];
  }
  if (cacheA) cache(resA);
  if (cacheB) cache(resB);
  if (checkNeighbors && !areNeighbors(resA, resB)) return 0;
  
  // get point clouds for each amino acid from the A/B sets and find interacting rotamer pairs
  fastmap<rotamerID*, fastmap<rotamerID*, bool> > clashing;
  DecoratedProximitySearch<rotamerID*>* cloudA;
  DecoratedProximitySearch<rotamerID*>* cloudB;
  for (string resA_aa : aaAllowedA) {
    if (aaNames.find(resA_aa) == aaNames.end()) MstUtils::error("amino acid with the name: "+resA_aa+" not in list of allowable amino acids");
    if (rotamerHeavySC[resA].count(resA_aa) != 0) cloudA = rotamerHeavySC[resA][resA_aa];
    else continue;
    if (cloudA == NULL) continue;
    
    for (string resB_aa : aaAllowedB) {
      if (aaNames.find(resB_aa) == aaNames.end()) MstUtils::error("amino acid with the name: "+resB_aa+" not in list of allowable amino acids");
      if (rotamerHeavySC[resB].count(resB_aa) != 0) cloudB = rotamerHeavySC[resB][resB_aa];
      else continue;
      
      if (cloudB == NULL) continue;

      // check if the point clouds representing the two rotamer trees even overlap
      if (!cloudA->overlaps(*cloudB,contDist)) continue;
      
      // if so, find rotamer pairs that clash
      // NOTE: this is the slow part of contact finding; could perhaps speed up by
      // storing smaller ProximitySearch object (not decorated) for each rotamer
      vector<rotamerID*> p;
      for (int ai = 0; ai < cloudA->pointSize(); ai++) {
        cloudB->getPointsWithin(cloudA->getPoint(ai), 0, contDist, &p);
        if (p.size() == 0) continue;
        rotamerID* rID = cloudA->getPointTag(ai);
        for (int i = 0; i < p.size(); i++) clashing[rID][p[i]] = true;
      }
    }
  }
  
  // compute contact degree
  mstreal cd = 0.0;
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
  
  // in case there are no available rotamer pairs, avoid division and just set to 0.0
  mstreal denom = weightOfAvailableRotamers(resA,aaAllowedA) * weightOfAvailableRotamers(resB,aaAllowedB);
  if (denom == 0.0) cd = 0.0;
  else cd /= denom;
  if (no_aa_restriction) {
    degrees[resA][resB] = cd;
    degrees[resB][resA] = cd;
  }

  return cd;
}

contactList ConFind::getContacts(Residue* res, mstreal cdcut, contactList* list) {
  // this way, the contact computing code lives only in one place (minimal cost)
  return getContacts(vector<Residue*>(1, res), cdcut, list);
}

vector<Residue*> ConFind::getContactingResidues(Residue* res, mstreal cdcut) {
  // this way, the contact computing code lives only in one place (minimal cost)
  contactList list = getContacts(vector<Residue*>(1, res), cdcut);
  vector<Residue*> partners(list.size(), NULL);
  collProbUpdateOn(res);
  for (int i = 0; i < list.size(); i++) partners[i] = list.residueB(i);
  return partners;
}

contactList ConFind::getContacts(Structure& S, mstreal cdcut, contactList* list) {
  contactList L;
  if (list == NULL) list = &L;
  vector<Residue*> allRes = S.getResidues();
  getContacts(allRes, cdcut, list);
  
  return *list;
}

contactList ConFind::getContacts(const vector<Residue*>& residues, mstreal cdcut, contactList* list) {
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

contactList ConFind::getConstrainedContacts(const vector<Residue *> &residues, mstreal cdcut, contactList* list) {
  contactList L;
  if (list == NULL) list = &L;
  
  cache(residues);
  
  for (int i = 0; i < residues.size(); i++) {
    Residue* resi = residues[i];
    vector<Residue*> neighborhood = getNeighbors(resi);
    for (int j = 0; j < neighborhood.size(); j++) {
      Residue* resj = neighborhood[j];
      if (resi != resj) {
        mstreal cd;
        for (string aa: aaNames) {
          cd = contactDegree(resi, resj, false, true, false, {aa}, aaNames);
          if (cd > cdcut) {
            list->addContact(resi, resj, cd, "", true, {aa}, aaNames);
          }
        }
      }
    }
  }
  return L;
}

contactList ConFind::getInterference(const vector<Residue*>& residues, mstreal incut, contactList* list) {
  cache(residues);
  contactList L;
  if (list == NULL) list = &L;
  set<Residue*> wanted;
  for (int i = 0; i < residues.size(); i++) wanted.insert(residues[i]);

  // because interference is directional, need to check if the desired residues
  // are involved in either direction
  for (auto itA = interference.begin(); itA != interference.end(); ++itA) {
    fastmap<Residue*, fastmap<string,mstreal>>& interB = itA->second;
    bool wantA = (wanted.find(itA->first) != wanted.end());
    for (auto itB = interB.begin(); itB != interB.end(); ++itB) {
      if (wantA || (wanted.find(itB->first) != wanted.end())) {
        mstreal in = interferenceValue(itA->first, itB->first);
        if (in >= incut) {
          list->addContact(itA->first, itB->first, in, "", true);
        }
      }
    }
  }
  return *list;
}

contactList ConFind::getInterference(const Structure& S, mstreal incut, contactList* list) {
  contactList L;
  if (list == NULL) list = &L;
  vector<Residue*> allRes = S.getResidues();
  getInterference(allRes, incut, list);

  return *list;
}

contactList ConFind::getInterfering(const vector<Residue*>& residues, mstreal incut, contactList* list) {
  cache(residues);
  contactList L;
  if (list == NULL) list = &L;
  set<Residue*> wanted;
  for (int i = 0; i < residues.size(); i++) wanted.insert(residues[i]);

  // check only in one direction, where one the specified residues is the source
  for (int i = 0; i < residues.size(); i++) {
    Residue* resA = residues[i];
    if (interference.find(resA) == interference.end()) continue;
    fastmap<Residue*, fastmap<string, mstreal>>& interB = interference[resA];
    for (auto itB = interB.begin(); itB != interB.end(); ++itB) {
      mstreal in = interferenceValue(resA, itB->first);
      if (in >= incut) {
        list->addContact(resA, itB->first, in, "", true);
      }
    }
  }
  return *list;
}

contactList ConFind::getInterfering(const Structure& S, mstreal incut, contactList* list) {
  contactList L;
  if (list == NULL) list = &L;
  vector<Residue*> allRes = S.getResidues();
  getInterference(allRes, incut, list);

  return *list;
}

mstreal ConFind::interferenceValue(Residue* resA, Residue* resB, set<string> aaAllowed) {
  mstreal in = 0.0;
  if (interference[resA].count(resB) == 0) return in;
  if (aaAllowed.empty()) aaAllowed = aaNames;
  
  // sum up interference of resB backbone on resA sidechain (considering only allowed amino acids)
  for (string aa : aaAllowed) {
    if (aaNames.find(aa) == aaNames.end()) MstUtils::error("amino acid with the name: "+aa+" not in list of allowable amino acids");
    in += interference[resA][resB][aa];
  }
  
  // normalize the interference by the weight of available rotamers at resA
  if (weightOfAvailableAminoAcids(aaAllowed) == 0.0) in = 0.0;
  else in /= weightOfAvailableAminoAcids(aaAllowed);
  return in;
}

mstreal ConFind::bbInteraction(Residue *resA, Residue *resB) {
  /* Get pointers to all of the ResA/ResB backbone atoms and then for each pair of atoms
   between the two sets calculate distance. Report the distance of the closest pair of atoms.
   */
  vector<Atom*> resA_bb = RotamerLibrary::getBackbone(resA);
  vector<Atom*> resB_bb = RotamerLibrary::getBackbone(resB);
  mstreal min_dist = 0; //initialized to keep the compiler happy
  mstreal curr_dist;
  for (int i = 0; i < resA_bb.size(); i++) {
    for (int j = 0; j < resB_bb.size(); j++) {
      if (resA_bb[i] == NULL || resB_bb[j] == NULL) continue; //in the case that an O is missing from the backbone
      curr_dist = resA_bb[i]->distance(resB_bb[j]);
      if (i == 0 || curr_dist < min_dist) min_dist = curr_dist;
    }
  }
  return min_dist;
}

contactList ConFind::getBBInteraction(Residue* res, mstreal dcut, int ignoreFlanking, contactList* list) {
  return getBBInteraction(vector<Residue*>(1, res), dcut, ignoreFlanking, list);
}

contactList ConFind::getBBInteraction(Structure& S, mstreal dcut, int ignoreFlanking, contactList* list) {
  contactList L;
  if (list == NULL) list = &L;
  vector<Residue*> allRes = S.getResidues();
  return getBBInteraction(allRes, dcut, ignoreFlanking, list);
}

contactList ConFind::getBBInteraction(const vector<Residue*>& residues, mstreal dcut, int ignoreFlanking, contactList* list) {
  fastmap<Residue*, fastmap<Residue*, bool> > checked;
  
  contactList L;
  if (list == NULL) list = &L;
  for (int i = 0; i < residues.size(); i++) {
    Residue* resi = residues[i];
    vector<Residue*> neighborhood = getNeighbors(resi);
    for (int j = 0; j < neighborhood.size(); j++) {
      Residue* resj = neighborhood[j];
      // Check if the residues are 1) on the same chain and 2) within the ignore range
      bool ignore = false;
      if (resi->getChain() == resj->getChain()) {
        int distance = abs(resi->getResidueIndex()-resj->getResidueIndex());
        if (distance <= ignoreFlanking) ignore = true;
      }
      if ((!ignore) && (resi != resj) && (checked[resi].find(resj) == checked[resi].end())) {
        checked[resj][resi] = true;
        mstreal dist = bbInteraction(resi, resj);
        if (dist < dcut) {
            list->addContact(resi, resj, dist);
        }
      }
    }
  }
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

mstreal ConFind::weightOfAvailableRotamers(Residue* res, set<string> available_aa) {
  mstreal weight = 0;
  if (survivingRotamers.find(res) == survivingRotamers.end()) MstUtils::error("residue not cached: " + MstUtils::toString(res), "ConFind::weightOfAvailableRotamers");
  vector<rotamerID*>& rots = survivingRotamers[res];
  for (int i = 0; i < rots.size(); i++) {
    rotamerID& rot = *(rots[i]);
    if (available_aa.find(rot.aminoAcid()) == available_aa.end()) continue;
    weight += aaProp[rot.aminoAcid()] * rotLib->rotamerProbability(rot);
  }
  return weight;
}

mstreal ConFind::weightOfAvailableAminoAcids(set<string> available_aa) {
  mstreal weight = 0;
  for (string aa : available_aa) {
    weight += aaProp[aa];
  }
  return weight/100.0;
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
  if (collProb.find(res) == collProb.end()) {
    MstUtils::error("residue not cached", "ConFind::computeFreedom");
  }

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
    vector<Residue*> neigh = getNeighbors(residues[k]);
    for (int i = 0; i < neigh.size(); i++) within[neigh[i]] = true;
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
