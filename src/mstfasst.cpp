#include "mstfasst.h"

/* --------- FASST::optList --------- */
void FASST::optList::setOptions(const vector<mstreal>& _costs, bool add) {
  // sort costs and keep track of indices to know the rank-to-index mapping
  costs = _costs;
  rankToIdx.resize(_costs.size());
  for (int i = 0; i < rankToIdx.size(); i++) rankToIdx[i] = i;
  sort(rankToIdx.begin(), rankToIdx.end(), [this](int i, int j) { return costs[i] < costs[j]; });
  for (int i = 0; i < rankToIdx.size(); i++) costs[i] = _costs[rankToIdx[i]];

  // sort the rank-to-index mapping and keep track of indices to know the index-to-rank mapping
  idxToRank.resize(costs.size());
  for (int i = 0; i < idxToRank.size(); i++) idxToRank[i] = i;
  sort(idxToRank.begin(), idxToRank.end(), [this](int i, int j) { return rankToIdx[i] < rankToIdx[j]; });

  // either include or exclude all options to start off, as instructed
  isIn.resize(costs.size(), add);
  if (add) {
    bestCostRank = 0;
    numIn = costs.size();
  } else {
    bestCostRank = -1;
    numIn = 0;
  }
}

void FASST::optList::addOption(int k) {
  if (!isIn[k]) {
    numIn++;
    isIn[k] = true;
    // update best, if the newly inserted option is better
    if ((bestCostRank < 0) || (costs[bestCostRank] > costs[idxToRank[k]])) {
      bestCostRank = idxToRank[k];
    }
  }
}

bool FASST::optList::consistencyCheck() {
  int nn = 0;
  for (int i = 0; i < isIn.size(); i++) {
    if (isIn[i]) nn++;
  }
  if (nn != numIn) {
    return false;
  }
  if (numIn == 0) {
    if (bestCostRank != -1) return false;
  } else {
    if ((bestCostRank < 0) || (bestCostRank >= costs.size())) return false;
    if (!isIn[rankToIdx[bestCostRank]]) return false;
  }
  return true;
}

void FASST::optList::removeOption(int k) {
  // if ((k < 0) || (k >= costs.size())) MstUtils::error("out-of-range index specified: " + MstUtils::toString(k), "FASST::optList::removeOption(int)");
  if (isIn[k]) {
    numIn--;
    isIn[k] = false;
    // find new best, if this was the best
    if (bestCostRank == idxToRank[k]) {
      int oldBestRank = bestCostRank;
      bestCostRank = -1;
      for (int r = oldBestRank + 1; r < costs.size(); r++) {
        if (isIn[rankToIdx[r]]) { bestCostRank = r; break; }
      }
    }
  }
}

void FASST::optList::removeOptions(int b, int e) {
  if (b < 0) b = 0;
  if (e >= totNumOptions()) e = totNumOptions() - 1;
  for (int k = b; k <= e; k++) removeOption(k);
}

void FASST::optList::removeAllOptions() {
  for (int i = 0; i < isIn.size(); i++) isIn[i] = false;
  numIn = 0;
  bestCostRank = -1;
}

void FASST::optList::copyIn(const FASST::optList& opt) {
  isIn = opt.isIn;
  bestCostRank = opt.bestCostRank;
  numIn = opt.numIn;
}

void FASST::optList::intersectOptions(const vector<int>& opts) {
  vector<bool> isOK(isIn.size(), false);
  for (int i = 0; i < opts.size(); i++) isOK[opts[i]] = true;
  for (int i = 0; i < isIn.size(); i++) {
    if (isIn[i] && !isOK[i]) removeOption(i);
  }
}

void FASST::optList::intersectOptions(const vector<bool>& isOK) {
  for (int i = 0; i < isIn.size(); i++) {
    if (isIn[i] && !isOK[i]) removeOption(i);
  }
}

void FASST::optList::constrainLE(int idx) {
  for (int i = MstUtils::max(idx+1, 0); i < isIn.size(); i++) removeOption(i);
}

void FASST::optList::constrainGE(int idx) {
  int UB = MstUtils::min(idx, (int) isIn.size());
  for (int i = 0; i < UB; i++) removeOption(i);
}

void FASST::optList::constrainRange(int idxLow, int idxHigh) {
  constrainGE(idxLow);
  constrainLE(idxHigh);
}

/* --------- fasstSearchOptions ---------- */
void fasstSearchOptions::setMaxNumMatches(int _max) {
  maxNumMatches = _max;
  if (!areNumMatchConstraintsConsistent()) MstUtils::error("invalid combination of match number constraints: [min, max, sufficient] = [" + MstUtils::toString(minNumMatches) + ", " + MstUtils::toString(maxNumMatches) + ", " + MstUtils::toString(suffNumMatches) + "]", "FASST::setMaxNumMatches");
}

void fasstSearchOptions::setMinNumMatches(int _min) {
  minNumMatches = _min;
  if (!areNumMatchConstraintsConsistent()) MstUtils::error("invalid combination of match number constraints: [min, max, sufficient] = [" + MstUtils::toString(minNumMatches) + ", " + MstUtils::toString(maxNumMatches) + ", " + MstUtils::toString(suffNumMatches) + "]", "FASST::setMinNumMatches");
}

void fasstSearchOptions::setSufficientNumMatches(int _suff) {
  suffNumMatches = _suff;
  if (!areNumMatchConstraintsConsistent()) MstUtils::error("invalid combination of match number constraints: [min, max, sufficient] = [" + MstUtils::toString(minNumMatches) + ", " + MstUtils::toString(maxNumMatches) + ", " + MstUtils::toString(suffNumMatches) + "]", "FASST::setSufficientNumMatches");
}

bool fasstSearchOptions::areNumMatchConstraintsConsistent() const {
  if (isMaxNumMatchesSet() && isMinNumMatchesSet() && (minNumMatches > maxNumMatches)) return false;
  if (isMaxNumMatchesSet() && isSufficientNumMatchesSet() && (maxNumMatches < suffNumMatches)) return false;
  if (isMinNumMatchesSet() && isSufficientNumMatchesSet() && (minNumMatches > suffNumMatches)) return false;
  return true;
}

void fasstSearchOptions::setMinGap(int i, int j, int gapLim) {
  if (gapLim < 0) MstUtils::error("gap constraints must be defined in the non-negative direction", "FASST::setMinGap");
  minGap[i][j] = gapLim;
  minGapSet[i][j] = true;
  gapConstSet = true;
}

void fasstSearchOptions::setMaxGap(int i, int j, int gapLim) {
  if (gapLim < 0) MstUtils::error("gap constraints must be defined in the non-negative direction", "FASST::setMinGap");
  maxGap[i][j] = gapLim;
  maxGapSet[i][j] = true;
  gapConstSet = true;
}

void fasstSearchOptions::resetGapConstraints(int numQuerySegs) {
  minGap.resize(numQuerySegs, vector<int>(numQuerySegs, 0));
  maxGap = minGap;
  minGapSet.resize(numQuerySegs, vector<bool>(numQuerySegs, false));
  maxGapSet = minGapSet;
  gapConstSet = false;
}

bool fasstSearchOptions::validateGapConstraints(int numQuerySegs) const {
  if ((numQuerySegs != minGap.size()) || (numQuerySegs != maxGap.size())) {
    MstUtils::error("gap constraints inconsistent with number of segments in query", "fasstSearchOptions::areGapsConsistentConsistent");
  }
  return true;
}

bool fasstSearchOptions::validateSearchRequest(int numQuerySegs) const {
  return validateGapConstraints(numQuerySegs);
}

/* --------- FASST --------- */
FASST::FASST() {
  recLevel = 0;
  opts.setRMSDCutoff(1.0);
  setSearchType(searchType::FULLBB);
  querySize = 0;
  updateGrids = false;
  gridSpacing = 15.0;
}

FASST::~FASST() {
  // need to delete atoms only on the lowest level of recursion, because at
  // higher levels we point to the same atoms
  if (targetMasks.size()) targetMasks.back().deletePointers();
  for (int i = 0; i < targetStructs.size(); i++) {
    if (targetStructs[i]) delete targetStructs[i];
    else targets[i].deletePointers();
  }
  for (int i = 0; i < ps.size(); i++) delete ps[i];
}

void FASST::setCurrentRMSDCutoff(mstreal cut, int p) {
  rmsdCut = cut;
  residualCut = rmsdCut*rmsdCut*querySize;
  rPrior = p;
  if (p >= 0) {
    if (p >= rmsdCutTemp.size()) rmsdCutTemp.resize(p + 1, -1);
    rmsdCutTemp[p] = cut;
  } else {
    rmsdCutDef = cut;
  }
}

void FASST::resetCurrentRMSDCutoff(int p) {
  if (p != rPrior) {
    for (int i = rPrior; (i > p) && (i >= 0); i--) rmsdCutTemp[i] = -1; // wipe out all RMSD at levels below the given one
    for (; p >= 0; p--) { // find the first priority level at/before the given one that has its RMSD set
      if (rmsdCutTemp[p] >= 0) break;
    }
    rmsdCut = (p < 0) ? rmsdCutDef : rmsdCutTemp[p];
    residualCut = rmsdCut*rmsdCut*querySize;
    rPrior = p;
  }
}

/* This function sets up the query, in the process deciding which part of the
 * query is really searchable (e.g., backbone). Various search types an be added
 * in the future to provide search capabilities over different parts of the
 * structure. The current implementation stipulates that the searchable part
 * involve the same number of atoms for each residue, so that atoms can be
 * stored in flat arrays for efficiency. But these could presumably be
 * interpreted differently. E.g., only some residue matches (parents of atoms)
 * may be accepted or some of the atoms can be dummy/empty ones. */
void FASST::setQuery(const string& pdbFile, bool autoSplitChains) {
  queryStruct.reset();
  queryStruct.readPDB(pdbFile);
  if (autoSplitChains) queryStruct = queryStruct.reassignChainsByConnectivity();
  processQuery();
}

void FASST::setQuery(const Structure& Q, bool autoSplitChains) {
  queryStruct = Q;
  if (autoSplitChains) queryStruct = queryStruct.reassignChainsByConnectivity();
  processQuery();
}

void FASST::processQuery() {
  querySize = 0;
  // auto-splitting segments by connectivity, but can do differently
  if (queryStruct.chainSize() == 0) MstUtils::error("query should not be an empty structure", "FASST::processQuery");
  query.resize(queryStruct.chainSize());
  for (int i = 0; i < queryStruct.chainSize(); i++) {
    query[i].resize(0);
    if (!parseChain(queryStruct[i], &(query[i]))) {
      MstUtils::error("could not set query, because some atoms for the specified search type were missing", "FASST::processQuery");
    }
    MstUtils::assertCond(query[i].size() > 0, "query contains empty segment(s)", "FASST::processQuery");
    querySize += query[i].size();
  }
  currAlignment.resize(query.size(), -1);

  // re-order query segments by length (longest first)
  queryOrig = query;
  qSegOrd.resize(query.size());
  for (int i = 0; i < qSegOrd.size(); i++) qSegOrd[i] = i;
  sort(qSegOrd.begin(), qSegOrd.end(), [this](size_t i, size_t j) {return query[i].size() > query[j].size();});
  for (int i = 0; i < qSegOrd.size(); i++) query[i] = queryOrig[qSegOrd[i]];

  // the distance from the centroid of each segment and the centroid of the
  // previous segments considered together
  centToCentDist.resize(query.size());
  CartesianPoint C(0, 0, 0);
  int N = 0;
  for (int L = 0; L < query.size(); L++) {
    CartesianPoint ci = query[L].getGeometricCenter();
    int n = query[L].size();
    C = (C*N + ci*n)/(N + n);
    centToCentDist[L].resize(n);
    for (int i = L + 1; i < query.size(); i++) {
      centToCentDist[L][i] = C.distance(query[i].getGeometricCenter());
    }
    N += n;
  }

  // set query masks (subsets of query segments involved at each recursion level)
  queryMasks.clear(); queryMasks.resize(query.size());
  for (int i = 0; i < query.size(); i++) {
    for (int L = i; L < query.size(); L++) {
      for (int j = 0; j < query[i].size(); j++) {
        queryMasks[L].push_back(query[i][j]);
      }
    }
  }
  updateGrids = true;

  // set gap constraints structure
  opts.resetGapConstraints(query.size());
}

AtomPointerVector FASST::getQuerySearchedAtoms() const {
  int len = 0, k = 0;
  for (int i = 0; i < queryOrig.size(); i++) len += queryOrig[i].size();
  AtomPointerVector atoms(len, NULL);
  for (int i = 0; i < queryOrig.size(); i++) {
    for (int j = 0; j < queryOrig[i].size(); j++) {
      atoms[k] = queryOrig[i][j]; k++;
    }
  }
  return atoms;
}

void FASST::addTarget(const string& pdbFile, short memSave) {
  Structure* targetStruct = new Structure(pdbFile, "QUIET");
  targetSource.push_back(targetInfo(pdbFile, targetFileType::PDB, 0, memSave));
  addTargetStructure(targetStruct, memSave);
}

void FASST::addTarget(const Structure& T, short memSave) {
  Structure* targetStruct = new Structure(T);
  targetSource.push_back(targetInfo(T.getName(), targetFileType::STRUCTURE, 0, memSave));
  addTargetStructure(targetStruct, memSave);
}

void FASST::addTargetStructure(Structure* targetStruct, short memSave) {
  if (memSave == 1) stripSidechains(*targetStruct);
  targetStructs.push_back(targetStruct);
  targets.push_back(AtomPointerVector());
  targSeqs.push_back(Sequence());
  AtomPointerVector& target = targets.back();
  Sequence& seq = targSeqs.back();
  // we don't care about the chain topology of the target, so append all residues
  for (int i = 0; i < targetStruct->chainSize(); i++) {
    bool foundAll = parseChain(targetStruct->getChain(i), &target, &seq);
    if ((memSave == 2) && !foundAll) MstUtils::error("for targets added under strict memory savings, all residues must be searchable", "FASST::addTargetStructure(Structure*, short)");
  }
  MstUtils::assertCond(target.size() > 0, "empty target named '" + targetStruct->getName() + "'", "FASST::addTargetStructure");
  seq.setName(targetStruct->getName());

  // orient the target structure in common frame (remember the transform)
  tr.push_back(TransformFactory::translate(-target.getGeometricCenter()));
  tr.back().apply(*targetStruct);

  // update extent for when will be creating proximity search objects
  mstreal _xlo, _ylo, _zlo, _xhi, _yhi, _zhi;
  ProximitySearch::calculateExtent(*targetStruct, _xlo, _ylo, _zlo, _xhi, _yhi, _zhi);
  if (targetStructs.size() == 1) {
    xlo = _xlo; ylo = _ylo; zlo = _zlo;
    xhi = _xhi; yhi = _yhi; zhi = _zhi;
  } else {
    xlo = min(xlo, _xlo); ylo = min(ylo, _ylo); zlo = min(zlo, _zlo);
    xhi = max(xhi, _xhi); yhi = max(yhi, _yhi); zhi = max(zhi, _zhi);
  }
  updateGrids = true;

  // chain lengths
  tightvector<int> chainLens(targetStruct->chainSize(), 0);
  for (int i = 0; i < chainLens.size(); i++) chainLens[i] = targetStruct->getChain(i).residueSize();
  targetChainLen.push_back(chainLens);

  // make stripped copies of atoms and destroy original target
  if (memSave == 2) {
    target = target.clone();
    for (int i = 0; i < target.size(); i++) target[i]->stripInfo();
    targetStructs.back() = NULL;
    delete targetStruct;
  }
}

void FASST::addTargets(const vector<string>& pdbFiles, short memSave) {
  for (int i = 0; i < pdbFiles.size(); i++) addTarget(pdbFiles[i], memSave);
}

bool FASST::parseChain(const Chain& C, AtomPointerVector* searchable, Sequence* seq) {
  bool foundAll = true;
  for (int i = 0; i < C.residueSize(); i++) {
    Residue& res = C.getResidue(i);
    AtomPointerVector bb;
    for (int k = 0; k < searchableAtomTypes.size(); k++) {
      Atom* a = NULL;
      for (int kk = 0; kk < searchableAtomTypes[k].size(); kk++) {
        if ((a = res.findAtom(searchableAtomTypes[k][kk], false)) != NULL) break;
      }
      if (a == NULL) { foundAll = false; break; }
      bb.push_back(a);
    }
    if (bb.size() == searchableAtomTypes.size()) {
      if (searchable != NULL) searchable->insert(searchable->end(), bb.begin(), bb.end());
      else {
        res.replaceAtoms(vector<Atom*>(), MstUtils::setdiff(res.getAtoms(), bb));
      }
      if (seq != NULL) seq->appendResidue(res.getName());
    }
  }
  return foundAll;
}

void FASST::addResidueStringProperties(int ti, const string& propType, const vector<string>& propVals) {
  if ((ti < 0) || (ti >= targetStructs.size())) MstUtils::error("requested target out of range: " + MstUtils::toString(ti), "FASST::addResidueStringProperties");
  int N = targetStructs[ti]->residueSize();
  if (N != propVals.size()) MstUtils::error("size of properties vector (" + MstUtils::toString(propVals.size()) + ") inconsistent with number of residues ("+ MstUtils::toString(N) +") for target: " + MstUtils::toString(ti), "FASST::addResidueStringProperties");
  resStringProperties[propType][ti] = propVals;
}

void FASST::addResidueProperties(int ti, const string& propType, const vector<mstreal>& propVals) {
  if ((ti < 0) || (ti >= targetStructs.size())) MstUtils::error("requested target out of range: " + MstUtils::toString(ti), "FASST::addResidueProperties");
  int N = targetStructs[ti]->residueSize();
  if (N != propVals.size()) MstUtils::error("size of properties vector inconsistent with number of residues for target: " + MstUtils::toString(ti), "FASST::addResidueProperties");
  resProperties[propType][ti] = propVals;
}

void FASST::addResiduePairProperties(int ti, const string& propType, const map<int, map<int, mstreal> >& propVals) {
  resPairProperties[propType][ti] = propVals;
}

void FASST::addResidueRelationship(int ti, const string& propType, int ri, int tj, int rj) {
  resRelProperties[propType][resAddress(ti, ri)].push_back(resAddress(tj, rj));
}

map<int, vector<FASST::resAddress>> FASST::getResidueRelationships(int ti, const string& propType) {
  simpleMap<resAddress, tightvector<resAddress>>& relMap = resRelProperties[propType];
  map<int, vector<resAddress>> ret;
  int beg = relMap.getLowerBound(resAddress(ti, 0));
  for (int i = beg; i < relMap.size(); i++) {
    resAddress ri = relMap.key(i);
    if (ri.targIndex() != ti) break;
    ret[ri.resIndex()] = relMap.value(i);
  }
  return ret;
}

bool FASST::hasResidueProperty(int ti, const string& propType, int ri) {
  return ((resProperties.find(propType) != resProperties.end()) &&
          (resProperties[propType].find(ti) != resProperties[propType].end()) &&
          (resProperties[propType][ti].size() > ri) && (ri >= 0));
}

bool FASST::hasResidueStringProperty(int ti, const string& propType, int ri) {
  return ((resStringProperties.find(propType) != resStringProperties.end()) &&
          (resStringProperties[propType].find(ti) != resStringProperties[propType].end()) &&
          (resStringProperties[propType][ti].size() > ri) && (ri >= 0));
}

mstreal FASST::getResidueProperty(int ti, const string& propType, int ri) {
  return hasResidueProperty(ti, propType, ri) ? resProperties[propType][ti][ri] : 0.0;
}

string FASST::getResidueStringProperty(int ti, const string& propType, int ri) {
  return hasResidueStringProperty(ti, propType, ri) ? resStringProperties[propType][ti][ri] : "";
}

bool FASST::hasResiduePairProperties(int ti, const string& propType, int ri) {
  return ((resPairProperties.find(propType) != resPairProperties.end()) &&
          (resPairProperties[propType].find(ti) != resPairProperties[propType].end()) &&
          (resPairProperties[propType][ti].find(ri) != resPairProperties[propType][ti].end()));
}

mstreal FASST::isResiduePairPropertyPopulated(const string& propType) {
  return (resPairProperties.find(propType) != resPairProperties.end());
}

mstreal FASST::isResidueRelationshipPopulated(const string& propType) {
  return (resRelProperties.find(propType) != resRelProperties.end());
}

map<int, mstreal> FASST::getResiduePairProperties(int ti, const string& propType, int ri) {
  return hasResiduePairProperties(ti, propType, ri) ? resPairProperties[propType][ti][ri] : map<int, mstreal>();
}

void FASST::writeDatabase(const string& dbFile) {
  fstream ofs; MstUtils::openFile(ofs, dbFile, fstream::out | fstream::binary, "FASST::writeDatabase");
  MstUtils::writeBin(ofs, 'V'); MstUtils::writeBin(ofs, (int) 1); // format version
  for (int ti = 0; ti < targetStructs.size(); ti++) {
    if (targetStructs[ti] == NULL) MstUtils::error("cannot write a database, in which full structures are not populated", "FASST::writeDatabase");
    MstUtils::writeBin(ofs, 'S'); // marks the start of a structure section
    targetStructs[ti]->writeData(ofs);
    for (auto p = resProperties.begin(); p != resProperties.end(); ++p) {
      if ((p->second).find(ti) != (p->second).end()) {
        vector<mstreal>& vals = (p->second)[ti];
        MstUtils::writeBin(ofs, 'P'); // marks the start of a residue property section
        MstUtils::writeBin(ofs, (string) p->first);
        MstUtils::assertCond(targetStructs[ti]->residueSize() == vals.size(), "the number of residue properties and residues does not agree for database entry", "FASST::writeDatabase(const string&)");
        for (int ri = 0; ri < vals.size(); ri++) MstUtils::writeBin(ofs, vals[ri]);
      }
    }
    for (auto p = resStringProperties.begin(); p != resStringProperties.end(); ++p) {
      if ((p->second).find(ti) != (p->second).end()) {
        vector<string>& vals = (p->second)[ti];
        MstUtils::writeBin(ofs, 'N'); // marks the start of a residue string property section
        MstUtils::writeBin(ofs, (string) p->first);
        MstUtils::assertCond(targetStructs[ti]->residueSize() == vals.size(), "the number of residue string properties and residues does not agree for database entry", "FASST::writeDatabase(const string&)");
        for (int ri = 0; ri < vals.size(); ri++) MstUtils::writeBin(ofs, vals[ri]);
      }
    }
    for (auto p = resPairProperties.begin(); p != resPairProperties.end(); ++p) {
      if ((p->second).find(ti) != (p->second).end()) {
        map<int, map<int, mstreal> >& vals = (p->second)[ti];
        MstUtils::writeBin(ofs, 'I'); // marks the start of a residue pair interaction property section
        MstUtils::writeBin(ofs, (string) p->first);
        MstUtils::writeBin(ofs, (int) vals.size());
        for (auto i = vals.begin(); i != vals.end(); ++i) {
          MstUtils::writeBin(ofs, (int) i->first);
          MstUtils::writeBin(ofs, (int) (i->second).size());
          for (auto j = (i->second).begin(); j != (i->second).end(); ++j) {
            MstUtils::writeBin(ofs, (int) j->first);
            MstUtils::writeBin(ofs, (mstreal) j->second);
          }
        }
      }
    }
  }

  // in the new version, we write all pair relation properties at the end
  for (auto p = resRelProperties.begin(); p != resRelProperties.end(); ++p) {
    simpleMap<resAddress, tightvector<resAddress>>& resRelProperty = p->second;
    MstUtils::writeBin(ofs, 'R'); // marks the start of a residue relation property section
    MstUtils::writeBin(ofs, (string) p->first);
    MstUtils::writeBin(ofs, (int) resRelProperty.size());
    for (int i = 0; i < resRelProperty.size(); i++) {
      MstUtils::writeBin(ofs, resRelProperty.key(i).targIndex());  // ti
      MstUtils::writeBin(ofs, resRelProperty.key(i).resIndex()); // ri
      tightvector<resAddress>& relatedList = resRelProperty.value(i);
      MstUtils::writeBin(ofs, (int) relatedList.size());
      for (int j = 0; j < relatedList.size(); j++) {
        MstUtils::writeBin(ofs, relatedList[j].targIndex());  // tj
        MstUtils::writeBin(ofs, relatedList[j].resIndex()); // rj
      }
    }
  }
  ofs.close();
}

void FASST::readDatabase(const string& dbFile, short memSave) {
  fstream ifs; MstUtils::openFile(ifs, dbFile, fstream::in | fstream::binary, "FASST::readDatabase");
  char sect; string name; mstreal val; string sval;
  int ver = 0;
  int ti = numTargets();
  MstUtils::readBin(ifs, sect);
  if (sect == 'V') {
    MstUtils::readBin(ifs, ver);
    MstUtils::readBin(ifs, sect);
  }
  if (sect != 'S') MstUtils::error("first section must be a structure one, while reading database file " + dbFile, "FASST::readDatabase(const string&)");
  while (ifs.peek() != EOF) {
    Structure* targetStruct = new Structure();
    streampos loc = ifs.tellg();
    targetStruct->readData(ifs);
    targetSource.push_back(targetInfo(dbFile, targetFileType::BINDATABASE, loc, memSave));
    int L = targetStruct->residueSize();
    addTargetStructure(targetStruct, memSave);
    while (ifs.peek() != EOF) {
      MstUtils::readBin(ifs, sect);
      if (sect == 'P') {
        MstUtils::readBin(ifs, name);
        vector<mstreal>& vals = resProperties[name][ti];
        vals.resize(L, 0);
        for (int i = 0; i < L; i++) {
          MstUtils::readBin(ifs, val);
          vals[i] = val;
        }
      } else if (sect == 'N') {
        MstUtils::readBin(ifs, name);
        vector<string>& vals = resStringProperties[name][ti];
        vals.resize(L);
        for (int i = 0; i < L; i++) {
          MstUtils::readBin(ifs, sval);
          vals[i] = sval;
        }
      } else if (sect == 'I') {
        MstUtils::readBin(ifs, name);
        map<int, map<int, mstreal> >& vals = resPairProperties[name][ti];
        int ri, rj, N, n; mstreal cd;
        MstUtils::readBin(ifs, N);
        for (int i = 0; i < N; i++) {
          MstUtils::readBin(ifs, ri);
          MstUtils::readBin(ifs, n);
          for (int j = 0; j < n; j++) {
            MstUtils::readBin(ifs, rj);
            MstUtils::readBin(ifs, cd);
            vals[ri][rj] = cd;
          }
        }
      } else if (sect == 'R') {
        MstUtils::readBin(ifs, name);
        simpleMap<resAddress, tightvector<resAddress>>& resRelProperty = resRelProperties[name];
        switch (ver) {
          case 0: {
            int ri, tj, rj, N, n1, n2;
            MstUtils::readBin(ifs, N);
            for (int i = 0; i < N; i++) {
              MstUtils::readBin(ifs, tj);
              MstUtils::readBin(ifs, n1);
              for (int j = 0; j < n1; j++) {
                MstUtils::readBin(ifs, ri);
                MstUtils::readBin(ifs, n2);
                tightvector<resAddress>& relatedList = resRelProperty[resAddress(ti, ri)];
                int off = relatedList.size();
                relatedList.resize(off + n2);
                for (int k = 0; k < n2; k++) {
                  MstUtils::readBin(ifs, rj);
                  relatedList[off + k] = resAddress(tj, rj);
                }
              }
            }
            break;
          }
          case 1: {
            // in the new version, we read them all at once
            resAddress ri, rj; int N, n;
            MstUtils::readBin(ifs, N);
            for (int i = 0; i < N; i++) {
              MstUtils::readBin(ifs, ri.targIndex());
              MstUtils::readBin(ifs, ri.resIndex());
              tightvector<resAddress>& relatedList = resRelProperty[ri];
              MstUtils::readBin(ifs, n);
              int off = relatedList.size();
              relatedList.resize(off + n);
              for (int j = 0; j < n; j++) {
                MstUtils::readBin(ifs, rj.targIndex());
                MstUtils::readBin(ifs, rj.resIndex());
                relatedList[off + j] = rj;
              }
            }
            break;
          }
          default:
            MstUtils::error("unknown database version " + MstUtils::toString(ver));
        }
      } else if (sect == 'S') {
        break;
      } else {
        MstUtils::error("unknown section type" + MstUtils::toString(sect) + ", while reading database file " + dbFile, "FASST::readDatabase(const string&)");
      }
    }
    ti++;
  }
  ifs.close();
}

void FASST::setSearchType(searchType _searchType) {
  type = _searchType;
  switch(type) {
    case searchType::CA:
      searchableAtomTypes =  {{"CA"}};
      break;
    case searchType::FULLBB:
      searchableAtomTypes =  {{"N", "NT"}, {"CA"}, {"C"}, {"O", "OT1", "OT2", "OXT"}};
      break;
    default:
      MstUtils::error("uknown search type '" + MstUtils::toString(type) + "' specified", "FASST::setSearchType");
  }
  atomsPerRes = searchableAtomTypes.size();
}

void FASST::stripSidechains(Structure& S) {
  for (int i = 0; i < S.chainSize(); i++) {
    parseChain(S[i]);
  }
}

void FASST::rebuildProximityGrids() {
  if (xlo == xhi) { xlo -= gridSpacing/2; xhi += gridSpacing/2; }
  if (ylo == yhi) { ylo -= gridSpacing/2; yhi += gridSpacing/2; }
  if (zlo == zhi) { zlo -= gridSpacing/2; zhi += gridSpacing/2; }
  int N = int(ceil(max(max((xhi - xlo), (yhi - ylo)), (zhi - zlo))/gridSpacing));
  for (int i = 0; i < ps.size(); i++) delete(ps[i]);
  ps.clear(); ps.resize(query.size(), NULL);
  for (int i = 0; i < query.size(); i++) {
    ps[i] = new ProximitySearch(xlo, ylo, zlo, xhi, yhi, zhi, N);
  }
  updateGrids = false;
}

void FASST::prepForSearch(int ti) {
  recLevel = 0;
  AtomPointerVector& target = targets[ti];
  if ((query.size() == 0) || (target.size() == 0)) {
    MstUtils::error("query and target must be set before starting search", "FASST::prepForSearch");
  }

  // if the search database has changed since the last search, update
  if (updateGrids) rebuildProximityGrids();

  // align every segment onto every admissible location on the target
  mstreal xc, yc, zc;
  segmentResiduals.resize(query.size());
  vector<vector<bool> > okAlignments(query.size());
  for (int i = 0; i < query.size(); i++) {
    bool seqConst = options().sequenceConstraintsSet() && options().getSequenceConstraints()->isSegmentConstrained(qSegOrd[i]);
    ps[i]->dropAllPoints();
    AtomPointerVector& seg = query[i];
    int Na = atomToResIdx(target.size()) - atomToResIdx(seg.size()) + 1; // number of possible alignments
    // make the default bad, so alignments skipped due to sequence constraints
    // get sorted to the bottom of the options list before they are removed
    segmentResiduals[i].resize(MstUtils::max(Na, 0), 9999.0);
    if (seqConst) {
      okAlignments[i].resize(segmentResiduals[i].size());
      options().getSequenceConstraints()->evalConstraint(qSegOrd[i], targSeqs[ti], okAlignments[i]);
    }
    AtomPointerVector targSeg(query[i].size(), NULL);
    for (int j = 0; j < Na; j++) {
      if (seqConst && !okAlignments[i][j]) { // save on calculating RMSDs for disallowed segment alignments
        if (query.size() > 1) ps[i]->addPoint(0, 0, 0, j); // add a dummy point, so point indexing is preserved
        continue;
      }
      // NOTE: can save on this in several ways:
      // 1. the centroid calculation is effectively already done inside RMSDCalculator::bestRMSD
      // 2. updating just one atom involves a simple centroid adjustment, rather than recalculation
      // 3. is there a speedup to be gained from re-calculating RMSD with one atom updated only?
      int off = resToAtomIdx(j);
      for (int k = 0; k < query[i].size(); k++) targSeg[k] = target[off + k];
      // AtomPointerVector targSeg = target.subvector(resToAtomIdx(j), resToAtomIdx(j) + query[i].size());
      segmentResiduals[i][j] = RC.bestResidual(query[i], targSeg);
      if (query.size() > 1) {
        targSeg.getGeometricCenter(xc, yc, zc);
        ps[i]->addPoint(xc, yc, zc, j);
      }
    }
  }

  // initialize remOptions; all options are available at top level
  remOptions.clear(); remOptions.resize(query.size());
  for (int L = 0; L < query.size(); L++) {
    remOptions[L].resize(query.size());
    for (int i = L; i < query.size(); i++) {
      // remOptions[L][i], where i < L, will hold an empty list of options, but
      // it won't be used and it is more convenient to access without thinking
      // of offsets (a tiny bit of memory waste for convenience)
      remOptions[L][i].setOptions(segmentResiduals[i], L == 0);
      if (!okAlignments[i].empty()) remOptions[L][i].intersectOptions(okAlignments[i]);
    }
  }

  // current residual starts with 0, since nothing is aligned yet
  currResidual = 0;
  currResiduals.resize(query.size(), 0);
  currRemBound = boundOnRemainder(true);

  // make room for target atom masks at different recursion levels
  if (targetMasks.size()) targetMasks.back().deletePointers();
  targetMasks.clear(); targetMasks.resize(query.size());
  for (int i = 0; i < query.size(); i++) {
    for (int L = i; L < query.size(); L++) {
      for (int j = 0; j < query[i].size(); j++) {
        if (L > i) {
          // for any repeated atoms on this level, make sure to point to the
          // same atoms as the mask on the previous level
          targetMasks[L].push_back(targetMasks[L-1][targetMasks[L-1].size() - query[i].size() + j]);
        } else {
          // and only create new atoms for the last segment
          targetMasks[L].push_back(new Atom(*(query[i][j])));
        }
      }
    }
  }

  // room for centroids of aligned sub-structure, at each recursion level
  currCents.clear();
  currCents.resize(query.size(), CartesianPoint(0, 0, 0));

  // mark chain beginning and end indices (if gap constraints are present)
  if (opts.gapConstraintsExist() || (opts.getRedundancyCut() != 1)) {
    fillTargetChainInfo(ti);
  } else {
    targChainBeg.resize(0); targChainEnd.resize(0);
  }
}

void FASST::fillTargetChainInfo(int ti) {
  AtomPointerVector& target = targets[ti];
  targChainBeg.resize(atomToResIdx(target.size()));
  targChainEnd.resize(atomToResIdx(target.size()));

  tightvector<int>& chainLengths = targetChainLen[ti];
  int ri = 0, cb = 0;
  for (int i = 0; i < chainLengths.size(); i++) {
    for (int j = 0; j < chainLengths[i]; j++, ri++) {
      targChainBeg[ri] = cb;
      targChainEnd[ri] = cb + chainLengths[i] - 1;
    }
    cb += chainLengths[i];
  }
}

mstreal FASST::boundOnRemainder(bool compute) {
  if (compute) {
    currRemBound = 0;
    for (int i = recLevel + 1; i < query.size(); i++) currRemBound += remOptions[recLevel][i].bestCost();
  }
  return currRemBound;
}

mstreal FASST::centToCentTol(int i) {
  mstreal remRes = residualCut - currResidual - currRemBound;
  if (remRes < 0) return -1.0;
  return sqrt((remRes * (queryMasks[recLevel].size() + query[i].size())) / (queryMasks[recLevel].size() * query[i].size()));
}

fasstSolutionSet FASST::search() {
  // auto begin = chrono::high_resolution_clock::now();
  // int prepTime = 0;
  int numSegs = query.size();
  opts.validateSearchRequest(numSegs);
  if (opts.isMinNumMatchesSet()) setCurrentRMSDCutoff(INFINITY);
  else setCurrentRMSDCutoff(opts.getRMSDCutoff());
  solutions.init(numSegs);
  bool redSet = opts.isRedundancyCutSet() || opts.isRedundancyPropertySet();
  bool doRedBar = redSet && (numSegs > 1); // should we apply special "barrier" RMSD cutoffs to partial matches that are
                                           // already known to be redundant to something in the current list of solutions?
  vector<int> segLen(numSegs); // number of residues in each query segment
  for (int i = 0; i < numSegs; i++) segLen[i] = atomToResIdx(query[i].size());
  vector<mstreal> ccTol(numSegs, -1.0);
  for (currentTarget = 0; currentTarget < targets.size(); currentTarget++) {
    // auto beginPrep = chrono::high_resolution_clock::now();
    if (doRedBar) {
      resetCurrentRMSDCutoff(); // if it was previously temporarily set
      solutions.resetAlignRedBarrierData(targets[currentTarget].size());
    }
    prepForSearch(currentTarget);
    // auto endPrep = chrono::high_resolution_clock::now();
    // prepTime += chrono::duration_cast<std::chrono::microseconds>(endPrep-beginPrep).count();
    vector<int> okLocations, badLocations;
    okLocations.reserve(targets[currentTarget].size()); badLocations.reserve(targets[currentTarget].size());
    while (true) {
      // Have to do three things:
      // 1. pick the best choice (from available ones) for the current segment,
      // remove it from the list of options, and move onto the next recursion level
      if (remOptions[recLevel][recLevel].empty()) {
        currAlignment[recLevel] = -1;
        if (recLevel > 0) {
          recLevel--;
          continue;
        } else {
          break; // search exhausted
        }
      }
      currAlignment[recLevel] = remOptions[recLevel][recLevel].bestChoice();
      remOptions[recLevel][recLevel].removeOption(currAlignment[recLevel]);

      // if redundancy removal is set, and there are matches in the current list
      // of solutions that are redundant with the current partial solution, any
      // full realization of the current partial solution will only be accepted
      // if they it improves upon the best RMSD of any of these redundant solu-
      // tions. So, temporarily lower the current RMSD threshold, if applicable,
      // but keep track of which segment in the current alignment this was due
      // to (which segment had the redundancy), so that this chane can be unwound
      // when this segment's alignment changes in the partial solution.
      if (doRedBar) {
        if (rmsdPriority() >= recLevel) resetCurrentRMSDCutoff(recLevel - 1);
        mstreal barrier = solutions.alignRedBarrier(qSegOrd[recLevel], currAlignment[recLevel]);
        if (getCurrentRMSDCutoff() > barrier) setCurrentRMSDCutoff(barrier, recLevel);
      }

      // 2. compute the total residual from the current alignment
      mstreal curBound = currentAlignmentResidual(true) + boundOnRemainder(true);
      if (curBound > residualCut) continue;
      // if (query.size() > 1) updateQueryCentroids();

      // 3. update update remaining options for subsequent segments based on the
      // newly made choice. The set of options on the next recursion level is a
      // subset of the set of options on the previous level.
      int remSegs = numSegs - (recLevel + 1);
      if (remSegs > 0) {
        bool levelExhausted = false;
        int nextLevel = recLevel + 1;
        // copy remaining options from the previous recursion level. This way,
        // we can compute bounds on this level and can do set intersections to
        // further narrow this down
        for (int i = nextLevel; i < numSegs; i++) {
          remOptions[nextLevel][i].copyIn(remOptions[nextLevel-1][i]);
          // except that segments cannot overlap, so remove from consideration
          // all alignments that overlap with the segments that was just placed
          remOptions[nextLevel][i].removeOptions(currAlignment[recLevel] - segLen[i] + 1,
                                                 currAlignment[recLevel] + segLen[recLevel] - 1);
        }
        // if any gap constraints exist, limit options at this recursion level accordingly
        if (opts.gapConstraintsExist()) {
          for (int j = 0; j < nextLevel; j++) {
            for (int i = nextLevel; i < numSegs; i++) {
              remOptions[nextLevel][i].constrainRange(targChainBeg[currAlignment[j]], targChainEnd[currAlignment[j]]);
              if (opts.minGapConstrained(qSegOrd[i], qSegOrd[j]))
                remOptions[nextLevel][i].constrainLE(currAlignment[j] - opts.getMinGap(qSegOrd[i], qSegOrd[j]) - segLen[i]);
              if (opts.maxGapConstrained(qSegOrd[i], qSegOrd[j]))
                remOptions[nextLevel][i].constrainGE(currAlignment[j] - opts.getMaxGap(qSegOrd[i], qSegOrd[j]) - segLen[i]);
              if (opts.minGapConstrained(qSegOrd[j], qSegOrd[i]))
                remOptions[nextLevel][i].constrainGE(currAlignment[j] + opts.getMinGap(qSegOrd[j], qSegOrd[i]) + segLen[j]);
              if (opts.maxGapConstrained(qSegOrd[j], qSegOrd[i]))
                remOptions[nextLevel][i].constrainLE(currAlignment[j] + opts.getMaxGap(qSegOrd[j], qSegOrd[i]) + segLen[j]);
              // if (minGapSet[qSegOrd[i]][qSegOrd[j]]) remOptions[nextLevel][i].constrainLE(currAlignment[j] - minGap[qSegOrd[i]][qSegOrd[j]] - segLen[i]);
              // if (maxGapSet[qSegOrd[i]][qSegOrd[j]]) remOptions[nextLevel][i].constrainGE(currAlignment[j] - maxGap[qSegOrd[i]][qSegOrd[j]] - segLen[i]);
              // if (minGapSet[qSegOrd[j]][qSegOrd[i]]) remOptions[nextLevel][i].constrainGE(currAlignment[j] + minGap[qSegOrd[j]][qSegOrd[i]] + segLen[j]);
              // if (maxGapSet[qSegOrd[j]][qSegOrd[i]]) remOptions[nextLevel][i].constrainLE(currAlignment[j] + maxGap[qSegOrd[j]][qSegOrd[i]] + segLen[j]);
              if (remOptions[nextLevel][i].empty()) { levelExhausted = true; break; }
            }
            if (levelExhausted) break;
          }
          if (levelExhausted) continue;
        }
        mstreal di, de, d, dePrev, eps = 10E-8;
        CartesianPoint& currCent = currCents[recLevel];
        for (int c = 0; true; c++) {
          bool updated = false;
          for (int i = nextLevel; i < numSegs; i++) {
            FASST::optList& remSet = remOptions[nextLevel][i];
            de = centToCentTol(i);
            if (de < 0) { levelExhausted = true; break; }
            di = centToCentDist[recLevel][i];
            dePrev = ((c == 0) ? -1 : ccTol[i]);
            int numLocs = remSet.size();

            // If the set of options for the current segment was arrived at,
            // in part, by limiting center-to-center distances, then we will
            // tighten that list by removing options that are outside of the
            // range allowed at this recursion level. Otherwise, we will do a
            // general proximity search given the current tolerance and will
            // tighten the list that way.
            if (dePrev < 0) {
              okLocations.resize(0);
              ps[i]->pointsWithin(currCent, max(di - de, 0.0), di + de, &okLocations);
              remSet.intersectOptions(okLocations);
            } else if (dePrev - de > eps) {
              badLocations.resize(0);
              ps[i]->pointsWithin(currCent, max(di - dePrev, 0.0), di - de - eps, &badLocations);
              ps[i]->pointsWithin(currCent, di + de + eps, di + dePrev, &badLocations);
              for (int k = 0; k < badLocations.size(); k++) {
                remSet.removeOption(badLocations[k]);
              }
            }
            ccTol[i] = de;
            if (numLocs != remSet.size()) {
              // this both updates the bound and checks that there are still
              // feasible solutions left
              if ((remSet.empty()) || (currResidual + boundOnRemainder(true) > residualCut)) {
                levelExhausted = true; break;
              }
              updated = true;
            }
          }
          if (levelExhausted) break;
          if (!updated) break;
        }
        if (levelExhausted) continue;
        recLevel = nextLevel;
      } else {
        // if at the lowest recursion level already, then record the solution
        fasstSolution sol(currAlignment, sqrt(currResidual/querySize), currentTarget, currentTransform(), segLen, qSegOrd);
        bool inserted = false;
        if (opts.isRedundancyCutSet()) {
          addSequenceContext(sol);
          inserted = solutions.insert(sol, opts.getRedundancyCut());
        } else if (opts.isRedundancyPropertySet()) {
          inserted = solutions.insert(sol, getRedundancyPropertyMap());
        } else {
          inserted = solutions.insert(sol);
        }

        if (doRedBar && !inserted) {
          for (int rL = 0; rL < numSegs; rL++) {
            mstreal barrier = solutions.alignRedBarrier(qSegOrd[rL], currAlignment[rL]);
            if (getCurrentRMSDCutoff() > barrier) setCurrentRMSDCutoff(barrier, rL);
          }
        }

        if (opts.isSufficientNumMatchesSet() && (solutions.size() == opts.getSufficientNumMatches())) return solutions;
        if (opts.isMaxNumMatchesSet() && (solutions.size() > opts.getMaxNumMatches())) {
          solutions.erase(--solutions.end());
          setCurrentRMSDCutoff(solutions.worstRMSD());
        } else if (opts.isMinNumMatchesSet() && (solutions.size() > opts.getMinNumMatches()) && (rmsdCut > opts.getRMSDCutoff())) {
          if (solutions.worstRMSD() > opts.getRMSDCutoff()) solutions.erase(--solutions.end());
          setCurrentRMSDCutoff(MstUtils::max(solutions.worstRMSD(), opts.getRMSDCutoff()));
        }
      }
    }
  }
  // auto end = chrono::high_resolution_clock::now();
  // int searchTime = chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
  // prepTime = prepTime/1000;
  // cout << "prep time " << prepTime << " ms" << std::endl;
  // cout << "total time " << searchTime << " ms" << std::endl;
  // cout << "prep time was " << (100.0*prepTime/searchTime) << " % of the total" << std::endl;

  solutions.clearTempData();
  return solutions;
}

mstreal FASST::currentAlignmentResidual(bool compute, bool setTransform) {
  if (compute) {
    if ((query.size() == 1) && !setTransform) {
      // this is a special case, because will not need to calculate centroid
      // locations for subsequent sub-queries
      currResidual = segmentResiduals[0][currAlignment[0]];
    } else {
      // fill up sub-alignment with target atoms
      AtomPointerVector& queryMask = queryMasks[recLevel];
      AtomPointerVector& targetMask = targetMasks[recLevel];
      int N = targetMask.size();
      int n = query[recLevel].size();
      int dN = N - n;
      int currPos = currAlignment[recLevel];
      int si = resToAtomIdx(currPos);
      AtomPointerVector& target = targets[currentTarget];
      for (int i = 0; i < n; i++) {
        targetMask[dN + i]->setCoor(target[si + i]->getX(), target[si + i]->getY(), target[si + i]->getZ());
      }
      currResidual = RC.bestResidual(targetMask, queryMask, setTransform);
      currResiduals[recLevel] = currResidual;
      if (query.size() > 1) {
        if (recLevel == 0) {
          currCents[recLevel] = ps[recLevel]->getPoint(currPos);
        } else {
          /* The following is A LOT faster than the equivalent:
           * currCents[recLevel] = (currCents[recLevel - 1] * (N - n) +  ps[recLevel]->getPoint(currAlignment[recLevel]) * n) / N;
           * (the above causes a temporary creation of a CartesianPoint object) */
          CartesianPoint& prevC = currCents[recLevel - 1];
          CartesianPoint& currC = ps[recLevel]->getPoint(currPos);
          for (int i = 0; i < 3; i++) {
            currCents[recLevel][i] = (prevC[i] * dN +  currC[i] * n) / N;
          }
        }
      }
    }
  }
  return currResidual;
}

Transform FASST::currentTransform() {
  currentAlignmentResidual(true, true);
  return Transform(RC.lastRotation(), RC.lastTranslation());
}

void FASST::getMatchStructure(const fasstSolution& sol, Structure& match, bool detailed, matchType type, bool algn) {
  vector<Structure> matches;
  fasstSolutionSet solSet(sol);
  getMatchStructures(solSet, matches, detailed, type, algn);
  match = matches[0];
}

Structure FASST::getMatchStructure(const fasstSolution& sol, bool detailed, matchType type, bool algn) {
  Structure match; getMatchStructure(sol, match, detailed, type, algn); return match;
}

void FASST::getMatchStructures(fasstSolutionSet& sols, vector<Structure>& matches, bool detailed, matchType type, bool algn) {
  // hash solutions by the target they come from, to visit each target only once
  map<int, vector<int> > solsFromTarget;
  for (int i = 0; i < sols.size(); i++) {
    const fasstSolution& sol = sols[i];
    int idx = sol.getTargetIndex();
    if ((idx < 0) || (idx >= targets.size())) {
      MstUtils::error("supplied FASST solution is pointing to an out-of-range target", "FASST::getMatchStructures");
    }
    solsFromTarget[idx].push_back(i);
  }

  // visit each target
  Structure dummy; RMSDCalculator rc;
  AtomPointerVector matchAtoms;
  matches.resize(sols.size(), Structure());
  for (auto it = solsFromTarget.begin(); it != solsFromTarget.end(); ++it) {
    int idx = it->first;
    Structure* targetStruct = targetStructs[idx];
    AtomPointerVector& target = targets[idx];
    bool reread = (detailed && targetSource[idx].memSave) || (targetStruct == NULL); // should we re-read the target structure?
    Transform& transf = tr[idx];
    if (reread) {
      dummy.reset();
      // re-read structure
      if (targetSource[idx].type == targetFileType::PDB) {
        dummy.readPDB(targetSource[idx].file, "QUIET");
      } else if (targetSource[idx].type == targetFileType::BINDATABASE) {
        fstream ifs; MstUtils::openFile(ifs, targetSource[idx].file, fstream::in | fstream::binary, "FASST::getMatchStructures");
        ifs.seekg(targetSource[idx].loc);
        dummy.readData(ifs);
        ifs.close();
      } else if (targetSource[idx].type == targetFileType::STRUCTURE) {
        MstUtils::error("cannot produce a detailed match if target was initialized from object", "FASST::getMatchStructures");
      } else {
        MstUtils::error("don't know how to re-read target of this type", "FASST::getMatchStructures");
      }
      transf.apply(dummy);
      targetStruct = &dummy;
      if (!detailed) stripSidechains(*targetStruct);
    }

    // visit each solution from this target
    vector<int>& solIndices = it->second;
    for (int i = 0; i < solIndices.size(); i++) {
      int solIndex = solIndices[i];
      const fasstSolution& sol = sols[solIndex];
      Structure& match = matches[solIndex];
      match.setName(targetStruct->getName());

      // cut out the part of the target Structure that will constitute the returned match
      if (type == matchType::FULL) {
        match = *targetStruct;
      } else {
        vector<int> resIndices = getMatchResidueIndices(sol, type);
        for (auto ri = resIndices.begin(); ri != resIndices.end(); ri++) {
          if (targetStructs[idx] == NULL) {
            match.addResidue(&(targetStruct->getResidue((*ri))));
          } else {
            Residue* res = target[resToAtomIdx(*ri)]->getResidue();
            if (reread) res = &(targetStruct->getResidue(res->getResidueIndex()));
            match.addResidue(res);
          }
        }
      }

      // align matching region onto query, transforming the match itself
      if (algn) {
        sol.getTransform().apply(match);
      } else {
        transf.inverse().apply(match);
      }
    }
  }
}

vector<Sequence> FASST::getMatchSequences(fasstSolutionSet& sols, matchType type) {
  vector<Sequence> seqs(sols.size());
  for (int i = 0; i < sols.size(); i++) {
    const fasstSolution& sol = sols[i];
    int idx = sol.getTargetIndex();
    MstUtils::assertCond((idx >= 0) && (idx < targSeqs.size()), "supplied FASST solution is pointing to an out-of-range target", "FASST::getMatchSequences");
    vector<int> alignment = sol.getAlignment();
    seqs[i].setName(targSeqs[idx].getName());

    // isolate out the part of the target Sequence that will constitute the returned match
    vector<int> resIndices = getMatchResidueIndices(sol, type);
    for (auto ri = resIndices.begin(); ri != resIndices.end(); ri++) seqs[i].appendResidue(targSeqs[idx][*ri]);
  }
  return seqs;
}

Sequence FASST::getMatchSequence(const fasstSolution& sol, matchType type) {
  fasstSolutionSet solSet(sol);
  vector<Sequence> seqs = getMatchSequences(solSet, type);
  return seqs[0];
}

vector<mstreal> FASST::getResidueProperties(const fasstSolution& sol, const string& propType, matchType type) {
  fasstSolutionSet solSet; solSet.insert(sol);
  vector<vector<mstreal> > props = getResidueProperties(solSet, propType, type);
  return props[0];
}

vector<vector<mstreal> > FASST::getResidueProperties(fasstSolutionSet& sols, const string& propType, matchType type) {
  vector<vector<mstreal> > props(sols.size());
  for (int i = 0; i < sols.size(); i++) {
    const fasstSolution& sol = sols[i];
    int idx = sol.getTargetIndex();
    MstUtils::assertCond((idx >= 0) && (idx < targets.size()), "supplied FASST solution is pointing to an out-of-range target", "FASST::getMatchSequences");
    AtomPointerVector& target = targets[idx];
    if (!isResiduePropertyDefined(propType, idx)) {
      MstUtils::error("target with index " + MstUtils::toString(idx) + " does not have property type " + propType, "FASST::getResidueProperties(fasstSolutionSet&, const string&, matchType)");
    }
    vector<mstreal>& propVals = resProperties[propType][idx];
    vector<int> resIndices = getMatchResidueIndices(sol, type);
    props[i].resize(resIndices.size()); int ii = 0;
    for (auto ri = resIndices.begin(); ri != resIndices.end(); ri++, ii++) {
      // if we have the full structure, then we have the ability to differentiate
      // between the original structure and the part that is searched over (e.g.,
      // some residues could be skipped for any reason). if instead we are in a
      // mode where the original structure was not saved, that implies that all
      // residues in the original structure made it to the portion being searched
      // over (otherwise, discarding the original would have been caught as an error)
      props[i][ii] = (targetStructs[idx] != NULL) ? propVals[target[resToAtomIdx(*ri)]->getResidue()->getResidueIndex()] : propVals[*ri];
    }
  }
  return props;
}

bool FASST::isResiduePropertyDefined(const string& propType, int ti) {
  return ((resProperties.find(propType) != resProperties.end()) || (resProperties[propType].find(ti) != resProperties[propType].end()));
}

bool FASST::isResiduePropertyDefined(const string& propType) {
  return (resProperties.find(propType) != resProperties.end());
}

bool FASST::isResidueStringPropertyDefined(const string& propType, int ti) {
  return ((resStringProperties.find(propType) != resStringProperties.end()) || (resStringProperties[propType].find(ti) != resStringProperties[propType].end()));
}

bool FASST::isResidueStringPropertyDefined(const string& propType) {
  return (resStringProperties.find(propType) != resStringProperties.end());
}

vector<int> FASST::getMatchResidueIndices(const fasstSolution& sol, matchType type) {
  vector<int> residueIndices;
  vector<int> alignment = sol.getAlignment();
  switch(type) {
    case matchType::REGION: {
      for (int k = 0; k < alignment.size(); k++) {
        int si = alignment[k];
        int L = sol.segLength(k);
        for (int ri = si; ri < si + L; ri++) residueIndices.push_back(ri);
      }
      break;
    }
    case matchType::WITHGAPS: {
      set<int> toInclude;
      for (int i = 0; i < sol.numSegments(); i++) {
        // each individual segments should be included
        for (int k = alignment[i]; k < alignment[i] + sol.segLength(i); k++) toInclude.insert(k);
        // then some gaps also
        for (int j = 0; j < sol.numSegments(); j++) {
          if (i == j) continue;
          // if gap constrained, fill in between these two segments
          if (opts.gapConstrained(i, j)) {
            if (alignment[i] > alignment[j]) MstUtils::error("solution not consistent with current gap constraints", "FASST::getMatchResidueIndices");
            for (int k = alignment[i] + sol.segLength(i); k < alignment[j]; k++) {
              toInclude.insert(k);
            }
          }
        }
      }
      for (auto it = toInclude.begin(); it != toInclude.end(); it++) {
        residueIndices.push_back(*it);
      }
      break;
    }
    case matchType::FULL: {
      int idx = sol.getTargetIndex();
      MstUtils::assertCond((idx >= 0) && (idx < targetStructs.size()), "supplied FASST solution is pointing to an out-of-range target", "FASST::getMatchResidueIndices");
      for (int i = 0; i < getTargetResidueSize(idx); i++) residueIndices.push_back(i);
      break;
    }
    default:
      MstUtils::error("unknown match output type", "FASST::getMatchResidueIndices");
  }

  return residueIndices;
}

vector<mstreal> FASST::matchRMSDs(fasstSolutionSet& sols, const AtomPointerVector& query, bool update) {
  vector<mstreal> rmsds(sols.size(), 0);
  if (sols.size() == 0) return rmsds;
  AtomPointerVector match(query.size(), NULL);
  RMSDCalculator rc;
  for (int c = 0; c < sols.size(); c++) {
    fasstSolution& sol = sols[c];
    int idx = sol.getTargetIndex();
    AtomPointerVector& target = targets[idx];
    int k = 0;
    for (int i = 0; i < sol.numSegments(); i++) {
      int si = sol[i];
      int L = sol.segLength(i);
      if (k + resToAtomIdx(L) > match.size()) MstUtils::error("solution alignment size is larger than the specified query", "FASST::matchRMSDs(const fasstSolutionSet& sols, AtomPointerVector& query)");
      for (int ri = si; ri < si + L; ri++) {
        for (int ai = resToAtomIdx(ri); ai < resToAtomIdx(ri) + atomsPerRes; ai++) {
          match[k] = target[ai]; k++;
        }
      }
    }
    if (k != match.size()) MstUtils::error("solution alignment size is smaller than the specified query", "FASST::matchRMSDs(const fasstSolutionSet& sols, AtomPointerVector& query)");
    rmsds[c] = rc.bestRMSD(match, query, update);
    if (update) {
      sol.setTransform(Transform(rc.lastRotation(), rc.lastTranslation()));
      sol.setRMSD(rmsds[c]);
    }
  }
  return rmsds;
}


string FASST::toString(const fasstSolution& sol) {
  stringstream ss;
  int ti = sol.getTargetIndex();
  ss << std::setprecision(6) << std::fixed << sol.getRMSD() << " " << ((targetStructs[ti] == NULL) ? MstUtils::toString(ti) : targetStructs[ti]->getName()) << " [" << MstUtils::vecToString(sol.getAlignment(), ", ") << "]";
  return ss.str();
}

void FASST::addSequenceContext(fasstSolution& sol) {
  int currentTarget = sol.getTargetIndex();
  Sequence& targSeq = targSeqs[currentTarget];
  if (targChainEnd.empty()) FASST::fillTargetChainInfo(currentTarget);
  sol.addSequenceContext(targSeqs[currentTarget], opts.getContextLength(), targChainBeg, targChainEnd);
}

void FASST::addSequenceContext(fasstSolutionSet& sols) {
  map<int, vector<int> > solsFromTarget; int i = 0;
  for (auto it = sols.begin(); it != sols.end(); ++it, ++i) solsFromTarget[it->getTargetIndex()].push_back(i);
  for (auto it = solsFromTarget.begin(); it != solsFromTarget.end(); ++it) {
    FASST::fillTargetChainInfo(it->first);
    vector<int>& solInds = it->second;
    for (int i = 0; i < solInds.size(); i++) {
      addSequenceContext(sols[solInds[i]]);
    }
  }
}


/* --------- fasstSolution --------- */
fasstSolution::fasstSolution(const vector<int>& _alignment, mstreal _rmsd, int _target, const Transform& _tr, const vector<int>& _segLengths, vector<int> segOrder) {
  alignment = _alignment; rmsd = _rmsd; targetIndex = _target; segLengths = _segLengths; tr = _tr;
  if (!segOrder.empty()) {
    for (int i = 0; i < _alignment.size(); i++) {
      alignment[segOrder[i]] = _alignment[i];
      segLengths[segOrder[i]] = _segLengths[i];
    }
  }
  context = NULL;
}

fasstSolution::fasstSolution(const fasstSolution& _sol) {
  alignment = _sol.alignment; rmsd = _sol.rmsd; targetIndex = _sol.targetIndex;
  segLengths = _sol.segLengths; tr = _sol.tr;
  if (_sol.context != NULL) context = new solContext(*(_sol.context));
  else context = NULL;
}

fasstSolution::fasstSolution(const fasstSolutionAddress& addr, const vector<int> segLen) {
  alignment = addr.alignment;
  targetIndex = addr.targetIndex;
  segLengths = segLen;
  context = NULL;
  rmsd = 0.0;
}

void fasstSolution::setSeqContext(const vector<Sequence>& _segSeq, const vector<Sequence>& _nSeq, const vector<Sequence>& _cSeq) {
  if (context == NULL) context = new solContext();
  context->segSeq = _segSeq; context->nSeq = _nSeq; context->cSeq = _cSeq;
}

void fasstSolution::setStructContext(const vector<AtomPointerVector>& _segStr, const vector<AtomPointerVector>& _nStr, const vector<AtomPointerVector>& _cStr) {
  if (context == NULL) context = new solContext();
  context->segStr = _segStr; context->nStr = _nStr; context->cStr = _cStr;
}

void fasstSolution::write(ostream& _os) const {
  MstUtils::writeBin(_os, rmsd);
  MstUtils::writeBin(_os, targetIndex);
  MstUtils::writeBin(_os, (int) alignment.size());
  for (int i = 0; i < alignment.size(); i++) {
    MstUtils::writeBin(_os, alignment[i]);
    MstUtils::writeBin(_os, segLengths[i]);
  }
  tr.write(_os);
  if (context == NULL) MstUtils::writeBin(_os, false);
  else { MstUtils::writeBin(_os, true); context->write(_os); }
}

void fasstSolution::read(istream& _is) {
  MstUtils::readBin(_is, rmsd);
  MstUtils::readBin(_is, targetIndex);
  int len; MstUtils::readBin(_is, len);
  alignment.resize(len); segLengths.resize(len);
  for (int i = 0; i < alignment.size(); i++) {
    MstUtils::readBin(_is, alignment[i]);
    MstUtils::readBin(_is, segLengths[i]);
  }
  tr.read(_is);
  bool hasCont; MstUtils::readBin(_is, hasCont);
  if (context != NULL) { delete(context); context = NULL; }
  if (hasCont) {
    context = new solContext();
    context->read(_is);
  }
}

void fasstSolution::addSequenceContext(const Sequence& targSeq, int contLen, const vector<int>& targChainBeg, const vector<int>& targChainEnd) {
  vector<Sequence> segs(alignment.size()), ntPad(alignment.size()), ctPad(alignment.size());
  for (int i = 0; i < alignment.size(); i++) {
    int Li = segLength(i);
    segs[i].resize(Li);
    for (int ri = alignment[i]; ri < alignment[i] + Li; ri++) {
      segs[i][ri - alignment[i]] = targSeq[ri];
    }
    if (contLen > segs[i].size()) {
      // pad on both ends, to accommodate different scenarios
      int pad = contLen - segs[i].size();
      ntPad[i].resize(pad, SeqTools::gapIdx()); ctPad[i].resize(pad, SeqTools::gapIdx());
      for (int ri = 0; ri < pad; ri++) {
        if (alignment[i] + Li + ri > targChainEnd[alignment[i]]) break; // stop if reach chain C-terminus
        ctPad[i][ri] = targSeq[alignment[i] + Li + ri];
      }
      for (int ri = 0; ri < pad; ri++) {
        if (alignment[i] - ri - 1 < targChainBeg[alignment[i]]) break; // stop if reach chain N-terminus
        ntPad[i][ntPad[i].size() - ri - 1] = targSeq[alignment[i] - ri - 1];
      }
    }
  }
  setSeqContext(segs, ntPad, ctPad);
}

void fasstSolution::addSequenceContext(const Structure& target, int contLen) {
  vector<int> targChainBeg(target.residueSize()), targChainEnd(target.residueSize());
  int off = 0;
  for (int ci = 0; ci < target.chainSize(); ci++) {
    for (int i = 0; i < target[ci].residueSize(); i++) {
      targChainBeg[off + i] = off;
      targChainEnd[off + i] = off + target[ci].residueSize() - 1;
    }
    off += target[ci].residueSize();
  }
  addSequenceContext(Sequence(target), contLen, targChainBeg, targChainEnd);
}


/* --------- fasstSolutionAddress --------- */
void fasstSolutionAddress::write(ostream& _os) const {
  MstUtils::writeBin(_os, targetIndex);
  MstUtils::writeBin(_os, (int) alignment.size());
  for (int i = 0; i < alignment.size(); i++) MstUtils::writeBin(_os, alignment[i]);
}

void fasstSolutionAddress::read(istream& _is) {
  MstUtils::readBin(_is, targetIndex);
  int len; MstUtils::readBin(_is, len);
  alignment.resize(len);
  for (int i = 0; i < alignment.size(); i++) MstUtils::readBin(_is, alignment[i]);
}

/* --------- fasstSolutionSet --------- */
fasstSolutionSet::fasstSolutionSet(const fasstSolutionSet& sols) {
  *this = sols;
}

fasstSolutionSet::fasstSolutionSet(const fasstSolution& sol) {
  updated = false;
  insert(sol);
}

fasstSolutionSet::fasstSolutionSet(const vector<fasstSolutionAddress>& addresses, const vector<int>& segLengths) {
  updated = false;
  for (int i = 0; i < addresses.size(); i++) insert(fasstSolution(addresses[i], segLengths));
}

fasstSolutionSet& fasstSolutionSet::operator=(const fasstSolutionSet& sols) {
  updated = false;
  solsSet.clear(); solsVec.clear();
  // NOTE: I currently do not believe that solsByCenRes should be returned to the user,
  // but if I change my mind later, this will need to be uncommented; also see fasstSolutionSet::insert(const fasstSolution&, mstreal)
  // solsByCenRes.resize(sols.solsByCenRes.size());
  for (auto it = sols.begin(); it != sols.end(); ++it) this->insert(*it);
  return *this;
}

bool fasstSolutionSet::insert(const fasstSolution& sol, mstreal redundancyCut) {
  // apply a redudancy filter, if needed
  fasstSolution* toRemove = NULL;
  bool algnBarInfoSet = isAlignRedBarrierDataSet();
  if ((redundancyCut < 1) && sol.seqContextDefined()) {
    const vector<Sequence>& segSeqs = sol.segmentSeqs();
    const vector<Sequence>& nTermPad = sol.nTermContext();
    const vector<Sequence>& cTermPad = sol.cTermContext();
    // compare this solution to each previously accepted solution
    for (auto it = solsSet.begin(); it != solsSet.end(); ++it) {
      fasstSolution* psol = (fasstSolution*) &(*it);
      const vector<Sequence>& segSeqsPrev = psol->segmentSeqs();
      const vector<Sequence>& nTermPadPrev = psol->nTermContext();
      const vector<Sequence>& cTermPadPrev = psol->cTermContext();
      // compare the contexts of each segment:
      for (int i = 0; i < sol.numSegments(); i++) {
        int numID = 0, numTot = 0;
        /* The total length of the alignment we want to have, at least some L,
         * is the length the segment itself, Li, plus whatever extra context
         * (from either end), if needed, to get to at least L. Since both N- and
         * C-terminal paddings are of the same length, equal to max(0, L - Li),
         * by constructoin, we can deduce max(L, Li) as Li + max(0, L - Li). */
        int contextLength = segSeqs[i].size() + nTermPad[i].size();

        // first the segment itself
        for (int k = 0; k < segSeqs[i].size(); k++) {
          numTot++;
          if (segSeqs[i][k] == segSeqsPrev[i][k]) numID++;
        }

        // then alternate expanding in C- and N-terminal directions
        bool nEnd = false, cEnd = false;
        for (int k = 0; k < nTermPad[i].size(); k++) {
          cEnd = cEnd || (cTermPad[i][k] == SeqTools::gapIdx()) || (cTermPadPrev[i][k] == SeqTools::gapIdx());
          if (!cEnd && (cTermPad[i][k] != SeqTools::gapIdx())) {
            numTot++;
            if (cTermPad[i][k] == cTermPadPrev[i][k]) numID++;
          }
          if (numTot >= contextLength) break;
          int kn = nTermPad[i].size() - k - 1;
          nEnd = nEnd || (nTermPad[i][kn] == SeqTools::gapIdx()) || (nTermPadPrev[i][kn] == SeqTools::gapIdx());
          if ((!nEnd) && (nTermPad[i][kn] != SeqTools::gapIdx())) {
            numTot++;
            if (nTermPad[i][kn] == nTermPadPrev[i][kn]) numID++;
          }
          if (numTot >= contextLength) break;
          if (cEnd && nEnd) break;
        }
        if (isWithinSeqID(contextLength, redundancyCut, numTot, numID)) {
          // psol and sol are redundant. Which should go?
          if (psol->getRMSD() <= sol.getRMSD()) { // sol got trumped by a better previous solution
            if (algnBarInfoSet && (algnRedBar[i][sol[i]] > psol->getRMSD())) {
              algnRedBar[i][sol[i]] = psol->getRMSD();
              algnRedBarSource[i][psol].insert(sol[i]);
            }
            return false;
          } else {
            // This previous solution gets trumped by sol, but we need to find the
            // best of the similar solutions to see if any trump sol. If none trump
            // sol, then the best of the similar previous solutions will be removed.
            // NOTE: why remove the best??? maybe the worst? or the one with the clostst RMSD?
            if ((toRemove == NULL) || (toRemove->getRMSD() > psol->getRMSD())) toRemove = psol;
            break; // no need to check other segments
          }
        }
      }
    }
  }
  if (toRemove != NULL) {
    if (algnBarInfoSet) {
      for (int i = 0; i < sol.numSegments(); i++) {
        if (algnRedBarSource[i].find(toRemove) != algnRedBarSource[i].end()) {
          set<int>& wasTrumping = algnRedBarSource[i][toRemove];
          for (int j : wasTrumping) algnRedBar[i][j] = INFINITY;
        }
        algnRedBarSource[i].erase(toRemove);
      }
    }
    erase(*toRemove);
  }
  pair<set<fasstSolution>::iterator, bool> ins = solsSet.insert(sol);
  // NOTE: I currently do not believe that solsByCenRes should be returned to the user,
  // but if I change my mind later, this will need to be uncommented; also see fasstSolutionSet::operator=
  // fasstSolution* inserted = (fasstSolution*) &(*(ins.first));
  // for (int i = 0; i < sol.numSegments(); i++) solsByCenRes[i][sol.segCentralResidue(i)].insert(inserted);
  updated = true;
  return true;
}

bool fasstSolutionSet::insert(const fasstSolution& sol, simpleMap<resAddress, tightvector<resAddress>>& relMap) {
  // apply a redudancy filter based on a pre-computed map of inter-residue relationships
  fasstSolution* toRemove = NULL;
  bool algnBarInfoSet = isAlignRedBarrierDataSet();
  for (int i = 0; i < sol.numSegments(); i++) {
    resAddress ri = sol.segCentralResidue(i);
    tightvector<resAddress> simList;
    if (relMap.find(ri) >= 0) simList = relMap[ri];
    for (int j = 0; j <= simList.size(); j++) {
      // make sure to also look for solutions involving the same central residue
      // as in the current solution (in case redundancy with self is not included
      // in the map, which would be wasteful)
      resAddress& rj = (j < simList.size()) ? simList[j] : ri;
      if (solsByCenRes[i].find(rj) == solsByCenRes[i].end()) continue;
      set<fasstSolution*>& simSols = solsByCenRes[i][rj];
      for (auto it = simSols.begin(); it != simSols.end(); ++it) {
        fasstSolution* psol = (fasstSolution*) *it;
        if (psol->getRMSD() <= sol.getRMSD()) {
          if (algnBarInfoSet && (algnRedBar[i][sol[i]] > psol->getRMSD())) {
            algnRedBar[i][sol[i]] = psol->getRMSD();
            algnRedBarSource[i][psol].insert(sol[i]);
          }
          return false; // sol got trumped by a better previous solution
        } else {
          // This previous solution gets trumped by sol, but we need to find the
          // best of the similar solutions to see if any trump sol. If none trump
          // sol, then the best of the similar previous solutions will be removed.
          // NOTE: why remove the best??? maybe the worst? or the one with the clostst RMSD?
          if ((toRemove == NULL) || (toRemove->getRMSD() > psol->getRMSD())) toRemove = psol;
        }
      }
    }
  }

  if (toRemove != NULL) {
    if (algnBarInfoSet) {
      for (int i = 0; i < sol.numSegments(); i++) {
        if (algnRedBarSource[i].find(toRemove) != algnRedBarSource[i].end()) {
          set<int>& wasTrumping = algnRedBarSource[i][toRemove];
          for (int j : wasTrumping) algnRedBar[i][j] = INFINITY;
        }
        algnRedBarSource[i].erase(toRemove);
      }
    }
    erase(*toRemove);
  }
  pair<set<fasstSolution>::iterator, bool> ins = solsSet.insert(sol);
  fasstSolution* inserted = (fasstSolution*) &(*(ins.first));
  for (int i = 0; i < sol.numSegments(); i++) solsByCenRes[i][sol.segCentralResidue(i)].insert(inserted);
  updated = true;
  return true;
}

void fasstSolutionSet::erase(fasstSolution& sol) {
  int ti = sol.getTargetIndex();
  if (!solsByCenRes.empty()) {
    for (int i = 0; i < sol.numSegments(); i++) {
      solsByCenRes[i][sol.segCentralResidue(i)].erase(&sol);
      if (solsByCenRes[i][sol.segCentralResidue(i)].empty()) solsByCenRes[i].erase(sol.segCentralResidue(i));
    }
  }
  solsSet.erase(sol);
  updated = true;
}

set<fasstSolution>::iterator fasstSolutionSet::erase(const set<fasstSolution>::iterator it) {
  int ti = it->getTargetIndex();
  fasstSolution* solPtr = (fasstSolution*) &(*it);
  fasstSolution& sol = *solPtr;
  for (int i = 0; i < sol.numSegments(); i++) {
    solsByCenRes[i][sol.segCentralResidue(i)].erase(&sol);
    if (solsByCenRes[i][sol.segCentralResidue(i)].empty()) solsByCenRes[i].erase(sol.segCentralResidue(i));
  }
  updated = true;
  return solsSet.erase(it);
}

fasstSolution& fasstSolutionSet::operator[] (int i) {
  // if solutions changed since last time, first copy to vector
  if (updated) {
    solsVec.resize(solsSet.size(), NULL); int k = 0;
    for (auto it = solsSet.begin(); it != solsSet.end(); ++it, ++k) solsVec[k] = (fasstSolution*) &(*it);
    updated = false;
  }
  return *(solsVec[i]);
}

bool fasstSolutionSet::isWithinSeqID(int L0, mstreal cut, int numTot, int numID) {
  // if alignment is not long enough, interpolate the cutoff
  if (numTot < L0) {
    int lowL = 5; mstreal lowCut = 1.0;
    if (numTot <= lowL) { cut = lowCut; }
    else {
      mstreal a = (cut - lowCut)/(L0 - lowL);
      cut = a*(numTot - lowL) + lowCut;
    }
  }
  return (numID >= numTot*cut);
}

vector<fasstSolution*> fasstSolutionSet::orderByDiscovery() {
  vector<fasstSolution*> sols(solsSet.size(), NULL);
  int k = 0;
  for (auto it = solsSet.begin(); it != solsSet.end(); ++it, ++k) {
    sols[k] = (fasstSolution*) &(*it);
  }
  sort(sols.begin(), sols.end(), fasstSolution::foundBefore);
  return sols;
}

vector<fasstSolutionAddress> fasstSolutionSet::extractAddresses() const {
  vector<fasstSolutionAddress> addresses(size()); int i = 0;
  for (auto it = begin(); it != end(); ++it, ++i) addresses[i] = it->getAddress();
  return addresses;
}

void fasstSolutionSet::write(ostream &_os) const {
  MstUtils::writeBin(_os, (int) solsSet.size());
  for (auto it = solsSet.begin(); it != solsSet.end(); ++it) it->write(_os);
}

void fasstSolutionSet::read(istream &_is) {
  updated = true;
  int len; MstUtils::readBin(_is, len);
  solsSet.clear();
  for (int i = 0; i < len; i++) {
    fasstSolution sol;
    sol.read(_is);
    insert(sol);
  }
}


/* --------- fasstSeqConstSimple --------- */
void fasstSeqConstSimple::evalConstraint(int segIdx, const Sequence& target, vector<bool>& alignments) {
  vector<int>& segPositions = positions[segIdx];
  if (segPositions.empty()) {
    for (int i = 0; i < alignments.size(); i++) alignments[i] = true;
    return;
  }
  vector<set<res_t> > segAminoAcids = aminoAcids[segIdx];
  for (int i = 0; i < alignments.size(); i++) {
    bool ok = true;
    for (int j = 0; j < segPositions.size(); j++) {
      int k = i + segPositions[j];
      ok = ok && (segAminoAcids[j].find(target[k]) != segAminoAcids[j].end());
    }
    alignments[i] = ok;
  }
}
