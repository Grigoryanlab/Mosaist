#include "mstfasst.h"

/* --------- FASST::optList --------- */
void FASST::optList::setOptions(const vector<mstreal>& _costs, bool add) {
  // sort costs and keep track of indices to know the rank-to-index mapping
  costs = _costs;
  rankToIdx.resize(_costs.size());
  for (int i = 0; i < rankToIdx.size(); i++) rankToIdx[i] = i;
  sort(rankToIdx.begin(), rankToIdx.end(), [this](int i, int j) { return costs[i] < costs[j]; });
  for (int i = 0; i < rankToIdx.size(); i++) costs[i] = _costs[rankToIdx[i]];

  // sort the the rank-to-index mapping and keep track of indices to know the index-to-rank mapping
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
  // if ((k < 0) || (k >= costs.size())) MstUtils::error("out-of-range index specified: " + MstUtils::toString(k), "FASST::optList::insertOption(int)");
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
  vector<bool> isGiven(isIn.size(), false);
  for (int i = 0; i < opts.size(); i++) isGiven[opts[i]] = true;
  for (int i = 0; i < isIn.size(); i++) {
    if (isIn[i] && !isGiven[i]) removeOption(i);
  }
}

void FASST::optList::constrainLE(int idx) {
  for (int i = idx+1; i < isIn.size(); i++) removeOption(i);
}

void FASST::optList::constrainGE(int idx) {
  for (int i = 0; i < MstUtils::min(idx, (int) isIn.size()); i++) removeOption(i);
}

void FASST::optList::constrainRange(int idxLow, int idxHigh) {
  constrainGE(idxLow);
  constrainLE(idxHigh);
}

/* --------- FASST --------- */
FASST::FASST() {
  recLevel = 0;
  setRMSDCutoff(1.0);
  setSearchType(searchType::FULLBB);
  querySize = 0;
  updateGrids = false;
  gridSpacing = 15.0;
  memSave = false;
  maxNumMatches = minNumMatches = suffNumMatches = -1;
  gapConstSet = false;
  contextLength = 30;
  redundancyCut = 1.0;
}

FASST::~FASST() {
  // need to delete atoms only on the lowest level of recursion, because at
  // higher levels we point to the same atoms
  if (targetMasks.size()) targetMasks.back().deletePointers();
  for (int i = 0; i < targetStructs.size(); i++) delete targetStructs[i];
  for (int i = 0; i < ps.size(); i++) delete ps[i];
}

void FASST::setCurrentRMSDCutoff(mstreal cut) {
  rmsdCut = cut;
  residualCut = cut*cut*querySize;
}

void FASST::setMaxNumMatches(int _max) {
  maxNumMatches = _max;
  if (!areNumMatchConstraintsConsistent()) MstUtils::error("invalid combination of match number constraints: [min, max, sufficient] = [" + MstUtils::toString(minNumMatches) + ", " + MstUtils::toString(maxNumMatches) + ", " + MstUtils::toString(suffNumMatches) + "]", "FASST::setMaxNumMatches");
}

void FASST::setMinNumMatches(int _min) {
  minNumMatches = _min;
  if (!areNumMatchConstraintsConsistent()) MstUtils::error("invalid combination of match number constraints: [min, max, sufficient] = [" + MstUtils::toString(minNumMatches) + ", " + MstUtils::toString(maxNumMatches) + ", " + MstUtils::toString(suffNumMatches) + "]", "FASST::setMinNumMatches");
}

void FASST::setSufficientNumMatches(int _suff) {
  suffNumMatches = _suff;
  if (!areNumMatchConstraintsConsistent()) MstUtils::error("invalid combination of match number constraints: [min, max, sufficient] = [" + MstUtils::toString(minNumMatches) + ", " + MstUtils::toString(maxNumMatches) + ", " + MstUtils::toString(suffNumMatches) + "]", "FASST::setSufficientNumMatches");
}

bool FASST::areNumMatchConstraintsConsistent() {
  if (isMaxNumMatchesSet() && isMinNumMatchesSet() && (minNumMatches > maxNumMatches)) return false;
  if (isMaxNumMatchesSet() && isSufficientNumMatchesSet() && (maxNumMatches < suffNumMatches)) return false;
  if (isMinNumMatchesSet() && isSufficientNumMatchesSet() && (minNumMatches > suffNumMatches)) return false;
  return true;
}

void FASST::setMinGap(int i, int j, int gapLim) {
  if (gapLim < 0) MstUtils::error("gap constraints must be defined in the non-negative direction", "FASST::setMinGap");
  minGap[i][j] = gapLim;
  minGapSet[i][j] = true;
  gapConstSet = true;
}

void FASST::setMaxGap(int i, int j, int gapLim) {
  if (gapLim < 0) MstUtils::error("gap constraints must be defined in the non-negative direction", "FASST::setMinGap");
  maxGap[i][j] = gapLim;
  maxGapSet[i][j] = true;
  gapConstSet = true;
}

void FASST::resetGapConstraints() {
  minGap.resize(query.size(), vector<int>(query.size(), 0));
  maxGap = minGap;
  minGapSet.resize(query.size(), vector<bool>(query.size(), false));
  maxGapSet = minGapSet;
  gapConstSet = false;
}

bool FASST::validateSearchRequest() {
  if (query.size() != minGap.size()) {
    MstUtils::error("gap constraints inconsistent with number of segments in query", "FASST::validateSearchRequest()");
  }
  return true;
}

/* This function sets up the query, in the process deciding which part of the
 * query is really searchable (e.g., backbone). Various search types an be added
 * in the future to provide search capabilities over different parts of the
 * structure. The current implementation stipulates that the searchable part
 * involve the same number of atoms for each residue, so that atoms can be
 * stored in flat arrays for efficiency. But these could presumably be
 * interpreted differently. E.g., only some residue matches (parents of atoms)
 * may be accepted or some of the atoms can be dummy/empty ones. */
void FASST::setQuery(const string& pdbFile) {
  queryStruct.reset();
  queryStruct.readPDB(pdbFile);
  processQuery();
}

void FASST::setQuery(const Structure& Q) {
  queryStruct = Q;
  processQuery();
}

void FASST::processQuery() {
  querySize = 0;
  // auto-splitting segments by connectivity, but can do diffeerntly
  queryStruct = queryStruct.reassignChainsByConnectivity();
  query.resize(queryStruct.chainSize());
  for (int i = 0; i < queryStruct.chainSize(); i++) {
    query[i].resize(0);
    if (!parseChain(queryStruct[i], query[i])) {
      MstUtils::error("could not set query, because some atoms for the specified search type were missing", "FASST::setQuery");
    }
    MstUtils::assert(query[i].size() > 0, "query contains empty segment(s)", "FASST::setQuery");
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
  centToCentDist.resize(query.size(), vector<mstreal>(query.size(), 0));
  CartesianPoint C(0, 0, 0);
  int N = 0;
  for (int L = 0; L < query.size(); L++) {
    CartesianPoint ci = query[L].getGeometricCenter();
    int n = query[L].size();
    C = (C*N + ci*n)/(N + n);
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
  resetGapConstraints();
}

void FASST::addTarget(const string& pdbFile) {
  Structure* targetStruct = new Structure(pdbFile, "QUIET");
  targetSource.push_back(targetInfo(pdbFile, targetFileType::PDB, 0, memSave));
  addTargetStructure(targetStruct);
}

void FASST::addTarget(const Structure& T) {
  Structure* targetStruct = new Structure(T);
  targetSource.push_back(targetInfo(T.getName(), targetFileType::STRUCTURE, 0, memSave));
  addTargetStructure(targetStruct);
}

void FASST::addTargetStructure(Structure* targetStruct) {
  if (memSave) stripSidechains(*targetStruct);
  targetStructs.push_back(targetStruct);
  targets.push_back(AtomPointerVector());
  targSeqs.push_back(Sequence());
  AtomPointerVector& target = targets.back();
  Sequence& seq = targSeqs.back();
  // we don't care about the chain topology of the target, so append all residues
  for (int i = 0; i < targetStruct->chainSize(); i++) {
    parseChain(targetStruct->getChain(i), target, &seq);
  }
  MstUtils::assert(target.size() > 0, "empty target named '" + targetStruct->getName() + "'", "FASST::addTargetStructure");
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
}

void FASST::addTargets(const vector<string>& pdbFiles) {
  for (int i = 0; i < pdbFiles.size(); i++) addTarget(pdbFiles[i]);
}

bool FASST::parseChain(const Chain& C, AtomPointerVector& searchable, Sequence* seq) {
  bool foundAll = true;
  for (int i = 0; i < C.residueSize(); i++) {
    Residue& res = C.getResidue(i);
    switch(type) {
      case searchType::CA:
      case searchType::FULLBB: {
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
          searchable.insert(searchable.end(), bb.begin(), bb.end());
          if (seq != NULL) seq->appendResidue(res.getName());
        }
        break;
      }
      default:
        MstUtils::error("uknown search type '" + MstUtils::toString(type) + "' specified", "FASST::parseStructure");
    }
  }
  return foundAll;
}

void FASST::addResidueProperties(int ti, const string& propType, const vector<mstreal>& propVals) {
  if ((ti < 0) || (ti >= targetStructs.size())) MstUtils::error("requested target out of range: " + MstUtils::toString(ti), "FASST::addResidueProperties");
  int N = targetStructs[ti]->residueSize();
  if (N != propVals.size()) MstUtils::error("size of properties vector inconsistent with number of residues for target: " + MstUtils::toString(ti), "FASST::addResidueProperties");
  vector<mstreal>& vals = resProperties[propType][ti];
  vals.resize(propVals.size());
  for (int i = 0; i < propVals.size(); i++) vals[i] = propVals[i];
}

void FASST::writeDatabase(const string& dbFile) {
  fstream ofs; MstUtils::openFile(ofs, dbFile, fstream::out | fstream::binary, "FASST::writeDatabase");
  for (int ti = 0; ti < targetStructs.size(); ti++) {
    MstUtils::writeBin(ofs, 'S'); // marks the start of a structure section
    targetStructs[ti]->writeData(ofs);
    for (auto p = resProperties.begin(); p != resProperties.end(); ++p) {
      if ((p->second).find(ti) != (p->second).end()) {
        vector<mstreal>& vals = (p->second)[ti];
        MstUtils::writeBin(ofs, 'P'); // marks the start of a residue property section
        MstUtils::writeBin(ofs, (string) p->first);
        MstUtils::assert(targetStructs[ti]->residueSize() == vals.size(), "the number of residue properties and residues does not agree for database entry", "FASST::writeDatabase(const string&)");
        for (int ri = 0; ri < vals.size(); ri++) MstUtils::writeBin(ofs, vals[ri]);
      }
    }
  }
  ofs.close();
}

void FASST::readDatabase(const string& dbFile) {
  fstream ifs; MstUtils::openFile(ifs, dbFile, fstream::in | fstream::binary, "FASST::readDatabase");
  char sect; string name; mstreal val; int ti = 0;
  MstUtils::readBin(ifs, sect);
  if (sect != 'S') MstUtils::error("first section must be a structure one, while reading database file " + dbFile, "FASST::readDatabase(const string&)");
  while (ifs.peek() != EOF) {
    Structure* targetStruct = new Structure();
    int loc = ifs.tellg();
    targetStruct->readData(ifs);
    targetSource.push_back(targetInfo(dbFile, targetFileType::BINDATABASE, loc, memSave));
    addTargetStructure(targetStruct);
    while (ifs.peek() != EOF) {
      MstUtils::readBin(ifs, sect);
      if (sect == 'P') {
        int N = targetStruct->residueSize();
        MstUtils::readBin(ifs, name);
        vector<mstreal>& vals = resProperties[name][ti];
        vals.resize(N, 0);
        for (int i = 0; i < vals.size(); i++) {
          MstUtils::readBin(ifs, val);
          vals[i] = val;
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
    Chain& C = S[i];
    for (int j = 0; j < C.residueSize(); j++) {
      Residue& R = C[j];
      int n = R.atomSize();
      for (int k = 0; k < n; k++) {
        char* name = R[k].getNameC();
        if (strlen(name)) { // check that it's okay to index atom name string
          switch (name[0]) {
            case 'N':
              if (!strcmp(name, "N") || !strcmp(name, "NT")) continue;
              break;
            case 'C':
              if (!strcmp(name, "C") || !strcmp(name, "CA")) continue;
              break;
            case 'O':
              if (!strcmp(name, "O") || !strcmp(name, "OT1") || !strcmp(name, "OT2") || !strcmp(name, "OXT")) continue;
              break;
          }
        }
        R.deleteAtom(k);
        k--; n--;
      }
      R.compactify();
    }
  }
}

void FASST::rebuildProximityGrids() {
  if (xlo == xhi) { xlo -= gridSpacing/2; xhi += gridSpacing/2; }
  if (ylo == yhi) { ylo -= gridSpacing/2; yhi += gridSpacing/2; }
  if (zlo == zhi) { zlo -= gridSpacing/2; zhi += gridSpacing/2; }
  int N = int(ceil(max(max((xhi - xlo), (yhi - ylo)), (zhi - zlo))/gridSpacing));
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
  for (int i = 0; i < query.size(); i++) {
    ps[i]->dropAllPoints();
    AtomPointerVector& seg = query[i];
    int Na = atomToResIdx(target.size()) - atomToResIdx(seg.size()) + 1; // number of possible alignments
    segmentResiduals[i].resize(MstUtils::max(Na, 0));
    for (int j = 0; j < Na; j++) {
      // NOTE: can save on this in several ways:
      // 1. the centroid calculation is effectively already done inside RMSDCalculator::bestRMSD
      // 2. updating just one atom involves a simple centroid adjustment, rather than recalculation
      // 3. is there a speedup to be gained from re-calculating RMSD with one atom updated only?
      AtomPointerVector targSeg = target.subvector(resToAtomIdx(j), resToAtomIdx(j) + query[i].size());
      segmentResiduals[i][j] = RC.bestResidual(query[i], targSeg);
      targSeg.getGeometricCenter(xc, yc, zc);
      if (query.size() > 1) ps[i]->addPoint(xc, yc, zc, j);
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
  if (gapConstraintsExist()) {
    targChainBeg.resize(atomToResIdx(target.size()), 0);
    for (int i = 1; i < targChainBeg.size(); i++) {
      if (target[resToAtomIdx(i)]->getChain() == target[resToAtomIdx(i-1)]->getChain()) targChainBeg[i] = targChainBeg[i-1];
      else targChainBeg[i] = i;
    }
    targChainEnd.resize(atomToResIdx(target.size()), atomToResIdx(target.size()) - 1);
    for (int i = targChainEnd.size() - 2; i >= 0; i--) {
      if (target[resToAtomIdx(i)]->getChain() == target[resToAtomIdx(i+1)]->getChain()) targChainEnd[i] = targChainEnd[i+1];
      else targChainEnd[i] = i;
    }
  } else {
    targChainBeg.resize(0); targChainEnd.resize(0);
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

void FASST::search() {
  // auto begin = chrono::high_resolution_clock::now();
  // int prepTime = 0;
  validateSearchRequest();
  if (isMinNumMatchesSet()) setCurrentRMSDCutoff(999.0);
  else setCurrentRMSDCutoff(rmsdCutRequested);
  solutions.clear();
  vector<int> segLen(query.size()); // number of residues in each query segment
  for (int i = 0; i < query.size(); i++) segLen[i] = atomToResIdx(query[i].size());
  vector<mstreal> ccTol(query.size(), -1.0);
  for (currentTarget = 0; currentTarget < targets.size(); currentTarget++) {
    // auto beginPrep = chrono::high_resolution_clock::now();
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

      // 2. compute the total residual from the current alignment
      mstreal curBound = currentAlignmentResidual(true) + boundOnRemainder(true);
      if (curBound > residualCut) continue;
      // if (query.size() > 1) updateQueryCentroids();

      // 3. update update remaining options for subsequent segments based on the
      // newly made choice. The set of options on the next recursion level is a
      // subset of the set of options on the previous level.
      int remSegs = query.size() - (recLevel + 1);
      if (remSegs > 0) {
        bool levelExhausted = false;
        int nextLevel = recLevel + 1;
        // copy remaining options from the previous recursion level. This way,
        // we can compute bounds on this level and can do set intersections to
        // further narrow this down
        for (int i = nextLevel; i < query.size(); i++) {
          remOptions[nextLevel][i].copyIn(remOptions[nextLevel-1][i]);
          // except that segments cannot overlap, so remove from consideration
          // all alignments that overlap with the segments that was just placed
          remOptions[nextLevel][i].removeOptions(currAlignment[recLevel] - segLen[i] + 1,
                                                 currAlignment[recLevel] + segLen[recLevel] - 1);
        }
        // if any gap constraints exist, limit options at this recursion level accordingly
        if (gapConstraintsExist()) {
          for (int j = 0; j < nextLevel; j++) {
            for (int i = nextLevel; i < query.size(); i++) {
              remOptions[nextLevel][i].constrainRange(targChainBeg[currAlignment[j]], targChainEnd[currAlignment[j]]);
              if (minGapSet[qSegOrd[i]][qSegOrd[j]]) remOptions[nextLevel][i].constrainLE(currAlignment[j] - minGap[qSegOrd[i]][qSegOrd[j]] - segLen[i]);
              if (maxGapSet[qSegOrd[i]][qSegOrd[j]]) remOptions[nextLevel][i].constrainGE(currAlignment[j] - maxGap[qSegOrd[i]][qSegOrd[j]] - segLen[i]);
              if (minGapSet[qSegOrd[j]][qSegOrd[i]]) remOptions[nextLevel][i].constrainGE(currAlignment[j] + minGap[qSegOrd[j]][qSegOrd[i]] + segLen[j]);
              if (maxGapSet[qSegOrd[j]][qSegOrd[i]]) remOptions[nextLevel][i].constrainLE(currAlignment[j] + maxGap[qSegOrd[j]][qSegOrd[i]] + segLen[j]);
              // if (minGapSet[qSegOrd[i]][qSegOrd[j]]) remOptions[nextLevel][i].constrainRange(targChainBeg[currAlignment[j]], currAlignment[j] - minGap[qSegOrd[i]][qSegOrd[j]] - segLen[i]);
              // if (maxGapSet[qSegOrd[i]][qSegOrd[j]]) remOptions[nextLevel][i].constrainRange(currAlignment[j] - maxGap[qSegOrd[i]][qSegOrd[j]] - segLen[i], targChainEnd[currAlignment[j]]);
              // if (minGapSet[qSegOrd[j]][qSegOrd[i]]) remOptions[nextLevel][i].constrainRange(currAlignment[j] + minGap[qSegOrd[j]][qSegOrd[i]] + segLen[j], targChainEnd[currAlignment[j]]);
              // if (maxGapSet[qSegOrd[j]][qSegOrd[i]]) remOptions[nextLevel][i].constrainRange(targChainBeg[currAlignment[j]], currAlignment[j] + maxGap[qSegOrd[j]][qSegOrd[i]] + segLen[j]);
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
          for (int i = nextLevel; i < query.size(); i++) {
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
        // if at the lowest recursion level already, then record the solution
        fasstSolution sol(currAlignment, sqrt(currResidual/querySize), currentTarget, solutions.size(), qSegOrd);
        if (redundancyCut < 1) addSequenceContext(sol, currentTarget, segLen, contextLength);
        solutions.insert(sol, redundancyCut);
        if (isSufficientNumMatchesSet() && solutions.size() == suffNumMatches) return;
        if (isMaxNumMatchesSet() && (solutions.size() > maxNumMatches)) {
          solutions.erase(--solutions.end());
          setCurrentRMSDCutoff(solutions.rbegin()->getRMSD());
        } else if (isMinNumMatchesSet() && (solutions.size() > minNumMatches) && (rmsdCut > rmsdCutRequested)) {
          if (solutions.rbegin()->getRMSD() > rmsdCutRequested) solutions.erase(--solutions.end());
          setCurrentRMSDCutoff(solutions.rbegin()->getRMSD());
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
}

mstreal FASST::currentAlignmentResidual(bool compute) {
  if (compute) {
    if (query.size() == 1) {
      // this is a special case, because will not need to calculate centroid
      // locations for subsequent sub-queries
      currResidual = segmentResiduals[0][currAlignment[0]];
    } else {
      // fill up sub-alignment with target atoms
      int N = targetMasks[recLevel].size();
      int n = query[recLevel].size();
      int si = resToAtomIdx(currAlignment[recLevel]);
      AtomPointerVector& target = targets[currentTarget];
      for (int i = 0; i < n; i++) {
        targetMasks[recLevel][N - n + i]->setCoor(target[si + i]->getX(), target[si + i]->getY(), target[si + i]->getZ());
      }
      currResidual = RC.bestResidual(queryMasks[recLevel], targetMasks[recLevel]);
      currResiduals[recLevel] = currResidual;
      if (recLevel == 0) {
        currCents[recLevel] = ps[recLevel]->getPoint(currAlignment[recLevel]);
      } else {
        currCents[recLevel] = (currCents[recLevel - 1] * (N - n) +  ps[recLevel]->getPoint(currAlignment[recLevel]) * n) / N;
      }
    }
  }
  return currResidual;
}

void FASST::getMatchStructure(const fasstSolution& sol, Structure& match, bool detailed, matchType type) {
  vector<Structure> matches;
  fasstSolutionSet solSet; solSet.insert(sol);
  getMatchStructures(solSet, matches, detailed, type);
  match = matches[0];
}

Structure FASST::getMatchStructure(const fasstSolution& sol, bool detailed, matchType type) {
  Structure match; getMatchStructure(sol, match, detailed, type); return match;
}

void FASST::getMatchStructures(fasstSolutionSet& sols, vector<Structure>& matches, bool detailed, matchType type) {
  // hash solutions by the target they come from, to visit each target only once
  map<int, vector<int> > solsFromTarget;
  for (int i = 0; i < sols.size(); i++) {
    fasstSolution& sol = sols[i];
    int idx = sol.getTargetIndex();
    if ((idx < 0) || (idx >= targets.size())) {
      MstUtils::error("supplied FASST solution is pointing to an out-of-range target", "FASST::getMatchStructures");
    }
    solsFromTarget[idx].push_back(i);
  }

  // flatten query into a single array of atoms
  AtomPointerVector queryAtoms;
  for (int k = 0; k < queryOrig.size(); k++) {
    for (int ai = 0; ai < queryOrig[k].size(); ai++) queryAtoms.push_back(queryOrig[k][ai]);
  }
  AtomPointerVector matchAtoms; matchAtoms.reserve(queryAtoms.size());

  // visit each target
  Structure dummy; RMSDCalculator rc;
  matches.resize(sols.size(), Structure());
  for (auto it = solsFromTarget.begin(); it != solsFromTarget.end(); ++it) {
    int idx = it->first;
    Structure* targetStruct = targetStructs[idx];
    AtomPointerVector& target = targets[idx];
    bool reread = detailed && targetSource[idx].memSave; // should we re-read the target structure?
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
      tr[idx].apply(dummy);
      targetStruct = &dummy;
    }

    // visit each solution from this target
    vector<int>& solIndices = it->second;
    for (int i = 0; i < solIndices.size(); i++) {
      int solIndex = solIndices[i];
      const fasstSolution& sol = sols[solIndex];
      Structure& match = matches[solIndex];
      match.setName(targetStruct->getName());
      vector<int> alignment = sol.getAlignment();
      if (alignment.size() != queryOrig.size()) {
        MstUtils::error("solution alignment size inconsistent with number of query segments", "FASST::getMatchStructures");
      }

      // excise match atoms that superimpose onto the query, in the right order
      matchAtoms.resize(0);
      for (int k = 0; k < alignment.size(); k++) {
        // copy matching residues from the original target structure
        int si = resToAtomIdx(alignment[k]);
        int L = queryOrig[k].size();
        if (si + L > target.size()) {
          MstUtils::error("solution points to atom outside of target range", "FASST::getMatchStructures");
        }
        for (int ai = si; ai < si + L; ai++) matchAtoms.push_back(target[ai]);
      }

      // cut out the part of the target Structure that will constitute the returned match
      if (type == matchType::FULL) {
        match = *targetStruct;
      } else {
        vector<int> resIndices = getMatchResidueIndices(sol, type);
        for (auto ri = resIndices.begin(); ri != resIndices.end(); ri++) {
          Residue* res = target[resToAtomIdx(*ri)]->getResidue();
          if (reread) res = &(targetStruct->getResidue(res->getResidueIndex()));
          match.addResidue(res);
        }
      }

      // align matching region onto query, transforming the match itself
      rc.align(matchAtoms, queryAtoms, match);
    }
  }
}

vector<Sequence> FASST::getMatchSequences(fasstSolutionSet& sols, matchType type) {
  vector<Sequence> seqs(sols.size());
  for (int i = 0; i < sols.size(); i++) {
    fasstSolution& sol = sols[i];
    int idx = sol.getTargetIndex();
    MstUtils::assert((idx >= 0) && (idx < targSeqs.size()), "supplied FASST solution is pointing to an out-of-range target", "FASST::getMatchSequences");
    vector<int> alignment = sol.getAlignment();
    seqs[i].setName(targSeqs[idx].getName());

    // isolate out the part of the target Sequence that will constitute the returned match
    vector<int> resIndices = getMatchResidueIndices(sol, type);
    for (auto ri = resIndices.begin(); ri != resIndices.end(); ri++) seqs[i].appendResidue(targSeqs[idx][*ri]);
  }
  return seqs;
}

Sequence FASST::getMatchSequence(const fasstSolution& sol, matchType type) {
  fasstSolutionSet solSet; solSet.insert(sol);
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
    fasstSolution& sol = sols[i];
    int idx = sol.getTargetIndex();
    MstUtils::assert((idx >= 0) && (idx < targets.size()), "supplied FASST solution is pointing to an out-of-range target", "FASST::getMatchSequences");
    AtomPointerVector& target = targets[idx];
    if (!isPropertyDefined(propType, idx)) {
      MstUtils::error("target with index " + MstUtils::toString(idx) + " does not have property type " + propType, "FASST::getResidueProperties(fasstSolutionSet&, const string&, matchType)");
    }
    vector<mstreal>& propVals = resProperties[propType][idx];
    vector<int> resIndices = getMatchResidueIndices(sol, type);
    props[i].resize(resIndices.size()); int ii = 0;
    for (auto ri = resIndices.begin(); ri != resIndices.end(); ri++, ii++) {
      // properties are defined for each residue in the original target structure
      Residue* res = target[resToAtomIdx(*ri)]->getResidue();
      props[i][ii] = propVals[res->getResidueIndex()];
    }
  }
  return props;
}

bool FASST::isPropertyDefined(const string& propType, int ti) {
  return ((resProperties.find(propType) != resProperties.end()) || (resProperties[propType].find(ti) != resProperties[propType].end()));
}

bool FASST::isPropertyDefined(const string& propType) {
  return (resProperties.find(propType) != resProperties.end());
}

vector<int> FASST::getMatchResidueIndices(const fasstSolution& sol, matchType type) {
  vector<int> residueIndices;
  vector<int> alignment = sol.getAlignment();
  switch(type) {
    case matchType::REGION: {
      for (int k = 0; k < alignment.size(); k++) {
        int si = alignment[k];
        int L = atomToResIdx(queryOrig[k].size());
        for (int ri = si; ri < si + L; ri++) residueIndices.push_back(ri);
      }
      break;
    }
    case matchType::WITHGAPS: {
      set<int> toInclude;
      for (int i = 0; i < queryOrig.size(); i++) {
        // each individual segments should be included
        for (int k = alignment[i]; k < alignment[i] + atomToResIdx(queryOrig[i].size()); k++) toInclude.insert(k);
        // then some gaps also
        for (int j = 0; j < queryOrig.size(); j++) {
          if (i == j) continue;
          // if gap constrained, fill in between these two segments
          if (gapConstrained(i, j)) {
            if (alignment[i] > alignment[j]) MstUtils::error("solution not consistent with current gap constraints", "FASST::getMatchResidueIndices");
            for (int k = alignment[i] + atomToResIdx(queryOrig[i].size()); k < alignment[j]; k++) {
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
      MstUtils::assert((idx >= 0) && (idx < targetStructs.size()), "supplied FASST solution is pointing to an out-of-range target", "FASST::getMatchResidueIndices");
      for (int i = 0; i < targetStructs[idx]->residueSize(); i++) residueIndices.push_back(i);
      break;
    }
    default:
      MstUtils::error("unknown match output type", "FASST::getMatchResidueIndices");
  }

  return residueIndices;
}

string FASST::toString(const fasstSolution& sol) {
  stringstream ss;
  ss << std::setprecision(6) << std::fixed << sol.getRMSD() << " " << targetStructs[sol.getTargetIndex()]->getName() << " [" << MstUtils::vecToString(sol.getAlignment(), ", ") << "]";
  return ss.str();
}

void FASST::addSequenceContext(fasstSolution& sol, int currentTarget, const vector<int>& segLengths, int contextLength) {
  Sequence& targSeq = targSeqs[currentTarget];
  vector<int> alignment = sol.getAlignment();
  vector<Sequence> segs(sol.numSegments()), ntPad(sol.numSegments()), ctPad(sol.numSegments());
  for (int i = 0; i < sol.numSegments(); i++) {
    segs[i].resize(segLengths[i]);
    for (int ri = sol[i]; ri < sol[i] + segLengths[i]; ri++) {
      segs[i][ri - alignment[i]] = targSeq[ri];
    }
    if (contextLength > segs[i].size()) {
      // pad on both ends, to accommodate different scenarios
      int pad = contextLength - segs[i].size();
      ntPad[i].resize(pad, SeqTools::gapIdx()); ctPad[i].resize(pad, SeqTools::gapIdx());
      for (int ri = 0; ri < pad; ri++) {
        if (sol[i] + segLengths[i] + ri > targChainEnd[sol[i]]) break; // stop if reach chain C-terminus
        ctPad[i][ri] = targSeq[sol[i] + segLengths[i] + ri];
      }
      for (int ri = 0; ri < pad; ri++) {
        if (sol[i] - ri - 1 < targChainBeg[sol[i]]) break; // stop if reach chain N-terminus
        ntPad[i][ntPad[i].size() - ri - 1] = targSeq[sol[i] - ri - 1];
      }
    }
  }
  sol.setSeqContext(segs, ntPad, ctPad);
}


/* --------- fasstSolution --------- */
fasstSolution::fasstSolution(const vector<int>& _alignment, mstreal _rmsd, int _target, int _foundOrder, vector<int> segOrder) {
  alignment = _alignment; rmsd = _rmsd; targetIndex = _target; foundOrder = _foundOrder;
  if (!segOrder.empty()) {
    for (int i = 0; i < _alignment.size(); i++) alignment[segOrder[i]] = _alignment[i];
  }
  context = NULL;
}

fasstSolution::fasstSolution(const fasstSolution& _sol) {
  alignment = _sol.alignment; rmsd = _sol.rmsd; targetIndex = _sol.targetIndex;
  foundOrder = _sol.foundOrder;
  if (_sol.context != NULL) context = new solContext(*(_sol.context));
  else context = NULL;
}

void fasstSolution::setSeqContext(const vector<Sequence>& _segSeq, const vector<Sequence>& _nSeq, const vector<Sequence>& _cSeq) {
  if (context == NULL) context = new solContext();
  context->segSeq = _segSeq; context->nSeq = _nSeq; context->cSeq = _cSeq;
}

void fasstSolution::setStructContext(const vector<AtomPointerVector>& _segStr, const vector<AtomPointerVector>& _nStr, const vector<AtomPointerVector>& _cStr) {
  if (context == NULL) context = new solContext();
  context->segStr = _segStr; context->nStr = _nStr; context->cStr = _cStr;
}

/* --------- fasstSolutionSet --------- */
bool fasstSolutionSet::insert(const fasstSolution& sol, mstreal redundancyCut) {
  // apply a redudancy filter, if needed
  if ((redundancyCut < 1) && sol.seqContextDefined()) {
    vector<Sequence> segSeqs = sol.segmentSeqs();
    vector<Sequence> nTermPad = sol.nTermContext();
    vector<Sequence> cTermPad = sol.cTermContext();
    // compare this solution to each previously accepted solution
    for (auto psol = solsSet.begin(); psol != solsSet.end(); ) {
      vector<Sequence> segSeqsPrev = psol->segmentSeqs();
      vector<Sequence> nTermPadPrev = psol->nTermContext();
      vector<Sequence> cTermPadPrev = psol->cTermContext();
      // compare the contexts of each segment:
      bool psolStays = true;
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
          if (psol->getRMSD() <= sol.getRMSD()) { return false; } // sol got trumped by a better previous solution
          else {
            // sol is staying, but psol will need to be removed
            psol = solsSet.erase(psol); psolStays = false; break; // no need to check other segments
          }
        }
      }
      if (psolStays) psol++;
    }
  }
  solsSet.insert(sol);
  updated = true;
  return true;
}

fasstSolution& fasstSolutionSet::operator[] (int i) {
  // if solutions changed since last time, first copy to vector
  if (updated) {
    solsVec.resize(solsSet.size(), NULL); int i = 0;
    for (auto it = solsSet.begin(); it != solsSet.end(); ++it, ++i) solsVec[i] = (fasstSolution*) &(*it);
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
