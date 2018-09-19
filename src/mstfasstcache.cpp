#include "mstfasstcache.h"

/* --------- cFASST --------- */
void cFASST::write(const string& filename) const {
  fstream out;
  MstUtils::openFile(out, filename, fstream::out | fstream::binary, "cFASST::write");
  write(out);
  out.close();
}

void cFASST::write(ostream &_os) const {
  MstUtils::writeBin(_os, maxNumResults);
  MstUtils::writeBin(_os, errTolPressure);
  MstUtils::writeBin(_os, maxNumPressure);
  MstUtils::writeBin(_os, sf);
  MstUtils::writeBin(_os, (int) cache.size());
  for (auto it = cache.begin(); it != cache.end(); ++it) (*it)->write(_os);
}

void cFASST::read(const string& filename) {
  fstream in;
  MstUtils::openFile(in, filename, fstream::in | fstream::binary, "cFASST::read");
  read(in);
  in.close();
}

void cFASST::read(istream &_is) { // read object from a binary stream
  clear();
  MstUtils::readBin(_is, maxNumResults);
  MstUtils::readBin(_is, errTolPressure);
  MstUtils::readBin(_is, maxNumPressure);
  MstUtils::readBin(_is, sf);
  int numResults; MstUtils::readBin(_is, numResults);
  for (int i = 0; i < numResults; i++) {
    cachedResult* result = new cachedResult();
    result->read(_is);
    cache.insert(result);
  }
}

vector<int> cFASST::getStructureTopology(const Structure& S) {
  vector<int> topo(S.chainSize(), 0);
  for (int i = 0; i < S.chainSize(); i++) topo[i] = S[i].residueSize();
  return topo;
}

void cFASST::clear() {
  for (auto it = cache.begin(); it != cache.end(); ++it) delete(*it);
  cache.clear();
}

fasstSolutionSet cFASST::search() {
  int n;
  if (isSufficientNumMatchesSet()) MstUtils::error("cannot cache when \"sufficient\" number of matches is set", "cFASST::search");
  if (isMinNumMatchesSet()) MstUtils::error("cannot cache when \"minimum\" number of matches is set", "cFASST::search");
  chrono::high_resolution_clock::time_point begin, end;
  int searchTime;
  fasstSolutionSet matches;

  fasstSearchOptions origOpts = options();
  mstreal cut = getRMSDCutoff(), rmsd;
  bool maxSet = isMaxNumMatchesSet();
  bool redSet = isRedundancyCutSet();
  bool redPropSet = isRedundancyPropertySet();
  int maxN = getMaxNumMatches();
  mstreal redCut = getRedundancyCut();
  string redProp = getRedundancyProperty();
  AtomPointerVector queryAtoms = getQuerySearchedAtoms();

  /* Old cached results should eventually "expire", so uniformly lower priority
   * slightly first. This way, cached results that have not been used in a while
   * will eventually have a lower priority than brand new searches and these
   * will then push out these old (aparently) useless results. */
  if (modifyPermission()) {
    for (auto it = cache.begin(); it != cache.end(); it++) {
      // NOTE: gets rid of const qualifier! This is safe to do only because I know
      // I will monotonically lower everybody's priority (order will not change).
      cachedResult* res = &(*(*it));
      res->elapsePriority(getMaxNumResults());
    }
  }

  /* See whether all matches for the current query, within the given cutoff, are
   * among the list of matches of some previously cached query. */
  vector<int> topo = cFASST::getStructureTopology(getQuery());
  auto bestComp = cache.end(); mstreal bestDist = -1, safeRadius = -1;
  if (readPermission()) {
    for (auto it = cache.begin(); it != cache.end(); ++it) {
      cachedResult* result = *it;
      // TODO: should consider permutations! Both when comparing and when calculating
      // best RMSD. isSameTopology should compare sets. Then, depending on which
      // permutation has the best RMSD, return the order, and remap query atoms.
      // Need to write combinatorialAlign(const Structure& A, const Structure& B),
      // which will align one set of chains onto another. OR
      // combinatorialAlign(const vector<AtomPointerVector>& A, const vector<AtomPointerVector>& B)
      if (!result->isSameTopology(topo)) continue;
      mstreal r = rc.bestRMSD(queryAtoms, result->getQuery());
      /* If a maximum number of matches is set, then we may not need to find ALL
       * of the matches below the given cutoff, so just find the closest cached
       * query that has hopes of having ANY matches within the cutoff. If no max
       * is set on the number of matches, however, we will need to find all matches
       * below the cutoff, so we must find a previously cached query guaranteed to
       * have ALL matches to the current query under the given cutoff. */
      if ((maxSet && (result->getRMSDCutoff() - r > 0)) || (!maxSet && (result->getRMSDCutoff() - r >= cut))) {
        // If max is set, we want as large of a safe radius as possible. If max is
        // not set, then all suitable queries are safe, so we want as few extra
        // fluff to search through as possible
        mstreal curSafeRadius = result->getRMSDCutoff() - r;
        if ((bestComp == cache.end()) || ((maxSet && (bestDist < curSafeRadius)) || (!maxSet && (bestDist > curSafeRadius)))) {
          bestComp = it;
          safeRadius = curSafeRadius;
          bestDist = safeRadius; // could optimize in terms of things other than safe radius
        }
      }
    }
  }

  // first try going through matches of a close query
  if (bestComp != cache.end()) {
    if (isVerbose()) begin = chrono::high_resolution_clock::now();
    // visit all matches of the most suitable cached result
    fasstSolutionSet sols((*bestComp)->getSolutions(), topo);
    if (redSet) addSequenceContext(sols);
    vector<mstreal> rmsds = matchRMSDs(sols, queryAtoms);
    for (int k = 0; k < sols.size(); k++) {
      // If max was set and i already found maxN solutions, i only care to insert
      // those solutions that are better than the last solution found. NOTE: insert
      // can incur the cost of redundancy removal with respect to all previous
      // matches, so can be a costly operation
      bool foundMax = (maxSet && (matches.size() >= maxN));
      if ((foundMax && (rmsds[k] < matches.worstRMSD())) || (!foundMax && (rmsds[k] <= cut))) {
        rmsd = sols[k].getRMSD(); sols[k].setRMSD(rmsds[k]);      // overwrite with RMSD relative to current query
        if (redPropSet) {
          matches.insert(sols[k], getRedundancyPropertyMap());
        } else if (redSet) {
          matches.insert(sols[k], getRedundancyCut());  // insert copies the solution
        } else {
          matches.insert(sols[k]);
        }
        if (maxSet && (matches.size() > maxN)) matches.erase(--matches.end());
        sols[k].setRMSD(rmsd);                                    // set RMSD back to the old value
      }
    }
    if (isVerbose()) {
      end = chrono::high_resolution_clock::now();
      searchTime = chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
      cout << "\tquick-search time " << searchTime << " ms to go over " << sols.size() << " solutions" << endl;
    }
  }

  /* We can skip doing an actual search (i.e., we are guaranteed to already have
   * all of the relevant matches) if all of the following are true:
   * 1. a suitable neighboring cached query was available
   * 2. one of the following (mutually exclusive) conditions is true:
   *    A. max was not set, in which case, the earlier requirement that the safe
   *       RMSD be no less than the desired cutoff means we have found all
   *       matches period.
   *    B. max was set ADN either i) we found that many matches, all within the
   *       safe RMSD cutoff (i.e., they are guaranteed to be the best matches
   *       overall) or ii) we possibly found fewer matches BUT the safe radius
   *       was no smaller than the cutoff, so these are all the matches that
   *       there are in the full database (under the cutoff).
   * */
  bool noNeedForNewSearch = (bestComp != cache.end()) &&
                            (!maxSet || ((safeRadius >= cut) || ((matches.size() >= maxN) && (matches.worstRMSD() < safeRadius))));
  if (!noNeedForNewSearch) {
    // if there was a suitable neighbor cached, but the list of matches was not
    // deep enough, put some pressure on tollerance parameters
    if (bestComp != cache.end()) {
      if ((*bestComp)->isLimitedByMaxNumMatches()) { incMaxNumPressure(); }
      else { incErrTolPressure(); }
    }
    if (isVerbose()) {
      cout << "\t\tFAILED, need to search (" << matches.size() << " matches were found)";
      if (maxSet) cout << " (" << maxN << " was the max) ";
      if (matches.size() > 0) cout << " (worst RMSD was " << matches.worstRMSD() << ", cutoff was " << cut << ", and safeRadius was " << safeRadius << ", search params: " << (*bestComp)->getSearchRMSDCutoff() << " / " << (*bestComp)->getSearchMaxNumMatches() << ")";
      cout << endl << "\tnew pressures/tollerance factors for maxN and RMSD are: " << getMaxNumPressure() << " and " << getErrTolPressure() << " / " << maxNumFactor() << " and " << errTolFactor() << endl;
      begin = chrono::high_resolution_clock::now();
    }
    // -- loosen search criteria a bit to extract maximal value from search
    if (maxSet) {
      setMaxNumMatches(MstUtils::max(int(maxN*maxNumFactor()), maxN + 1000));
    }
    options().unsetRedundancyCut();
    options().unsetRedundancyProperty();
    setRMSDCutoff(cut*errTolFactor());

    // -- perform the search and cache
    matches = ((FASST*) this)->search();
    if (modifyPermission()) {
      if (matches.size() > 0) {
        cachedResult* result = new cachedResult(queryAtoms, matches, getRMSDCutoff(), getMaxNumMatches(), topo);
        cache.insert(result);
        if (isVerbose()) {
          cout << "\t\tfound " << matches.size() << " matches, last RMSD " << matches.rbegin()->getRMSD() << ", cutoff was " << getRMSDCutoff() << endl;
          cout << "\t\tcache now has " << cache.size() << " elements" << endl;
        }
        if (cache.size() > maxNumResults) {
          auto leastUseful = --cache.end();
          if (isVerbose()) cout << "\t\t\t\tERASING entry with priority " << (*leastUseful)->getPriority() << endl;
          delete(*leastUseful);
          cache.erase(leastUseful); // bump off the least used cached result if reached limit
        }
      }
    }

    // apply all needed cutoffs
    fasstSolutionSet finalMatches;
    vector<fasstSolution*> reordered;
    if (strictEquivalence()) reordered = matches.orderByDiscovery();
    if (redPropSet) options().setRedundancyProperty(redProp);
    if (redSet) addSequenceContext(matches);
    for (int i = 0; i < matches.size(); i++) {
      fasstSolution& sol = strictEquivalence() ? *(reordered[i]) : matches[i];
      if (redPropSet) {
        finalMatches.insert(sol, getRedundancyPropertyMap());
      } else if (redSet) {
        finalMatches.insert(sol, redCut);
      } else {
        finalMatches.insert(sol);
      }
      // since these are matches found specifically by searching for this qiery,
      // they are ordered by decreasing RMSD and if we already have maxN, we are done
      if (maxSet && (finalMatches.size() == maxN)) break;
    }
    matches = finalMatches;

    // -- reset search setting to their old values
    setOptions(origOpts);
    if (isVerbose()) {
      end = chrono::high_resolution_clock::now();
      searchTime = chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
      cout << "\t\tregular-search time " << searchTime << " ms" << endl;
    }
  } else {
    // if there was a suitable neighbor cached and the list of matches was deep
    // enough, reduce pressure on the tollerance parameters: redcuce more on the
    // inactive parameter and less on the active one.
    if ((*bestComp)->isLimitedByMaxNumMatches()) {
      decMaxNumPressure(0.01);
      decErrTolPressure(0.02);
    } else {
      decMaxNumPressure(0.02);
      decErrTolPressure(0.01);
    }

    // up the priority of this cached result just used
    cachedResult* result = *bestComp;
    cache.erase(bestComp);
    result->upPriority();
    cache.insert(result);
    if (isVerbose()) {
      cout << "\tdone upping priority" << std::endl;
      cout << "\tSUCCEEDED, NO need to search!!!" << endl;
      if (matches.size() > 0) cout << "\t\tfound " << matches.size() << " matches, worst RMSD was " << matches.worstRMSD() << ", cutoff was " << cut << ", and safeRadius was " << safeRadius << endl;
      cout << endl << "\tnew tollerance factors for maxN and RMSD are: " << maxNumFactor() << " and " << errTolFactor() << endl;
    }
    // update matches to reflect spacial alignment onto the new query
    matchRMSDs(matches, queryAtoms, true);
  }

  return matches;
}


/* --------- cFASST::cachedResult --------- */
cFASST::cachedResult::cachedResult(const AtomPointerVector& q, const fasstSolutionSet& sols, mstreal cut, int max, vector<int> topo) {
  q.clone(query);
  solSet = sols.extractAddresses();
  priority = 1.0; topology = topo;
  searchRMSDcut = cut; searchMaxNumMatches = max;
  rmsdCut = sols.rbegin()->getRMSD();
}

cFASST::cachedResult::cachedResult(const cachedResult& r) {
  r.query.clone(query);
  solSet = r.solSet;
  priority = r.priority; topology = r.topology;
  searchRMSDcut = r.searchRMSDcut; searchMaxNumMatches = r.searchMaxNumMatches;
  rmsdCut = r.rmsdCut;
}

cFASST::cachedResult::~cachedResult() {
  query.deletePointers();
}

bool cFASST::cachedResult::isSameTopology(const vector<int>& compTopo) const {
  if (topology.size() != compTopo.size()) return false;
  for (int i = 0; i < topology.size(); i++) {
    if (topology[i] != compTopo[i]) return false;
  }
  return true;
}

void cFASST::cachedResult::write(ostream &_os) const {
  MstUtils::writeBin(_os, topology);
  MstUtils::writeBin(_os, priority);
  MstUtils::writeBin(_os, searchRMSDcut);
  MstUtils::writeBin(_os, rmsdCut);
  MstUtils::writeBin(_os, searchMaxNumMatches);
  query.write(_os);
  MstUtils::writeBin(_os, (int) solSet.size());
  for (int i = 0; i < solSet.size(); i++) solSet[i].write(_os);
  // solSet->write(_os);
}

void cFASST::cachedResult::read(istream &_is) {
  // MstUtils::readBin(_is, topology);
  // MstUtils::readBin(_is, priority);
  // MstUtils::readBin(_is, searchRMSDcut);
  // MstUtils::readBin(_is, searchMaxNumMatches);
  // if (solSet != NULL) delete(solSet);
  // solSet = new fasstSolutionSet();
  // solSet->read(_is);
  // query.deletePointers();
  // query.read(_is);

  MstUtils::readBin(_is, topology);
  MstUtils::readBin(_is, priority);
  MstUtils::readBin(_is, searchRMSDcut);
  MstUtils::readBin(_is, rmsdCut);
  MstUtils::readBin(_is, searchMaxNumMatches);
  query.deletePointers();
  query.read(_is);
  int numSols; MstUtils::readBin(_is, numSols);
  solSet.clear(); solSet.resize(numSols);
  for (int i = 0; i < numSols; i++) solSet[i].read(_is);
}
