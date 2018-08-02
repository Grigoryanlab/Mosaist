#include "mstfasstcache.h"

/* --------- fasstCache --------- */
void fasstCache::write(const string& filename) const {
  fstream out;
  MstUtils::openFile(out, filename, fstream::out | fstream::binary, "fasstCache::write");
  write(out);
  out.close();
}

void fasstCache::write(ostream &_os) const {
  MstUtils::writeBin(_os, maxNumResults);
  MstUtils::writeBin(_os, (int) cache.size());
  for (auto it = cache.begin(); it != cache.end(); ++it) (*it)->write(_os);
}

void fasstCache::read(const string& filename) {
  fstream in;
  MstUtils::openFile(in, filename, fstream::in | fstream::binary, "fasstCache::read");
  read(in);
  in.close();
}

void fasstCache::read(istream &_is) { // read object from a binary stream
  clear();
  MstUtils::readBin(_is, maxNumResults);
  int numResults; MstUtils::readBin(_is, numResults);
  for (int i = 0; i < numResults; i++) {
    fasstCachedResult* result = new fasstCachedResult();
    result->read(_is);
    cache.insert(result);
  }
}

vector<int> fasstCache::getStructureTopology(const Structure& S) {
  vector<int> topo(S.chainSize(), 0);
  for (int i = 0; i < S.chainSize(); i++) topo[i] = S[i].residueSize();
  return topo;
}

void fasstCache::clear() {
  for (auto it = cache.begin(); it != cache.end(); ++it) delete(*it);
  cache.clear();
}

fasstSolutionSet fasstCache::search(bool verb) {
  int n;
  if (S->isSufficientNumMatchesSet()) MstUtils::error("cannot cache when \"sufficient\" number of matches is set", "fasstCache::search");
  if (S->isMinNumMatchesSet()) MstUtils::error("cannot cache when \"minimum\" number of matches is set", "fasstCache::search");
  fasstSolutionSet matches;
  mstreal cut = S->getRMSDCutoff(), rmsd;
  bool maxSet = S->isMaxNumMatchesSet();
  bool redSet = S->isRedundancyCutSet();
  int maxN = S->getMaxNumMatches();
  mstreal redCut = S->getRedundancyCut();
  chrono::high_resolution_clock::time_point begin, end;
  int searchTime;
  AtomPointerVector queryAtoms = S->getQuerySearchedAtoms();

  /* Old cached results should eventually "expire", so uniformly lower priority
   * slightly first. This way, cached results that have not been used in a while
   * will eventually have a lower priority than brand new searches and these
   * will then push out these old (aparently) useless results. */
  for (auto it = cache.begin(); it != cache.end(); it++) {
    // NOTE: gets rid of const qualifier! This is safe to do only because I know
    // I will monotonically lower everybody's priority (order will not change).
    fasstCachedResult* res = &(*(*it));
    res->elapsePriority(getMaxNumResults());
    // if (verb) cout << " " << res->getPriority();
  }
  // if (verb) cout << endl;

  /* See whether all matches for the current query, within the given cutoff, are
   * among the list of matches of some previously cached query. */
  vector<int> topo = fasstCache::getStructureTopology(S->getQuery());
  auto bestComp = cache.end(); mstreal bestDist = -1, safeRadius = -1;
  for (auto it = cache.begin(); it != cache.end(); ++it) {
    fasstCachedResult* result = *it;
    if (!result->isSameTopology(topo)) continue;
    mstreal r = rc.bestRMSD(queryAtoms, result->getQuery());
    /* If a maximum number of matches is set, then we _may_ not need to find ALL
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

  // first try going through matches of a close query
  if (bestComp != cache.end()) {
    if (verb) begin = chrono::high_resolution_clock::now();
    // visit all matches of the most suitable cached result
    fasstSolutionSet& sols = *((*bestComp)->getSolutions());
    vector<mstreal> rmsds = S->matchRMSDs(sols, queryAtoms);
    for (int k = 0; k < sols.size(); k++) {
      // take all below the cutoff
      if (rmsds[k] <= cut) {
        rmsd = sols[k].getRMSD(); sols[k].setRMSD(rmsds[k]);      // overwrite with RMSD relative to current query
// TODO: if match does not already context info and redundancy is set, add context info!
        if (redSet) matches.insert(sols[k], S->getRedundancyCut());  // insert copies the solution
        else matches.insert(sols[k]);
        if (maxSet && (matches.size() > maxN)) matches.erase(--matches.end());
        sols[k].setRMSD(rmsd);                                    // set RMSD back to the old value
      }
    }
    if (verb) {
      end = chrono::high_resolution_clock::now();
      searchTime = chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
      if (verb) cout << "\tquick-search time " << searchTime << " ms" << std::endl;
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
    if (verb) {
      cout << "\t\tFAILED, need to search (" << matches.size() << " matches were found)...";
      if (matches.size() > 0) cout << " (worst RMSD was " << matches.worstRMSD() << ", cutoff was " << cut << ", and safeRadius was " << safeRadius << ", search params: " << (*bestComp)->getSearchRMSDCutoff() << " / " << (*bestComp)->getSearchMaxNumMatches() << ")";
      cout << endl << "\tnew pressures/tollerance factors for maxN and RMSD are: " << getMaxNumPressure() << " and " << getErrTolPressure() << " / " << maxNumFactor() << " and " << errTolFactor() << endl;
      begin = chrono::high_resolution_clock::now();
    }
    // -- loosen search criteria a bit to extract maximal value from search
    if (maxSet) {
      S->setMaxNumMatches(MstUtils::max(int(maxN*maxNumFactor()), maxN + 1000));
      // S->unsetMaxNumMatches();
    }
    // set redundancy value to above 1 in order to accumulate sequence context,
    // regardless of current requirements (for any future needs)
    S->setRedundancyCut(1.1); // TODO: this is probably wasteful; if redundancy is set, we can decorate only the final matches (only when needed); if not, then forget about it.
    S->setRMSDCutoff(cut*errTolFactor());

    // -- perform the search and cache
    matches = S->search();
    if (matches.size() > 0) {
      fasstCachedResult* result = new fasstCachedResult(queryAtoms, matches, S->getRMSDCutoff(), S->getMaxNumMatches(), topo);
      cache.insert(result);
      if (verb) {
        cout << "\t\tfound " << matches.size() << " matches, last RMSD " << matches.rbegin()->getRMSD() << ", cutoff was " << S->getRMSDCutoff() << endl;
        cout << "\t\tcache now has " << cache.size() << " elements" << endl;
      }
      if (cache.size() > maxNumResults) {
        auto leastUseful = --cache.end();
        if (verb) cout << "\t\t\t\tERASING entry with priority " << (*leastUseful)->getPriority() << endl;
        delete(*leastUseful);
        cache.erase(leastUseful); // bump off the least used cached result if reached limit
      }
    }

    // apply all needed cutoffs
    fasstSolutionSet finalMatches;
    // vector<fasstSolution*> orderedMatches = matches.orderByDiscovery();
    for (int i = 0; i < matches.size(); i++) {
      const fasstSolution& sol = matches[i];
// TODO: if match does not already context info and redundancy is set, add context info!
      if (redSet) finalMatches.insert(sol, redCut);
      else finalMatches.insert(sol);
      if (maxSet && (finalMatches.size() > maxN)) finalMatches.erase(--finalMatches.end());
    }
    matches = finalMatches;

    // -- reset search setting to their old values
    S->setRedundancyCut(redCut);
    S->setRMSDCutoff(cut);
    if (maxSet) S->setMaxNumMatches(maxN);
    if (verb) {
      end = chrono::high_resolution_clock::now();
      searchTime = chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
      cout << "\t\tregular-search time " << searchTime << " ms" << endl;
    }
  } else {
    // if there was a suitable neighbor cached and the list of matches was deep
    // enough, reduce pressure on the tollerance parameters: redcuce more on the
    // inactive parameter and less on the active one.
    if ((*bestComp)->isLimitedByMaxNumMatches()) {
      decMaxNumPressure(0.1);
      decErrTolPressure(0.2);
    } else {
      decMaxNumPressure(0.2);
      decErrTolPressure(0.1);
    }

    // up the priority of this cached result just used
    fasstCachedResult* result = *bestComp;
    cache.erase(bestComp);
    result->upPriority();
    cache.insert(result);
    if (verb) {
      cout << "\tdone upping priority" << std::endl;
      cout << "\tSUCCEEDED, NO need to search!!!" << endl;
      if (matches.size() > 0) cout << "\tworst RMSD was " << matches.worstRMSD() << ", cutoff was " << cut << ", and safeRadius was " << safeRadius << endl;
      cout << endl << "\tnew tollerance factors for maxN and RMSD are: " << maxNumFactor() << " and " << errTolFactor() << endl;
    }
    // update matches to reflect spacial alignment onto the new query
    S->matchRMSDs(matches, queryAtoms, true);
  }

  return matches;
}


/* --------- fasstCache::fasstCachedResult --------- */
fasstCache::fasstCachedResult::fasstCachedResult(const AtomPointerVector& q, const fasstSolutionSet& sols, mstreal cut, int max, vector<int> topo) {
  q.clone(query);
  solSet = new fasstSolutionSet(sols);
  priority = 1.0; topology = topo;
  searchRMSDcut = cut; searchMaxNumMatches = max;
}

fasstCache::fasstCachedResult::fasstCachedResult(const fasstCachedResult& r) {
  r.query.clone(query);
  solSet = new fasstSolutionSet(*(r.solSet));
  priority = r.priority; topology = r.topology;
  searchRMSDcut = r.searchRMSDcut; searchMaxNumMatches = r.searchMaxNumMatches;
}

fasstCache::fasstCachedResult::~fasstCachedResult() {
  query.deletePointers();
  delete(solSet);
}

bool fasstCache::fasstCachedResult::isSameTopology(const vector<int>& compTopo) const {
  if (topology.size() != compTopo.size()) return false;
  for (int i = 0; i < topology.size(); i++) {
    if (topology[i] != compTopo[i]) return false;
  }
  return true;
}

void fasstCache::fasstCachedResult::write(ostream &_os) const {
  MstUtils::writeBin(_os, topology);
  MstUtils::writeBin(_os, priority);
  MstUtils::writeBin(_os, searchRMSDcut);
  MstUtils::writeBin(_os, searchMaxNumMatches);
  solSet->write(_os);
  query.write(_os);
}

void fasstCache::fasstCachedResult::read(istream &_is) {
  MstUtils::readBin(_is, topology);
  MstUtils::readBin(_is, priority);
  MstUtils::readBin(_is, searchRMSDcut);
  MstUtils::readBin(_is, searchMaxNumMatches);
  if (solSet != NULL) delete(solSet);
  solSet = new fasstSolutionSet();
  solSet->read(_is);
  query.deletePointers();
  query.read(_is);
}
