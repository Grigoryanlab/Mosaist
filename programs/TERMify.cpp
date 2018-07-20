#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>
#include <typeinfo>
#include <stdio.h>

#include "msttypes.h"
#include "mstfasst.h"
#include "mstcondeg.h"
#include "mstfuser.h"
#include "mstcondeg.h"
#include "mstoptions.h"
#include "mstmagic.h"

using namespace std;
using namespace MST;

class fasstCache {
  public:
    fasstCache(FASST* s) { S = s; maxNumResults = 1000; }
    ~fasstCache() { clear(); }
    void clear(); // clears cache (removes all solutions)
    FASST* getFASST() const { return S; }
    fasstSolutionSet search(bool verb = false);
    int getMaxNumResults() const { return maxNumResults; }

    void write(const string& filename) const {
      fstream out;
      MstUtils::openFile(out, filename, fstream::out | fstream::binary, "fasstCache::write");
      write(out);
      out.close();
    }

    void write(ostream &_os) const { // write the object to a binary stream
      MstUtils::writeBin(_os, maxNumResults);
      MstUtils::writeBin(_os, (int) cache.size());
      for (auto it = cache.begin(); it != cache.end(); ++it) (*it)->write(_os);
    }

    void read(const string& filename) {
      fstream in;
      MstUtils::openFile(in, filename, fstream::in | fstream::binary, "fasstCache::read");
      read(in);
      in.close();
    }

    void read(istream &_is) { // write object from a binary stream
      clear();
      MstUtils::readBin(_is, maxNumResults);
      int numResults; MstUtils::readBin(_is, numResults);
      for (int i = 0; i < numResults; i++) {
        fasstCachedResult* result = new fasstCachedResult();
        result->read(_is);
        cache.insert(result);
      }
    }

  protected:
    static vector<int> getStructureTopology(const Structure& S) {
      vector<int> topo(S.chainSize(), 0);
      for (int i = 0; i < S.chainSize(); i++) topo[i] = S[i].residueSize();
      return topo;
    }

  private:
    class fasstCachedResult {
      friend class fasstCache;
      public:
        fasstCachedResult() { rmsdCut = 0.0; solSet = NULL; priority = 1.0; }
        fasstCachedResult(const AtomPointerVector& q, const fasstSolutionSet& sols, mstreal cut, vector<int> topo) {
          q.clone(query);
          solSet = new fasstSolutionSet(sols);
          priority = 1.0; topology = topo; rmsdCut = cut;
        }
        fasstCachedResult(const fasstCachedResult& r) {
          r.query.clone(query);
          solSet = new fasstSolutionSet(*(r.solSet));
          rmsdCut = r.rmsdCut; priority = r.priority; topology = r.topology;
        }
        ~fasstCachedResult() {
          query.deletePointers();
          delete(solSet);
        }
        void upPriority(mstreal del = 1.0) { priority += del; }
        void elapsePriority(int Thalf) { priority *= pow(0.5, 1.0/Thalf); }

        AtomPointerVector getQuery() const { return query; }
        fasstSolutionSet* getSolutions() const { return solSet; }
        mstreal getRMSDCutoff() const { return rmsdCut; }
        mstreal getPriority() const { return priority; }
        vector<int> getTopology() const { return topology; }
        bool isSameTopology(const vector<int>& compTopo) const {
          if (topology.size() != compTopo.size()) return false;
          for (int i = 0; i < topology.size(); i++) {
            if (topology[i] != compTopo[i]) return false;
          }
          return true;
        }

        friend bool operator<(const fasstCachedResult& ri, const fasstCachedResult& rj) {
          if (ri.priority != rj.priority) return (ri.priority > rj.priority);
          if (ri.solSet->size() != rj.solSet->size()) return ri.solSet->size() < rj.solSet->size();
          if (ri.query.size() != rj.query.size()) return ri.query.size() < rj.query.size();
          return &ri < &rj;
        }

        friend ostream& operator<<(ostream &_os, const fasstCachedResult& _res) {
          _os << "[" << MstUtils::vecToString(_res.topology) << "] " << _res.rmsdCut << "(" << _res.priority << ")" << endl;
          _os << _res.query << endl << *(_res.solSet);
          return _os;
        }

        void write(ostream &_os) const {
          MstUtils::writeBin(_os, topology);
          MstUtils::writeBin(_os, rmsdCut);
          MstUtils::writeBin(_os, priority);
          solSet->write(_os);
          query.write(_os);
        }

        void read(istream &_is) {
          MstUtils::readBin(_is, topology);
          MstUtils::readBin(_is, rmsdCut);
          MstUtils::readBin(_is, priority);
          if (solSet != NULL) delete(solSet);
          solSet = new fasstSolutionSet();
          solSet->read(_is);
          query.deletePointers();
          query.read(_is);
        }

      private:
        AtomPointerVector query;
        vector<int> topology;
        mstreal rmsdCut;
        fasstSolutionSet* solSet;
        mstreal priority;
    };

    struct compResults {
      bool operator() (const fasstCachedResult* lhs, const fasstCachedResult* rhs) const {
        return (*lhs < *rhs);
      }
    };

    int maxNumResults; // max number of searches to cache
    FASST* S;
    set<fasstCachedResult*, compResults> cache;
    RMSDCalculator rc;
};

void fasstCache::clear() {
  for (auto it = cache.begin(); it != cache.end(); ++it) delete(*it);
  cache.clear();
}

fasstSolutionSet fasstCache::search(bool verb) {
  int n;
  if (S->isSufficientNumMatchesSet()) MstUtils::error("cannot cache when \"sufficient\" number of matches is set", "fasstCache::search");
  fasstSolutionSet matches;
  mstreal cut = S->getRMSDCutoff(), rmsd;
  bool maxSet = S->isMaxNumMatchesSet();
  bool redSet = S->isRedundancyCutSet();
  bool minSet = S->isMinNumMatchesSet();
  int maxN = S->getMaxNumMatches();
  int minN = S->getMinNumMatches();
  mstreal redCut = S->getRedundancyCut();
  // for some of the logic below it is convenient to assume no funny business
  if (minSet && (minN <= 0)) MstUtils::error("min number of matches is set, but actual limit is not positive");
  if (maxSet && (maxN <= 0)) MstUtils::error("min number of matches is set, but actual limit is not positive");
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
    if (verb) cout << " " << res->getPriority();
  }
  if (verb) cout << endl;

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
      // take all below the cutoff, but at least minN and at most maxN (if set)
      if ((rmsds[k] <= cut) || (minSet && (matches.size() < minN))) {
        rmsd = sols[k].getRMSD(); sols[k].setRMSD(rmsds[k]);      // overwrite with RMSD relative to current query
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
   * 3. one of the following (mutually exclusive) conditions is true:
   *    A. min was not set
   *    B. min was set and we have at least that many matches, all of which are
   *       under the safe RMSD cutoff.
   * */
  bool noNeedForNewSearch = (bestComp != cache.end()) &&
                            (!maxSet || ((safeRadius >= cut) || ((matches.size() >= maxN) && (matches.worstRMSD() < safeRadius)))) &&
                            (!minSet || ((matches.size() >= minN) && (matches.worstRMSD() < safeRadius)));
  bool doNewSearch = !noNeedForNewSearch;
  if (doNewSearch) {
    if (verb) {
      if (verb) cout << "\t\tFAILED, need to search (" << matches.size() << " matches were found)...";
      if (matches.size() > 0) cout << " (worst RMSD was " << matches.worstRMSD() << ", cutoff was " << cut << ", and safeRadius was " << safeRadius << ")";
      cout << endl;
      begin = chrono::high_resolution_clock::now();
    }
    // -- loosen search criteria a bit to extract maximal value from search
    if (maxSet) {
      S->setMaxNumMatches(MstUtils::max(int(maxN*2 + 1), maxN + 1000));
      // S->unsetMaxNumMatches();
    }
    // set redundancy value to above 1 in order to accumulate sequence context,
    // regardless of current requirements (for any future needs)
    S->setRedundancyCut(1.1);
    S->setRMSDCutoff(1.1*cut);

    /* The minimum number of matches limit is costly to include in the search,
     * but is sometimes not necesary (i.e., you find at least that many matches
     * under the given RMSD cutoff anyway). So first try without it. */
    for (int c = 0; c < (minSet ? 2 : 1); c++) {
      if (minSet) {
        // if removing the min requirement does not fly, bring the requirement
        // back, by find a few more. In this way, if this type of search is
        // repeated with a slightly different query, we will be more likely to
        // succeed with a lookup
        if (c) S->setMinNumMatches(minN*2);
        else S->unsetMinNumMatches();
        if (c && verb) cout << "\t\t\tDID NOT WORK without min (found only " << matches.size() << " matches, while minN = " << minN << "), so repeating search..." << endl;
      }

      // -- perform the search and cache
      matches = S->search();
      if (matches.size() > 0) {
        fasstCachedResult* result = new fasstCachedResult(queryAtoms, matches, matches.rbegin()->getRMSD(), topo);
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
        bool haveMin = !minSet || (minSet && (finalMatches.size() >= minN));
        if (haveMin && (sol.getRMSD() > cut)) continue;
        if (redSet) finalMatches.insert(sol, redCut);
        else finalMatches.insert(sol);
        if (maxSet && (finalMatches.size() > maxN)) finalMatches.erase(--finalMatches.end());
      }
      matches = finalMatches;

      // did we get the min number of matches anyway?
      if (minSet && (c == 0) && (matches.size() >= minN)) {
        break;
      }
    }

    // -- reset search setting to their old values
    S->setRedundancyCut(redCut);
    S->setRMSDCutoff(cut);
    if (maxSet) S->setMaxNumMatches(maxN);
    if (minSet) S->setMinNumMatches(minN);
    if (verb) {
      end = chrono::high_resolution_clock::now();
      searchTime = chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
      cout << "\t\tregular-search time " << searchTime << " ms" << endl;
    }
  } else {
    // up the priority of this cached result just used
    fasstCachedResult* result = *bestComp;
    cache.erase(bestComp);
    result->upPriority();
    cache.insert(result);
    if (verb) {
      cout << "\tdone upping priority" << std::endl;
      cout << "\tSUCCEEDED, NO need to search!!!" << endl;
      if (matches.size() > 0) cout << "\tworst RMSD was " << matches.worstRMSD() << ", cutoff was " << cut << ", and safeRadius was " << safeRadius << endl;
    }
    // update matches to reflect spacial alignment onto the new query
    S->matchRMSDs(matches, queryAtoms, true);
  }

  return matches;
}


AtomPointerVector getBackbone(const Structure& S, vector<int> residues) {
  AtomPointerVector atoms;
  vector<string> bba = {"N", "CA", "C", "O"};
  for (int i = 0; i < residues.size(); i++) {
    Residue& res = S.getResidue(residues[i]);
    for (int j = 0; j < bba.size(); j++) {
      atoms.push_back(res.findAtom(bba[j]));
    }
  }
  return atoms;
}

mstreal getRadius(const Structure& S) {
  selector sel(S);
  AtomPointerVector atoms = sel.select("name N or name CA or name C or name O");
  mstreal rad = 0;
  for (int i = 0; i < atoms.size(); i++) {
    for (int j = i+1; j < atoms.size(); j++) {
      mstreal d = atoms[i]->distance(atoms[j]);
      if (d > rad) rad = d;
    }
  }
  return rad;
}

vector<Structure*> getMatches(fasstCache& C, Structure& frag, vector<int>& fragResIdx) {
  FASST* F = C.getFASST();
  F->setQuery(frag, false);
  // F->setRMSDCutoff(RMSDCalculator::rmsdCutoff(frag));
  F->setMaxNumMatches(200);
  F->setMinNumMatches(100);
// TODO: fix so that it _always_ reports the same results as regular searching (basically, up to removing duplicates in the right way)
if (0) {
  fasstSolutionSet matchesOld = F->search();
  fasstSolutionSet matchesNew = C.search(true);
  bool error = false;
  if (matchesOld.size() != matchesNew.size()) { error = true; }
  else {
    for (int i = 0; i < matchesOld.size(); i++) {
      if (((matchesOld[i] < matchesNew[i]) || (matchesNew[i] < matchesOld[i])) && (fabs(matchesOld[i].getRMSD()-matchesNew[i].getRMSD()) > 10E-8)) {
        cout << matchesOld[i] << endl << "vs" << endl << matchesNew[i] << endl;
        printf("RMSD difference = %e\n", matchesOld[i].getRMSD() - matchesNew[i].getRMSD());
        error = true; break;
      }
    }
  }
  if (error) {
    cout << "cutoff was " << F->getRMSDCutoff() << "\nold style:\n" << matchesOld << "\nnew style:\n" << matchesNew << endl;
    exit(1);
  }
}

  fasstSolutionSet matches = C.search(false);
  cout << "\tfound " << matches.size() << " matches" << endl;
  vector<Structure*> matchStructures;
  for (auto it = matches.begin(); it != matches.end(); ++it) {
    matchStructures.push_back(new Structure(F->getMatchStructure(*it, false, FASST::matchType::REGION)));
    Structure& match = *(matchStructures.back());
    MstUtils::assert(match.residueSize() == fragResIdx.size(), "unexpected match size");
    for (int k = 0; k < match.residueSize(); k++) {
      match.getResidue(k).setNum(fragResIdx[k]); // make residue numbers store indices into the original structure
    }
  }
  return matchStructures;
}

void addMatches(vector<Structure*>& matchStructures, vector<vector<Residue*> >& resTopo, fstream& matchOut, int ri = -1) {
  for (int i = 0; i < matchStructures.size(); i++) {
    if ((ri >= 0) && (i != ri)) continue;
    Structure& match = *(matchStructures[i]);
    if (matchOut.is_open()) {
      match.writePDB(matchOut);
      matchOut << "END" << endl;
    }
    for (int k = 0; k < match.residueSize(); k++) {
      Residue& res = match.getResidue(k);
      resTopo[res.getNum()].push_back(&res); // residue numbers store indices into the originating structure
    }
  }
}

bool mc(mstreal oldScore, mstreal newScore, mstreal kT) {
  return ((newScore < oldScore) || (MstUtils::randUnit() < exp((oldScore - newScore)/kT)));
}

fusionTopology getTopo(int L, vector<vector<Structure*> >& allMatches, vector<int>& picks, fstream& matchOut, Structure* global = NULL) {
  fusionTopology resTopo(L);
  for (int si = 0; si < allMatches.size(); si++) {
    Structure& match = *(allMatches[si][picks[si]]);
    resTopo.addFragment(match);
    if (matchOut.is_open()) { match.writePDB(matchOut); matchOut << "END" << endl; }
  }
  if (global != NULL) {
    MstUtils::assert(global->residueSize() == L, "the global target structure specified has an unexpected number of residues for the topology");
    resTopo.addFragment(*global, MstUtils::range(0, L), -10.0);
    if (matchOut.is_open()) { global->writePDB(matchOut); matchOut << "END" << endl; }
  }
  return resTopo;
}

fusionTopology getTopo(int L, vector<vector<Structure*> >& allMatches, vector<vector<int> >& picks, fstream& matchOut, Structure* global = NULL) {
  fusionTopology resTopo(L);
  for (int si = 0; si < allMatches.size(); si++) {
    for (int j = 0; j < picks[si].size(); j++) {
      Structure& match = *(allMatches[si][picks[si][j]]);
      resTopo.addFragment(match);
      if (matchOut.is_open()) { match.writePDB(matchOut); matchOut << "END" << endl; }
    }
  }
  if (global != NULL) {
    MstUtils::assert(global->residueSize() == L, "the global target structure specified has an unexpected number of residues for the topology");
    resTopo.addFragment(*global, MstUtils::range(0, L), -10.0);
    if (matchOut.is_open()) { global->writePDB(matchOut); matchOut << "END" << endl; }
  }
  return resTopo;
}

AtomPointerVector getCorrespondingAtoms(Structure& from, Structure& like) {
  MstUtils::assert(from.residueSize() == like.residueSize(), "the two structures must have the same number of residues", "getCorrespondingAtoms()");
  AtomPointerVector atoms;
  for (int ri = 0; ri < like.residueSize(); ri++) {
    Residue& fromRes = from.getResidue(ri);
    Residue& likeRes = like.getResidue(ri);
    for (int ai = 0; ai < likeRes.atomSize(); ai++) {
      atoms.push_back(fromRes.findAtom(likeRes[ai].getName()));
    }
  }
  return atoms;
}

mstreal totalScore(fusionScores& scoreObj, Structure& fused, AtomPointerVector& init, bool report = false) {
  RMSDCalculator rc;
  // mstreal a = scoreObj.getScore()/10;
  // mstreal a = scoreObj.getTotRMSDScore();
  // mstreal a = scoreObj.getScore();
  // mstreal b = rc.bestRMSD(init, fused.getAtoms());
  // if (report) cout << "fuser: " << scoreObj << "; total = " << a << " - " << b;
  if (report) cout << "fuser: " << scoreObj << "; RMSD from start = " << rc.bestRMSD(init, fused.getAtoms());
  // return a - b;
  return scoreObj.getScore();
}

int main(int argc, char** argv) {
  // TODO: debug Neilder-Meid optimization by comparing a simple case with Matlab
  // TODO: enable a setting in Fuser, whereby fully overlapping segments are scored
  //       (in terms of RMSD) via a weighted average, such that the lowest-RMSD
  //       segment makes a dominant contribution to the score
  MstOptions op;
  op.setTitle("Starting from some arbitrary conformation of a chain, iteratively build a TERM-based compact structure. Options:");
  op.addOption("p", "starting conformation PDB file.", true);
  op.addOption("rLib", "a path to an MST-formatter rotamer library.", true);
  op.addOption("d", "a database file with a list of PDB files.");
  op.addOption("b", "a binary database file. If both --d and --b are given, will overwrite this file with a corresponding binary database.");
  op.addOption("o", "output base name.", true);
  op.addOption("n", "pick this many matches for each TERM for fusion (default is 1). At each iteration, one match for a randomly-selected TERM will be randomly substituted.");
  op.addOption("r", "if specified, will randomly pick the matches for the first iteration. By default, will take the top match for each TERM or the top --n matches, if --n is specified.");
  op.addOption("m", "if specified, will try to move away from the original structure by including a term in the objective function that rewards large RMSDs to it. Default is no.");
  op.addOption("cyc", "number of cycles--i.e., number of times fresh TERMs are searched for (10 by default).");
  op.addOption("iter", "number of iterations per cycle (1 by default). At the start of each iteration, the overall structure is reinitialized to the current structure");
  op.addOption("f", "a quoted, space-separated list of 0-initiated residue integers to fix.");
  op.addOption("fs", "a selection string for residues to fix.");
  op.addOption("rad", "compactness radius. Default will be based on protein length.");
  if (op.isGiven("f") && op.isGiven("fs")) MstUtils::error("only one of --f or --fs can be given!");
  MstUtils::setSignalHandlers();
  op.setOptions(argc, argv);
  RMSDCalculator rc;
  Structure I(op.getString("p"));
  FASST F;
  fasstCache cache(&F);
  vector<int> fixed;
  F.setMemorySaveMode(true);
  if (op.isGiven("d")) {
    F.addTargets(MstUtils::fileToArray(op.getString("d")));
    if (op.isGiven("b")) {
      F.writeDatabase(op.getString("b"));
    }
  } else if (op.isGiven("b")) {
    F.readDatabase(op.getString("b"));
  } else {
    MstUtils::error("either --b or --d must be given!");
  }
  if (op.isGiven("f")) {
    fixed = MstUtils::splitToInt(op.getString("f"));
  }
  if (op.isGiven("fs")) {
    selector sel(I);
    vector<Residue*> fixedResidues = sel.selectRes(op.getString("fs"));
    cout << "fix selection gave " << fixedResidues.size() << " residues, fixing..." << endl;
    map<Residue*, int> indices = MstUtils::indexMap(I.getResidues());
    fixed.resize(fixedResidues.size());
    for (int i = 0; i < fixedResidues.size(); i++) fixed[i] = indices[fixedResidues[i]];
  }
  int numPerTERM = op.getInt("n", 1);

  F.setRedundancyCut(0.5);
  RotamerLibrary RL(op.getString("rLib"));
  int pmSelf = 2, pmPair = 1;
  int Ni = 1000;
  fusionScores bestScore, currScore;
  fstream out, shellOut, dummy;

  // TERMify loop
  Structure S = I.reassignChainsByConnectivity(); // fixed residues are already selected, but this does not change residue order
  RotamerLibrary::standardizeBackboneNames(S);
  mstreal R0 = getRadius(I);
  mstreal Rf = op.getReal("rad", pow(I.residueSize() * 1.0, 1.0/3)*5.0);
  int Ncyc = op.getInt("cyc", 10);
  MstUtils::openFile(out, op.getString("o") + ".traj.pdb", ios::out);
  out << "MODEL " << 0 << endl; S.writePDB(out); out << "ENDMDL" << endl;
  for (int c = 0; c < Ncyc; c++) {
    cout << "Cycle " << c+1 << "..." << endl;
    if (c == 0) {
      MstUtils::openFile(shellOut, op.getString("o") + ".shell.init.pdb", ios::out);
    } else if (c == Ncyc - 1) {
      MstUtils::openFile(shellOut, op.getString("o") + ".shell.pdb", ios::out);
    }

    /* --- Decorate the current conformation with TERMs --- */
    // first self TERMs
    cout << "Searching for self TERMs..." << endl;
    vector<vector<Structure*>> allMatches;
    for (int ci = 0; ci < S.chainSize(); ci++) {
      Chain& C = S[ci];
      for (int ri = 0; ri < C.residueSize(); ri++) {
        Structure frag; vector<int> fragResIdx;
        TERMUtils::selectTERM({&C[ri]}, frag, pmSelf, &fragResIdx);
        if (MstUtils::setdiff(fragResIdx, fixed).empty()) continue; // TERMs composed entirely of fixed residues have no impact
        cout << "TERM around " << C[ri] << endl;
        F.setRMSDCutoff(RMSDCalculator::rmsdCutoff(fragResIdx, S)); // account for spacing between residues from the same chain
        vector<Structure*> matches = getMatches(cache, frag, fragResIdx);
        allMatches.push_back(matches);
      }
    }

    // then pair TERMs
    cout << "Searching for pair TERMs..." << endl;
    ConFind cfd(&RL, S);
    contactList L = cfd.getContacts(S, 0.01);
    vector<pair<Residue*, Residue*> > contactList = L.getOrderedContacts();
    for (int k = 0; k < contactList.size(); k++) {
      Residue* resA = contactList[k].first;
      Residue* resB = contactList[k].second;
      Structure frag; vector<int> fragResIdx;
      TERMUtils::selectTERM({resA, resB}, frag, pmPair, &fragResIdx);
      if (MstUtils::setdiff(fragResIdx, fixed).empty()) continue; // TERMs composed entirely of fixed residues have no impact
      cout << "TERM around " << *resA << " x " << *resB << endl;
      F.setRMSDCutoff(RMSDCalculator::rmsdCutoff(fragResIdx, S)); // account for spacing between residues from the same chain
      vector<Structure*> matches = getMatches(cache, frag, fragResIdx);
      allMatches.push_back(matches);
    }

    // fuser options
    fusionParams opts; opts.setNumIters(Ni); opts.setVerbose(false);
    opts.setMinimizerType(fusionParams::gradDescent);
    opts.setRepFC(1);
    opts.setCompFC(0.1);
    mstreal compactnessRadius = Rf;
    // (R0 - Rf)*exp(-c/10.0) + Rf; // exponential scaling
    // (Rf*(c + 1) + R0*(Ncyc - c - 1))/Ncyc; // linear scaling
    opts.setCompRad(compactnessRadius);
    cout << "will be trying to combine structure to a radius of " << compactnessRadius << endl;

    /* --- pick a random combination of TERMs and fuse --- */
    vector<vector<int> > currPicks(allMatches.size()), bestPicks;
    vector<vector<Residue*> > resTopo(I.residueSize());
    for (int si = 0; si < allMatches.size(); si++) {
      for (int j = 0; j < MstUtils::min(numPerTERM, (int) allMatches[si].size()); j++) {
        if (op.isGiven("r")) {
          currPicks[si].push_back(MstUtils::randInt(allMatches[si].size()));
        } else {
          currPicks[si].push_back(j);
        }
      }
    }

    // if there are fixed residues, then their starting conformation (and thus)
    // the final one also) is chosen from the original structure
    if (!fixed.empty()) opts.setStartingStructure(S);

    /* --- do an MC simulation to find a good combo of TERMs --- */
    fusionTopology bestTopo, currTopo;
    Structure bestFused, currFused;
    fusionScores bestScore, currScore, propScore;
    mstreal kT = 0.001;
    for (int it = 0; it < op.getInt("iter", 1); it++) {
      vector<vector<int> > propPicks = currPicks;
      // make a "mutation"
      if (it != 0) {
        int si = MstUtils::randInt(allMatches.size());
        int mi = MstUtils::randInt(propPicks[si].size());
        propPicks[si][mi] = MstUtils::randInt(allMatches[si].size());
      }
      fusionTopology propTopo = getTopo(I.residueSize(), allMatches, propPicks, (it == 0) ? shellOut : dummy, op.isGiven("m") ? &S : NULL);
      propTopo.addFixedPositions(fixed);
      // this is to provide IC information for connections between any possible
      // fixed residues residues and flexible ones. NOTE: the weight is set to 0,
      // so this does not affect the objective function (not needed when there are
      // no fixed residues, but adding just in case)
      propTopo.addFragment(S, MstUtils::range(0, propTopo.length()), 0.0);
      Structure propFused = Fuser::fuse(propTopo, propScore, opts);
      AtomPointerVector init = getCorrespondingAtoms(S, propFused);
      cout << "\titeration " << it << " => "; totalScore(propScore, propFused, init, true); cout << endl;
      if ((it == 0) || mc(totalScore(currScore, currFused, init), totalScore(propScore, propFused, init), kT)) {
        cout << "\t\taccepted" << endl;
        currFused = propFused;
        currTopo = propTopo;
        currScore = propScore;
        currPicks = propPicks;
      }
      if ((it == 0) || (totalScore(propScore, propFused, init) < totalScore(bestScore, bestFused, init))) {
        cout << "\t\t\tnew best" << endl;
        bestFused = propFused;
        bestTopo = propTopo;
        bestScore = propScore;
        bestPicks = propPicks;
      }
    }

    out << "MODEL " << c+1 << endl;
    bestFused.writePDB(out);
    out << "ENDMDL" << endl;

    // align based on the fixed part, if anything was fixed (for ease of visualization)
    if (fixed.size() > 0) {
      AtomPointerVector before = getBackbone(S, fixed);
      AtomPointerVector after = getBackbone(bestFused, fixed);
      rc.align(after, before, bestFused);
    }

    /* --- write intermediate result and clean up--- */
    S = bestFused.reassignChainsByConnectivity();
    for (int si = 0; si < allMatches.size(); si++) {
      for (int mi = 0; mi < allMatches[si].size(); mi++) {
        delete(allMatches[si][mi]);
      }
    }
    if (shellOut.is_open()) shellOut.close();
  }
  out.close();
  S.writePDB(op.getString("o") + ".fin.pdb");
}
