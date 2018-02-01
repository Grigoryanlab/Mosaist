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

using namespace std;
using namespace MST;

class fasstCache {
  public:
    fasstCache(FASST* s) { S = s; maxNumResults = 1000; }
    ~fasstCache();
    FASST* getFASST() const { return S; }
    fasstSolutionSet search(bool verb = false);
    int getMaxNumResults() const { return maxNumResults; }

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

fasstCache::~fasstCache() {
  for (auto it = cache.begin(); it != cache.end(); ++it) delete(*it);
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
  mstreal redCut = S->getRedundancy();
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
      if (rmsds[k] <= cut) {
        rmsd = sols[k].getRMSD(); sols[k].setRMSD(rmsds[k]);      // overwrite with RMSD relative to current query
        if (redSet) matches.insert(sols[k], S->getRedundancy());  // insert copies the solution
        else matches.insert(sols[k]);
        if (maxSet && (matches.size() > maxN)) matches.erase(--matches.end());
        sols[k].setRMSD(rmsd);                                    // set RMSD back to the old value
      }
    }
    if (verb) {
      end = chrono::high_resolution_clock::now();
      searchTime = chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
      cout << "\tquick-search time " << searchTime << " ms" << std::endl;
    }
  }

  /* In order for us to forgo a raw search, we have to have found a suitable
   * neighboring cached query. If max number of matches had been set, we have to
   * have identified maxN matches to the current query within safe RMSD of it.
   * The matches being within the safe RMSD cutoff guarantees that no better
   * matches exist, and there being maxN of them means no better matches would
   * be needed anyway, since we have already git the maximum number allowed. If
   * maxN ahd not been set, then the test of suitability for a cached query
   * guarantees that all matches to the current query within the desired cutoff
   * are among the matches to the cached query. */
  bool doNewSearch = (bestComp == cache.end()) ||
                     (matches.size() == 0) ||               // test if there are any matches so that worsetRMSD() exists
                     (matches.worstRMSD() > safeRadius) ||  // if maxN is not set, this is guaranteed
                     (maxSet && (matches.size() < maxN));
  if (doNewSearch) {
    if (verb) {
      cout << "\t\tFAILED, need to search (" << matches.size() << " matches were found)...";
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
    S->pruneRedundancy(1.1);
    S->setRMSDCutoff(1.1*cut);

    /* The minimum number of matches limit is costly to include in the search,
     * but is sometimes not necesary (i.e., you find at least that many matches
     * under the given RMSD cutoff anyway). So first try without it. */
    for (int c = 0; c < (minSet ? 2 : 1); c++) {
      if (minSet) {
        if (c) S->setMinNumMatches(minN);
        else S->unsetMinNumMatches();
        if (c && verb) cout << "\t\t\tDID NOT WORK without min (found only " << matches.size() << " matches, while minN = " << minN << "), so repeating search..." << endl;
      }

      // -- perform the search and cache
      matches = S->search();
      if ((c == 0) && (matches.size() > 0)) {
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
    S->pruneRedundancy(redCut);
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

void selectAround(const vector<Residue*>& cenRes, int pm, Structure& frag, vector<int>& fragResIdx) {
  Structure* S = cenRes[0]->getChain()->getParent();
  vector<bool> selected(S->residueSize(), false);
  for (int i = 0; i < cenRes.size(); i++) {
    Residue& res = *(cenRes[i]);
    Chain* C = res.getChain();
    int ri = res.getResidueIndex();
    int Li = C->getResidue(C->residueSize() - 1).getResidueIndex(); // last residue index in the chain
    for (int k = ri - pm; k <= ri + pm; k++) {
      if ((k < 0) || (k > Li)) continue;
      selected[k] = true;
    }
  }
  Chain* newChain = frag.appendChain("A", true);
  for (int k = 0; k < selected.size(); k++) {
    if (selected[k]) {
      // where there is a break in the selection, start a new chain
      if ((newChain->residueSize() > 0) && ((!selected[k-1]) || (S->getResidue(k-1).getChain() != S->getResidue(k).getChain()))) {
        newChain = frag.appendChain("A", true);
      }
      newChain->appendResidue(new Residue(S->getResidue(k)));
      fragResIdx.push_back(k);
    }
  }
}

void addMatches(fasstCache& C, Structure& frag, const vector<int>& fragResIdx, vector<Structure*>& allMatches, vector<vector<Residue*> >& resTopo) {
//frag.writePDB("/tmp/pair.pdb"); cout << "fragment saved, RMSD = " << RMSDCalculator::rmsdCutoff(frag) << endl;
  FASST* F = C.getFASST();
  F->setQuery(frag, false);
  F->setRMSDCutoff(RMSDCalculator::rmsdCutoff(frag));
  F->setMaxNumMatches(10);
  F->setMinNumMatches(2);
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

  fasstSolutionSet matches = C.search(true);
  cout << "found " << matches.size() << " matches" << endl;
  for (auto it = matches.begin(); it != matches.end(); ++it) {
    allMatches.push_back(new Structure(F->getMatchStructure(*it, false, FASST::matchType::REGION)));
    Structure& match = *(allMatches.back());
    for (int k = 0; k < fragResIdx.size(); k++) {
      resTopo[fragResIdx[k]].push_back(&(match.getResidue(k)));
    }
  }
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
  op.addOption("rad", "compactness radius. Default will be based on protein length.");
  op.addOption("cyc", "number of iteration cycles (10 by default).");
  op.addOption("f", "a quoted, space-separated list of 0-initiated residue integers to fix.");
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
  F.pruneRedundancy(0.5);
  RotamerLibrary RL(op.getString("rLib"));
  int pmSelf = 2, pmPair = 1;
  int Ni = 1000;
  fusionScores scores;

  // TERMify loop
  Structure S = I;
  mstreal R0 = getRadius(I);
  mstreal Rf = op.getReal("rad", pow(I.residueSize() * 1.0, 1.0/3)*5.0);
  int Ncyc = op.getInt("cyc", 10);
  fstream out; MstUtils::openFile(out, op.getString("o") + ".traj.pdb", ios_base::out);
  out << "MODEL " << 0 << endl; S.writePDB(out); out << "ENDMDL" << endl;
  for (int c = 0; c < Ncyc; c++) {
    out << "MODEL " << c+1 << endl;
    cout << "Cycle " << c+1 << "..." << endl;
    /* --- Decorate the current conformation with TERMs --- */
    // first self TERMs
    cout << "Searching for self TERMs..." << endl;
    vector<vector<Residue*> > resTopo(I.residueSize());
    // by inserting the starting structure into the topology first, any fixed
    // regions will be fixed to these starting coordinates
    for (int ri = 0; ri < S.residueSize(); ri++) resTopo[ri].push_back(&(S.getResidue(ri)));
    vector<Structure*> allMatches;
    for (int ci = 0; ci < S.chainSize(); ci++) {
      Chain& C = S[ci];
      for (int ri = 0; ri < C.residueSize(); ri++) {
        Structure frag; vector<int> fragResIdx;
        selectAround({&C[ri]}, pmSelf, frag, fragResIdx);
// frag.writePDB("/tmp/self." + MstUtils::toString(ri) + ".pdb");
        cout << "TERM around " << C[ri] << endl;
        addMatches(cache, frag, fragResIdx, allMatches, resTopo);
      }
    }

    // then pair TERMs
    cout << "Searching for pair TERMs..." << endl;
    ConFind cfd(&RL, S);
    contactList L = cfd.getContacts(S, 0.01);
    vector<pair<Residue*, Residue*> > list = L.getOrderedContacts();
    for (int k = 0; k < list.size(); k++) {
      Residue* resA = list[k].first;
      Residue* resB = list[k].second;
      Structure frag; vector<int> fragResIdx;
      selectAround({resA, resB}, pmPair, frag, fragResIdx);
// frag.writePDB("/tmp/pair." + MstUtils::toString(k) + ".pdb");
      cout << "TERM around " << *resA << " x " << *resB << endl;
      addMatches(cache, frag, fragResIdx, allMatches, resTopo);
    }

// fstream ofs;
// MstUtils::openFile(ofs, "/tmp/allmatches.pdb", ios_base::out);
// S.writePDB(ofs);
// for (int i = 0; i < allMatches.size(); i++) allMatches[i]->writePDB(ofs);
// ofs.close();

    /* --- Fuse --- */
    fusionParams opts; opts.setNumIters(Ni); opts.setVerbose(false);
    opts.setGradientDescentFlag(true);
    opts.setRepFC(1);
    opts.setCompFC(0.1);
    mstreal compactnessRadius = Rf;
    // (R0 - Rf)*exp(-c/10.0) + Rf; // exponential scaling
    // (Rf*(c + 1) + R0*(Ncyc - c - 1))/Ncyc; // linear scaling

    opts.setCompRad(compactnessRadius);
    cout << "will be trying to combine structure to a radius of " << compactnessRadius << endl;
    Structure fused = Fuser::fuse(resTopo, scores, fixed, opts);
    cout << "\t" << scores << endl;
    // opts.setStartingStructure(fused);
    // opts.setGradientDescentFlag(false);
    // fused = Fuser::fuse(resTopo, scores, vector<int>(), opts);
    // cout << "\t" << scores << endl;

    // align based on the fixed part, if anything was fixed (for ease of visualization)
    if (fixed.size() > 0) {
      AtomPointerVector before = getBackbone(S, fixed);
      AtomPointerVector after = getBackbone(fused, fixed);
      rc.align(after, before, fused);
    }

    /* --- write intermediate result and clean up--- */
    fused.writePDB(out); out << "ENDMDL" << endl;
    S = fused;
    for (int mi = 0; mi < allMatches.size(); mi++) delete(allMatches[mi]);
  }
  out.close();
  S.writePDB(op.getString("o") + ".fin.pdb");
}
