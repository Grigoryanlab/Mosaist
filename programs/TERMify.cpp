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
    fasstSolutionSet search();

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
        fasstCachedResult() { rmsdCut = 0.0; solSet = NULL; priority = 0; }
        fasstCachedResult(const AtomPointerVector& q, fasstSolutionSet* sols, mstreal cut, vector<int> topo) {
          q.clone(query); rmsdCut = cut; solSet = sols; priority = 0; topology = topo;
        }
        fasstCachedResult(const fasstCachedResult& r) {
          query = r.query; rmsdCut = r.rmsdCut; solSet = r.solSet; priority = r.priority; topology = r.topology;
        }
        void upPriority() { priority++; }

        AtomPointerVector getQuery() const { return query; }
        fasstSolutionSet* getSolutions() const { return solSet; }
        mstreal getRMSDCutoff() const { return rmsdCut; }
        int getPriority() const { return priority; }
        vector<int> getTopology() const { return topology; }
        bool isSameTopology(const vector<int>& compTopo) const {
          if (topology.size() != compTopo.size()) return false;
          for (int i = 0; i < topology.size(); i++) {
            if (topology[i] != compTopo[i]) return false;
          }
          return true;
        }

        friend bool operator<(const fasstCachedResult& ri, const fasstCachedResult& rj) {
          if (ri.priority != rj.priority) return (ri.priority < rj.priority);
          if (ri.solSet->size() != rj.solSet->size()) return ri.solSet->size() < rj.solSet->size();
          if (ri.query.size() != rj.query.size()) return ri.query.size() < rj.query.size();
          return ri.query < rj.query;
        }

      private:
        AtomPointerVector query;
        vector<int> topology;
        mstreal rmsdCut;
        fasstSolutionSet* solSet;
        int priority;
    };
    int maxNumResults; // max number of searches to cache
    FASST* S;
    set<fasstCachedResult> cache;
    RMSDCalculator rc;
};

fasstCache::~fasstCache() {
  for (auto it = cache.begin(); it != cache.end(); ++it) {
    it->getQuery().deletePointers();
    delete(it->getSolutions());
  }
}

fasstSolutionSet fasstCache::search() {
  int n;
  if (S->isSufficientNumMatchesSet()) MstUtils::error("cannot cache when \"sufficient\" number of matches is set", "fasstCache::search");
  mstreal cut = S->getRMSDCutoff();
  fasstSolutionSet matches;

  /* See whether all matches for the current query, within the given cutoff, are
   * among the list of matches of some previously cached query. */
  vector<int> topo = fasstCache::getStructureTopology(S->getQuery());
  set<fasstCachedResult>::iterator bestComp = cache.end(); mstreal bestDist = -1;
  for (auto it = cache.begin(); it != cache.end(); ++it) {
    if (!it->isSameTopology(topo)) continue;
    mstreal r = rc.bestRMSD(S->getQuerySearchedAtoms(), it->getQuery());
    // if multiple cached options are available, pick the one with least extra
    if (it->getRMSDCutoff() - r >= cut) {
      if ((bestComp == cache.end()) || (bestDist > r)) {
        bestComp = it; bestDist = r;
      }
    }
  }
  if (bestComp != cache.end()) {
    // visit all matches of the most suitable cached result and find matches for current query
    const fasstSolutionSet& sols = *(bestComp->getSolutions());
    for (auto it = sols.begin(); it != sols.end(); ++it) {
      Structure match;
      S->getMatchStructure(*it, match, false, FASST::matchType::REGION, false);
      if (rc.bestRMSD(S->getQuerySearchedAtoms(), match.getAtoms()) <= cut) {
        matches.insert(*it);
      }
    }

    // up the priority of this cached result
    fasstCachedResult result(*bestComp);
    result.upPriority();
    cache.erase(bestComp); cache.insert(result);
  }

  /* If no suitable cached query was found OR if the best cached query did not
   * provide sufficient matches to match the minimum cutoff, will need to
   * initiate a new search, which will be cached. */
  bool doNewSearch = (bestComp == cache.end()) || (S->isMinNumMatchesSet() && (matches.size() < S->getMinNumMatches()));
  if (doNewSearch) {
    // -- loosen search criteria a bit to extract maximal value from search
    mstreal redCut; int maxN;
    bool maxSet = S->isMaxNumMatchesSet();
    bool redSet = S->isRedundancyCutSet();
    if (maxSet) {
      maxN = S->getMaxNumMatches();
      S->setMaxNumMatches(MstUtils::max(int(maxN*2 + 1), maxN + 500));
    }
    if (redSet) {
      redCut = S->getRedundancy();
      S->pruneRedundancy(1.1); // set to a value above 1 in order to accumulate sequence context within solutions
    }
    S->setRMSDCutoff(1.1*cut);

    // -- perform the search and cache
    matches = S->search();
    fasstCachedResult result(S->getQuerySearchedAtoms(),
                             new fasstSolutionSet(matches),
                             matches.rbegin()->getRMSD(),
                             fasstCache::getStructureTopology(S->getQuery()));
    cache.insert(result);
    if (cache.size() > maxNumResults) cache.erase(cache.begin()); // bump off the least used cached result if reached limit

    // -- reset search setting to their old values
    if (redSet) S->pruneRedundancy(redCut);
    S->setRMSDCutoff(cut);
    if (maxSet) S->setMaxNumMatches(maxN);
  }

  /* At this point, matches stores raw search results, whether from an actual
   * search or deduced from a previously cached result. But these may need to be
   * filtered in some way, if various limits were set in the FASST object. */
  // first, remove redundancy (in the same way it would be during the search)
  if (S->isRedundancyCutSet()) {
    fasstSolutionSet finalMatches;
    vector<fasstSolution*> orderedMatches = matches.orderByDiscovery();
    for (int i = 0; i < orderedMatches.size(); i++) {
      fasstSolution& sol = *(orderedMatches[i]);
      finalMatches.insert(sol,  S->getRedundancy());
    }
    matches = finalMatches;
  }

  // apply the RMSD filter if a looser raw search was run
  if (doNewSearch) {
    n = 0;
    for (auto it = matches.begin(); it != matches.end(); ) {
      if ((S->isMinNumMatchesSet() && (n < S->getMinNumMatches())) || (it->getRMSD() <= cut)) { n++; it++; }
      else { it = matches.erase(it); }
    }
  }

  // apply the max filter
  if (S->isMaxNumMatchesSet()) {
    if (matches.size() > S->getMaxNumMatches()) {
      auto it = matches.begin(); advance(it, S->getMaxNumMatches());
      while (it != matches.end()) it = matches.erase(it);
    }
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
  fasstSolutionSet matches = C.search();
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
        cout << "TERM around " << C[ri] << " ";
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
      cout << "TERM around " << *resA << " x " << *resB << " ";
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
