#include "msttypes.h"

using namespace MST;

struct compAlignments {
  // Compare by RMSD first, and if the same, then by alignment index. NOTE: for
  // a given segment, no two different alignments can have the same index.
  bool operator() (const pair<int, mstreal>& lhs, const pair<int, mstreal>& rhs) const {
    if (lhs.second == rhs.second) return lhs.first < rhs.first;
    return lhs.second < rhs.second;
  }
};

struct compSolutions {
  bool operator() (const pair<vector<int>, mstreal>& lhs, const pair<vector<int>, mstreal>& rhs) const {
    if (lhs.second == rhs.second) return lhs.first < rhs.first;
    int dim = lhs.first.size();
    for (int i = 0; i < dim - 1; i++) {
      if (lhs.first[i] != rhs.first[i]) return lhs.first[i] < rhs.first[i];
    }
    return lhs.first.back() < rhs.first.back();
  }
};

/* A class for storing a sub-set of a fixed set of options, each with a fixed
 * cost. The options are stored sorted by cost, so that the best option is
 * alwasys quickly availabe. Insertion is constant time. Deletion is constant
 * on average. */
class optList {
  public:
    /* Takes an unsorted list of costs and initializes all necessary interal
     * structures. Elements of the set of options will be known externally by
     * their index into this array (the world index), although internally costs
     * will also be stored (plus mappings between the two orders), so that the
     * best available option will always be known. If second argument given as
     * true, will start with all options being available. Otherwise, none will
     * be marked as available. */
    optList(const vector<mstreal>& _costs, bool add = false) { setOptions(_costs, add); }
    optList() { bestCostRank = -1; numIn = 0; }

    void setOptions(const vector<mstreal>& _costs, bool add = false);
    void addOption(int k);
    void removeAllOptions();
    void removeOption(int k);
    mstreal bestCost() { return costs[bestCostRank]; }
    int bestChoice() { return rankToIdx[bestCostRank]; }
    int totNumOptions() { return costs.size(); }
    int numOptions() { return numIn; }
    int size() { return numIn; }
    bool empty() { return (numIn == 0); }

  private:
    // costs, sorted in ascending order
    vector<mstreal> costs;

    // is each option (by world index) currently available?
    vector<bool> isIn;

    // back-and-forth mappings between rank by cost (ascending order) and the
    // flat index known to the word
    vector<int> idxToRank, rankToIdx;

    // of the available options, the lowest-cost rank
    int bestCostRank;

    // number of currently available options
    int numIn;
};

void optList::setOptions(const vector<mstreal>& _costs, bool add) {
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

void optList::addOption(int k) {
  // if ((k < 0) || (k >= costs.size())) MstUtils::error("out-of-range index specified: " + MstUtils::toString(k), "optList::insertOption(int)");
  if (!isIn[k]) {
    numIn++;
    isIn[k] = true;
    // update best, if the newly inserted option is better
    if ((bestCostRank < 0) || (costs[bestCostRank] > costs[idxToRank[k]])) {
      bestCostRank = idxToRank[k];
    }
  }
}

void optList::removeOption(int k) {
  // if ((k < 0) || (k >= costs.size())) MstUtils::error("out-of-range index specified: " + MstUtils::toString(k), "optList::removeOption(int)");
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

void optList::removeAllOptions() {
  for (int i = 0; i < isIn.size(); i++) isIn[i] = false;
  numIn = 0;
  bestCostRank = -1;
}

class MASTER {
  public:
    ~MASTER();
    MASTER() { recLevel = 0; setRMSDCutoff(1.0); setSearchType(2); }
    void setQuery(const string& pdbFile);
    void setTarget(const string& pdbFile);
    void prepForSearch();
    void setRMSDCutoff(mstreal cut);
    void setSearchType(int _searchType);
    bool visitNextOption();
    int numSolutions() { return solutions.size(); }

  protected:
    bool parseChain(const Chain& S, AtomPointerVector& searchable);
    mstreal currentAlignmentResidual(bool compute);   // computes the accumulated residual up to and including segment recLevel
    mstreal boundOnRemainder(bool compute);           // computes the lower bound expected from segments recLevel+1 and on
    int resToAtomIdx(int resIdx) { return resIdx * atomsPerRes; }
    int atomToResIdx(int atomIdx) { return atomIdx / atomsPerRes; }
    void fillTargetMask();
    mstreal centToCentTol(int i, int j);

  private:
    Structure queryStruct, targetStruct;
    vector<AtomPointerVector> query;   // just the part of the query that will be sought, split by segment
    AtomPointerVector target;          // just the part of the target structure that will be searched over

    // segmentResiduals[i][j] is the residual of the alignment of segment i, in which
    // its starting residue aligns with the residue index j in the target
    vector<vector<mstreal> > segmentResiduals;

    int recLevel;
    int searchType, atomsPerRes;
    vector<string> searchableAtoms;

    // remOptions[L][i] is a set of alignments for segment i (i >= L) at recursion level
    // L, which are stored sorted by their own residual, through the optList
    // data structure. Note that segments 0 through L-1 have already been place
    // at recursion level L.
    vector<vector<optList> > remOptions;

    // ccTol[i][j], j > i, is the acceptable tolerance (the delta) on the
    // center-to-center distance between segments i and j at recursion level j
    vector<vector<mstreal> > ccTol;

    // alignment indices for segments visited up to the current recursion level
    vector<int> currAlignment;

    // the residual of the above alignment, computed and stored
    mstreal currResidual;

    // the bound on the parts remaining to align, computed and stored
    mstreal currRemBound;

    // center-to-center distances between segments of the query
    vector<vector<mstreal> > centToCentDist;

    // set of solutions, sorted by RMSD
    set<pair<vector<int>, mstreal> > solutions; // TODO: need to store target info, order added, and implemenet comparison

    // ProximitySearch for finding nearby centroids of target segments
    vector<ProximitySearch> ps;

    // various solution constraints
    mstreal rmsdCut, residualCut;

    // Atom subsets needed at different recursion levels. So queryMasks[i] stores
    // all atoms of the first i+1 segments of the query combined. The same for
    // targetMasks, although (of course), the content of the latter will change
    // depending on the alignment
    vector<AtomPointerVector> queryMasks, targetMasks;

    RMSDCalculator RC;
};

int main(int argc, char *argv[]) {
  MASTER S;
  S.setQuery("/tmp/query.pdb");
  S.setTarget("/tmp/target.pdb");
  S.prepForSearch();
  S.setRMSDCutoff(2.0);

  while (S.visitNextOption()) {}
  cout << "found " << S.numSolutions() << " solutions" << endl;
}


MASTER::~MASTER() {
  // need to delete atoms only on the lowest level of recursion, because at
  // higher levels we point to the same atoms
  targetMasks.back().deletePointers();
}

void MASTER::setRMSDCutoff(mstreal cut) {
  rmsdCut = cut;
  residualCut = cut*cut*query.size();
}

/* This function sets up the query, in the process deciding which part of the
 * query is really searchable (e.g., backbone). Various search types an be added
 * in the future to provide search capabilities over different parts of the
 * structure. The current implementation stipulates that the searchable part
 * involve the same number of atoms for each residue, so that atoms can be
 * stored in flat arrays for efficiency. But these could presumably be
 * interpreted differently. E.g., only some residue matches (parents of atoms)
 * may be accepted or some of the atoms can be dummy/empty ones. */
void MASTER::setQuery(const string& pdbFile) {
  queryStruct.reset();
  queryStruct.readPDB(pdbFile);
  // auto-splitting segments by connectivity, but can do diffeerntly
  queryStruct = queryStruct.reassignChainsByConnectivity();
  query.resize(queryStruct.chainSize());
  for (int i = 0; i < queryStruct.chainSize(); i++) {
    query[i].resize(0);
    if (!parseChain(queryStruct[i], query[i])) {
      MstUtils::error("could not set query, because some atoms for the specified search type were missing", "MASTER::setQuery");
    }
    MstUtils::assert(query[i].size() > 0, "query contains empty segment(s)", "MASTER::setQuery");
  }
  currAlignment.resize(query.size(), -1);
}

void MASTER::setTarget(const string& pdbFile) {
  targetStruct.reset();
  targetStruct.readPDB(pdbFile);
  target.resize(0);
  // we don't care about the chain topology of the target, so append all residues
  for (int i = 0; i < targetStruct.chainSize(); i++) {
    parseChain(targetStruct[i], target);
  }
  MstUtils::assert(target.size() > 0, "empty target passed", "MASTER::setTarget");
  setRMSDCutoff(rmsdCut); // update the max residual
}

bool MASTER::parseChain(const Chain& C, AtomPointerVector& searchable) {
  bool foundAll = true;
  for (int i = 0; i < C.residueSize(); i++) {
    Residue& res = C.getResidue(i);
    switch(searchType) {
      case 1:
      case 2: {
        AtomPointerVector bb;
        for (int k = 0; k < searchableAtoms.size(); k++) {
          Atom* a = res.findAtom(searchableAtoms[k], false);
          if (a == NULL) { foundAll = false; break; }
          bb.push_back(a);
        }
        if (bb.size() == searchableAtoms.size()) searchable.insert(searchable.end(), bb.begin(), bb.end());
        break;
      }
      default:
        MstUtils::error("uknown search type '" + MstUtils::toString(searchType) + "' specified", "MASTER::parseStructure");
    }
  }
  return foundAll;
}

void MASTER::setSearchType(int _searchType) {
  searchType = _searchType;
  switch(searchType) {
    case 1:
      searchableAtoms =  {"CA"};
      break;
    case 2:
      searchableAtoms =  {"N", "CA", "C", "O"};
      break;
    default:
      MstUtils::error("uknown search type '" + MstUtils::toString(searchType) + "' specified", "MASTER::setSearchType");
  }
  atomsPerRes = searchableAtoms.size();
}

void MASTER::prepForSearch() {
  if ((query.size() == 0) || (target.size() == 0)) {
    MstUtils::error("query and target must be set before starting search", "MASTER::prepForSearch");
  }
  solutions.clear();

  // align every segment onto every admissible location on the target
  segmentResiduals.resize(query.size());
  ps.resize(query.size());
  for (int i = 0; i < query.size(); i++) {
    ps[i] = ProximitySearch(target, 3.0, false);
    AtomPointerVector& seg = query[i];
    int Na = atomToResIdx(target.size()) - atomToResIdx(seg.size()) + 1; // number of possible alignments
    segmentResiduals[i].resize(Na);
    for (int j = 0; j < Na; j++) {
      // NOTE: can save on this in several ways:
      // 1. the centroid calculation is effectively already done inside RMSDCalculator::bestRMSD
      // 2. updating just one atom involves a simple centroid adjustment, rather than recalculation
      // 3. is there a speedup to be gained from re-calculating RMSD with one atom updated only?
      AtomPointerVector targSeg = target.subvector(resToAtomIdx(j), resToAtomIdx(j) + query[i].size());
      segmentResiduals[i][j] = pow(RC.bestRMSD(query[i], targSeg), 2) * targSeg.size();
      ps[i].addPoint(targSeg.getGeometricCenter(), j);
    }
  }

  // initialize remOptions; all options are available at top level
  remOptions.resize(query.size());
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
  currRemBound = boundOnRemainder(true);

  // initialize center-to-center tolerances
  ccTol.resize(query.size());
  for (int i = 0; i < query.size(); i++) {
    ccTol[i].resize(query.size(), -1.0); // unnecessary (unused) entries will stay as -1
    for (int j = i+1; j < query.size(); j++) {
      ccTol[i][j] = centToCentTol(i, j);
    }
  }

  // make room for various atom masks
  queryMasks.resize(query.size());
  targetMasks.resize(query.size());
  for (int i = 0; i < query.size(); i++) {
    for (int L = i; L < query.size(); L++) {
      for (int j = 0; j < query[i].size(); j++) {
        queryMasks[L].push_back(query[i][j]);
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

  // center-to-center distances
  centToCentDist.resize(query.size(), vector<mstreal>(query.size(), 0));
  for (int i = 0; i < query.size(); i++) {
    for (int j = i+1; j < query.size(); j++) {
      centToCentDist[i][j] = query[i].getGeometricCenter().distance(query[j].getGeometricCenter());
      centToCentDist[j][i] = centToCentDist[i][j];
    }
  }
}

mstreal MASTER::boundOnRemainder(bool compute) {
  if (compute) {
    currRemBound = 0;
    for (int i = recLevel + 1; i < query.size(); i++) currRemBound += remOptions[recLevel][i].bestCost();
  }
  return currRemBound;
}

mstreal MASTER::centToCentTol(int i, int j) {
  // TODO: validate the bound
  return sqrt((((residualCut - currentAlignmentResidual(false) - boundOnRemainder(false))) * (query[i].size() + query[j].size())) / (query[i].size() * query[j].size()));
}

void MASTER::fillTargetMask() {
  int N = targetMasks[recLevel].size();
  int n = query[recLevel].size();
  int si = resToAtomIdx(currAlignment[recLevel]);
  for (int i = 0; i < n; i++) {
    targetMasks[recLevel][N - n + i]->setCoor(target[si + i]->getCoor());
  }
}

bool MASTER::visitNextOption() {
  // Have to do three things:
  // 1. pick the best choice (from available ones) for the current segment,
  // remove it from the list of options, and move onto the next recursion level
  if (remOptions[recLevel][recLevel].empty()) {
    currAlignment[recLevel] = -1;
    if (recLevel > 0) {
      for (int i = recLevel; i < query.size(); i++) remOptions[recLevel][recLevel].removeAllOptions();
      recLevel--;
      return false; // search exhausted
    }
    return true;
  }
  currAlignment[recLevel] = remOptions[recLevel][recLevel].bestChoice();
  fillTargetMask();
  remOptions[recLevel][recLevel].removeOption(currAlignment[recLevel]);

  // 2. compute the total residual from the current alignment
  mstreal curBound = currentAlignmentResidual(true) + boundOnRemainder(true);
  if (curBound > residualCut) {
    return true;
  }

  // 3. update update remaining options for subsequent segments based on the
  // newly made choice. The set of options on the next recursion level is a
  // subset of the set of options on the previous level.
  int remSegs = query.size() - (recLevel + 1);
  if (remSegs > 0) {
    recLevel++;
    mstreal dij, de, dU, dL, d;
    while (true) {
      bool updated = false;
      for (int i = recLevel; i < query.size(); i++) {
        // IDEA: write a function in RMSDCalculator that does optimization only over
        // rotation, not centering (just skips the centering part). Then, pre-center
        // all the query sub-structures at all recursion levels and store them. Will
        // not have to do this later during optimization. Finally, because you can use the function that does not
        // center, you can make sure that when you are scanning individual segments,
        // you are centering them once to compute the centroids and not again in the
        // RMSDCalculator.

        optList& remSet = remOptions[recLevel][recLevel];
        for (int j = 0; j < recLevel; j++) {                    // j runs over already placed segments
          CartesianPoint cj = ps[j].getPoint(currAlignment[j]); // the current centroid of segment j
          dij = centToCentDist[i][j];
          de = centToCentTol(i, j);
          mstreal& dePrev = ccTol[j][i];
          dU = dij - de; dL = dij + de;
          int numLocs = remSet.size();

          // The first time we do a proximity search to find some options for
          // segment i. But if we already have some options for the segment, we
          // need to do a set difference to remove no-longer-valid ones
          if (numLocs == 0) {
            vector<int> okLocations;
            ps[i].pointsWithin(cj, dij - de, dij + de, &okLocations);
            for (int k = 0; k < okLocations.size(); k++) {
              remSet.addOption(okLocations[k]);
            }
          } else {
            vector<int> badLocations; mstreal eps = 10E-8;
            ps[i].pointsWithin(cj, dij - dePrev, dij - de - eps, &badLocations);
            ps[i].pointsWithin(cj, dij + de + eps, dij + dePrev, &badLocations);
            for (int k = 0; k < badLocations.size(); k++) {
              remSet.removeOption(badLocations[k]);
            }
          }
          dePrev = dij;
          if (numLocs != remSet.size()) updated = true;
        }
      }
      if (!updated) break;
    }
  } else {
    // if at the lowest recursion level already, then record the solution
    solutions.insert(pair<vector<int>, mstreal>(currAlignment, currResidual));
  }

  return true;
}

mstreal MASTER::currentAlignmentResidual(bool compute) {
  if (compute) {
    if (recLevel == 0) {
      currResidual = segmentResiduals[recLevel][currAlignment[recLevel]];
    } else {
      mstreal rmsd = RC.bestRMSD(queryMasks[recLevel], targetMasks[recLevel]);
      currResidual = rmsd*rmsd;
    }
  }
  return currResidual;
}
