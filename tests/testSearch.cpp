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
    void removeOption(int k);
    mstreal bestCost() { return costs[bestCostRank]; }
    int bestChoice() { return rankToIdx[bestCostRank]; }
    int totNumOptions() { return costs.size(); }
    int numOptions() { return numIn; }
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
  rankToIdx.resize(_costs.size());
  for (int i = 0; i < rankToIdx.size(); i++) rankToIdx[i] = i;
  sort(rankToIdx.begin(), rankToIdx.end(), [&_costs](int i, int j) { return costs[i] < costs[j]; });
  costs.resize(_costs.size());
  for (int i = 0; i < rankToIdx.size(); i++) costs[i] = _costs[rankToIdx[i]];

  // sort the the rank-to-index mapping and keep track of indices to know the index-to-rank mapping
  idxToRank.resize(costs.size());
  for (int i = 0; i < idxToRank.size(); i++) idxToRank[i] = i;
  sort(idxToRank.begin(), idxToRank.end(), [&rankToIdx](int i, int j) { return costs[i] < costs[j]; });

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
      bestCostRank = -1;
      for (int i = bestCostRank; i < costs.size(); i++) {
        if (isIn[i]) { bestCostRank = i; break; }
      }
    }
    if ((bestCostRank < 0) || (costs[bestCostRank] > costs[idxToRank[k]])) {
      bestCostRank = idxToRank[k];
    }
  }
}


class MASTER {
  public:
    ~MASTER();
    MASTER() { recLevel = 0; setRMSDCutoff(1.0); setSearchType(2); ps = NULL; }
    void setQuery(const string& pdbFile);
    void setTarget(const string& pdbFile);
    void prepForSearch();
    void setRMSDCutoff(mstreal cut);
    void setSearchType(int _searchType);
    int numOptions() { return remOptions[curSeg].size(); }
    bool outOfOptions() { return remOptions[curSeg].empty(); }
    bool visitNextOption();

  protected:
    static bool parseChain(const Chain& S, AtomPointerVector& searchable);
    mstreal currentAlignmentResidual(bool recompute); // TODO
    template <class T>
    void inPlaceSetDifference(set<int, T> in, const vector<int>& notin);

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

    // remOptions[L][i] is a set of alignments for segment L+i at recursion level
    // L, which are stored sorted by their own residual, through the optList
    // data structure. Note that segments 0 through L-1 have already been place
    // at recursion level L.
    vector<vector<optList> > > remOptions;

    // ccTol[i][j], j > i, is the acceptable tolerance (the delta) on the
    // center-to-center distance between segments i and j at recursion level j
    vector<vector<mstreal> > ccTol;

    // alignment indices for segments visited up to the current recursion level
    vector<int> currAlignment;

    // the residual of the above alignment, computed and stored
    mstreal currResidual;

    // set of solutions, sorted by RMSD
    set<pair<vector<int, mstreal> > solutions; // TODO: need to store target info, order added, and implemenet comparison

    // ProximitySearch for finding nearby centroids of target segments
    vector<ProximitySearch> ps; // TODO: make a default and copy constructors for ProximitySearch

    // various solution constraints
    mstreal rmsdCut, residualCut;
};

int main(int argc, char *argv[]) {
  MASTER S;
  S.setQuery("query.pdb");
  S.setTarget("target.pdb");
  S.prepForSearch();
  S.setRMSDCutoff(1.0);

  while (!S.outOfOptions()) {
    S.visitNextOption();
  }
}


MASTER::~MASTER() {
  if (ps != NULL) delete(ps);
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
    if (!MASTER::parseChain(queryStruct[i], query[i])) {
      MstUtils::error("could not set query, because some atoms for the specified search type were missing", "MASTER::setQuery");
    }
  }
  currAlignment.resize(query.size(), -1);
  segmentResidualsInAlignment.resize(query.size(), -1);
}

void MASTER::setTarget(const string& pdbFile) {
  targetStruct.reset();
  targetStruct.readPDB(pdbFile);
  target.resize(0);
  // we don't care about the chain topology of the target, so append all residues
  for (int i = 0; i < targetStruct.chainSize(); i++) {
    MASTER::parseChain(targetStruct[i], target);
  }
  setRMSDCutoff(rmsdCut); // update the max residual
}

bool MASTER::parseChain(const Chain& C, AtomPointerVector& searchable) {
  bool foundAll = true;
  for (int i = 0; i < S.residueSize(); i++) {
    Residue& res = S.getResidue(i);
    switch(searchType) {
      case 1:
      case 2:
        AtomPointerVector bb;
        for (int k = 0; k < searchableAtoms.size(); k++) {
          Atom* a = res.findAtom(searchableAtoms[k], false);
          if (a == NULL) { foundAll = false; break; }
          bb.push_back(a);
        }
        if (bb.size() == bba.size()) searchable.insert(searchable.end(), bb.begin(), bb.end());
        break
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
      atomsPerRes = 4;
      break
    default:
      MstUtils::error("uknown search type '" + MstUtils::toString(searchType) + "' specified", "MASTER::setSearchType");
  }
  atomsPerRes = searchableAtoms.size();
}

void MASTER::prepForSearch() {
  if ((query.size() == 0) || (target.size() == 0)) {
    MstUtils::error("query and target must be set before starting search", "MASTER::prepForSearch");
  }
  solutions.resize(0);

  // align every segment onto every admissible location on the target
  RMSDCalculator rc;
  segmentResiduals.resize(query.size());
  ps.resize(query.size());
  for (int i = 0; i < query.size(); i++) {
    ps[i] = ProximitySearch(target, 3.0, false);
    vector<Residue>& seg = query[i];
    int Na = min(target.size(), target.size() - seg.size() + 1)/atomsPerRes; // number of possible alignments
    segmentResiduals[i].resize(Na);
    for (int j = 0; j < Na; j++) {
      // NOTE: can save on this in several ways:
      // 1. the centroid calculation is effectively already done inside RMSDCalculator::bestRMSD
      // 2. updating just one atom involves a simple centroid adjustment, rather than recalculation
      // 3. is there a speedup to be gained from re-calculating RMSD with one atom updated only?
      AtomPointerVector targSeg = target.subvector(j*atomsPerRes, j*atomsPerRes + query[i].size() - 1);
      segmentResiduals[i][j] = pow(rc.bestRMSD(query[i], targSeg), 2) * targSeg.size();
      ps[i].addPoint(targSeg.getGeometricCenter(), j);
    }
  }

  // initialize remOptions; all options are available at top level
  remOptions.resize(query.size());
  for (int L = 0; L < query.size(); L++) {
    remOptions[L].resize(query.size() - L);
    for (int i = L; i < query.size(); i++) {
      remOptions[L][i].setOptions(segmentResiduals[i], L == 0);
    }
  }

  // current residual starts with 0, since nothing is aligned yet
  currResidual = 0;

  // initialize center-to-center tolerances
  ccTol.resize(query.size());
  for (int i = 0; i < query.size(); i++) {
    ccTol[i].resize(query.size(), -1.0); // unnecessary (unused) entries will stay as -1
    for (int j = i+1; j < query.size(); j++) {
      ccTol[i][j] = centToCentTol(i, j);
    }
  }
}

mstreal MASTER::lowBoundOnRemainder() {
  mstreal bound;
  for (int i = level; i < query.size(); i++) bound += remOptions[recLevel][i].bestCost();
  return bound;
}

mstreal MASTER::centToCentTol(int i, int j) {
  // the second parameter to currentAlignmentResidual means do not recalculate
  // TODO: validate the bound
  return sqrt((((residualCut - currentAlignmentResidual(false) - lowBoundOnRemainder(false ??))) * (query[i].size() + query[j].size())) / (query[i].size() * query[j].size()));
}

bool MASTER::visitNextOption() {
  // Have to do three things:
  // 1. pick the best choice (from available ones) for the current segment, and
  // remove it from the list of options
  if (remOptions[recLevel][recLevel].empty()) {
    upRecLevel();
    return false;
  }
  currAlignment[recLevel] = remOptions[recLevel][recLevel].bestChoice();
  remOptions[recLevel][recLevel].removeOption(currAlignment[recLevel]);

  // 2. compute the total residual from the current alignment
  mstreal curBound = currentAlignmentResidual(true) + lowBoundOnRemainder(true);
  if (curBound > residualCut) {
    return false;
  }

  // 3. update remOptions[recLevel+1][recLevel+1...end] based on the
  // choice in 1. Note that the set of options on the next recursion level is a
  // subset of the set of options on the previous level.
  // TODO: MUST DECIDE WHEN TO DO A PROXIMITY SEARCH VERSUS WHEN TO DO A SUBSET
  recLevel++;
  int remSegs = query.size() - recLevel;
  if (remSegs > 0) {
    mstreal dij, de, dU, dL, d;
    while (true) {
      bool updated = false;
      for (int i = recLevel; i < query.size(); i++) {
        // IDEA: write a function in RMSDCalculator that does optimization only over
        // rotation, not centering (just skips the centering part). Then, pre-center
        // all the query sub-structures at all recursion levels and store them. Will
        // not have to do this later during optimization. Also, create ROOM for storing
        // target atoms at all recursion levels, but these will have to be assigned
        // each time you do total RMSD. Actually, just the last segment added will need
        // to be assigned. Finally, because you can use the function that does not
        // center, you can make sure that when you are scanning individual segments,
        // you are centering them once to compute the centroids and not again in the
        // RMSDCalculator.

        optList& remSet = remOptions[recLevel][i - recLevel];
        for (int j = 0; j < recLevel; j++) {                    // j runs over already placed segments
          CartesianPoint cj = ps[j].getPoint(currAlignment[j]); // the current centroid of segment j
          dij = centToCentDist[i][j];                        // TODO: need to populate this upon loading the query
          de = centToCentTol(i, j);
          mstreal& dePrev = ccTol[j][i];
          dU = dij - de; dL = dij + de;
          int numLocs = remSet.size();

          // The first time we do a proximity search to find some options for
          // segment i. But if we already have some options for the segment, we
          // need to do a set difference to remove no-longer-valid ones
          if (numLocs == 0) {
            vector<int> okLocations;
            ps[i].pointsWithin(cj, dij - de, dij + de, okLocations);
            for (int k = 0; k < okLocations.size(); k++) {
              remSet.addOption(okLocations[k]);
            }
          } else {
            vector<int> badLocations; mstreal eps = 10E-8;
            ps[i].pointsWithin(cj, dij - dePrev, dij - de - eps, badLocations);
            ps[i].pointsWithin(cj, dij + de + eps, dij + dePrev, badLocations);
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
    recLevel++;
  } else {
    // if at the lowest recursion level already, then record the solution
    solutions.insert(pair<vector<int, mstreal>(currAlignment, currResidual));
  }

  return true;
}

void MASTER::upRecLevel() {
  if (recLevel > 0) {
    currAlignment[recLevel] = -1;
    segmentResidualsInAlignment[recLevel] = -1;
    recLevel--;
  }
}

mstreal MASTER::currentAlignmentResidual(bool recompute) {
  if (!recompute) return currResidual;
  if (recLevel == 0) currResidual = segmentResiduals[recLevel][currAlignment[recLevel]];
  // JOIN THE APPROPRIATE segments
  // COMPUTE RMSD
  // TODO: ??????????????
  return currResidual;
}

template <class T>
void MASTER::inPlaceSetDifference(set<int, T> in, const vector<int>& notin) {
  std::set<int, T>::iterator it = in.begin();
  int i = 0;
  while ((it != in.end()) && (i != notin.size())) {
      if (*it < *it2) {
          set_1.erase(it1++);
      } else if (*it2 < *it1) {
          ++it2;
      } else { // *it1 == *it2
              ++it1;
              ++it2;
      }
  }
  // Anything left in set_1 from here on did not appear in set_2,
  // so we remove it.
  set_1.erase(it1, set_1.end());
}
