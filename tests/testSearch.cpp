#include "msttypes.h"
#include "mstoptions.h"
#include "msttransforms.h"
#include <list>
#include <chrono>

using namespace MST;

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
    // after this, the only remaining allowed options will be those that were
    // allowed before AND are in the specified list
    void intersectOptions(const vector<int>& opts);
    // keep only options with indices less than or equal to the given one
    void constrainLE(int idx);
    // keep only options with indices greater than or equal to the given one
    void constrainGE(int idx);
    // wipe all options
    void removeAllOptions();
    // copy the valid options from a different set of the same options
    void copyIn(const optList& opt);
    void removeOption(int k);
    mstreal bestCost() { return costs[bestCostRank]; }
    int bestChoice() { return rankToIdx[bestCostRank]; }
    int totNumOptions() { return costs.size(); }
    int numOptions() { return numIn; }
    int size() { return numIn; }
    bool empty() { return (numIn == 0); }
    bool consistencyCheck();

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

bool optList::consistencyCheck() {
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

void optList::copyIn(const optList& opt) {
  isIn = opt.isIn;
  bestCostRank = opt.bestCostRank;
  numIn = opt.numIn;
}

void optList::intersectOptions(const vector<int>& opts) {
  vector<bool> isGiven(isIn.size(), false);
  for (int i = 0; i < opts.size(); i++) isGiven[opts[i]] = true;
  for (int i = 0; i < isIn.size(); i++) {
    if (isIn[i] && !isGiven[i]) removeOption(i);
  }
}

void optList::constrainLE(int idx) {
  for (int i = idx+1; i < isIn.size(); i++) removeOption(i);
}

void optList::constrainGE(int idx) {
  for (int i = 0; i < MstUtils::min(idx, (int) isIn.size()); i++) removeOption(i);
}

class fasstSolution {
  public:
    fasstSolution(const vector<int>& _alignment, mstreal _rmsd, int _target, int _foundOrder) {
      alignment = _alignment; rmsd = _rmsd; targetIndex = _target; foundOrder = _foundOrder;
    }
    fasstSolution(const fasstSolution& _sol) {
      alignment = _sol.alignment; rmsd = _sol.rmsd; targetIndex = _sol.targetIndex; foundOrder = _sol.foundOrder;
    }

    mstreal getRMSD() const { return rmsd; }
    int getTargetIndex() const { return targetIndex; }
    vector<int> getAlignment() const { return alignment; }

    friend bool operator<(const fasstSolution& si, const fasstSolution& sj) {
      if (si.rmsd != sj.rmsd) return (si.rmsd < sj.rmsd);
      return si.foundOrder < sj.foundOrder;
    }

    friend ostream& operator<<(ostream &_os, const fasstSolution& _sol) {
      _os << std::setprecision(6) << std::fixed << _sol.rmsd << " " << _sol.targetIndex << " [" << MstUtils::vecToString(_sol.alignment, ", ") << "]";
      return _os;
    }

  private:
    vector<int> alignment;
    mstreal rmsd;
    int targetIndex, foundOrder;
};

/* FASST -- Fast Algorithm for Searching STructure */
class FASST {
  public:
    enum matchType { REGION = 1, FULL, WITHGAPS };
    enum searchType { CA = 1, FULLBB };
    enum targetFileType { PDB = 1, BINDATABASE };
    class targetInfo {
      public:
        targetInfo(const string& _file, targetFileType _type, int _index, bool _memSave) {
          file = _file; type = _type; index = _index; memSave = _memSave;
        }
        targetInfo(const targetInfo& I) {
          file = I.file; type = I.type; index = I.index; memSave = I.memSave;
        }
        string file;         // source file
        targetFileType type; // file type (PDB or database)
        int index;           // location info within file
        bool memSave;        // was the target read with memory save on?
    };
    ~FASST();
    FASST();
    void setQuery(const string& pdbFile);
    Structure getQuery() { return queryStruct; }
    void addTarget(const string& pdbFile);
    void addTargets(const vector<string>& pdbFiles);
    void setRMSDCutoff(mstreal cut) { rmsdCutRequested = cut; }
    void setSearchType(searchType _searchType);
    void setMemorySaveMode(bool _memSave) { memSave = _memSave; }
    void setGridSpacing(mstreal _spacing) { gridSpacing = _spacing; updateGrids = true; }
    void search();
    int numMatches() { return solutions.size(); }
    void setMaxNumMatches(int _max);
    void setMinNumMatches(int _min);
    int getMaxNumMatches() { return maxNumMatches; }
    int getMinNumMatches() { return minNumMatches; }
    bool isMaxNumMatchesSet() { return (maxNumMatches > 0); }
    bool isMinNumMatchesSet() { return (minNumMatches >= 0); }
    set<fasstSolution> getMatches() { return solutions; }

    bool gapConstrained(int i, int j) { return (minGapSet[i][j] || maxGapSet[i][j]); }
    void setMinGap(int i, int j, int gapLim); // target topology: [segment i] [gap of at list gapLim long] [segment j]
    void setMaxGap(int i, int j, int gapLim); // target topology: [segment i] [gap of at most gapLim long] [segment j]
    void resetGapConstraints();
    bool gapConstraintsExist() { return gapConstSet; }
    bool validateSearchRequest(); // make sure all user specified requirements are consistent
    void getMatchStructure(const fasstSolution& sol, Structure& match, bool detailed = false, matchType type = matchType::REGION);
    Structure getMatchStructure(const fasstSolution& sol, bool detailed = false, matchType type = matchType::REGION);
    void getMatchStructures(const vector<fasstSolution>& sols, vector<Structure>& matches, bool detailed = false, matchType type = matchType::REGION);
    // TODO
    // void writeDatabase(const string& dbFile);
    // void readDatabase(const string& dbFile);

  protected:
    void setCurrentRMSDCutoff(mstreal cut);
    void prepForSearch(int ti);
    bool parseChain(const Chain& S, AtomPointerVector& searchable);
    mstreal currentAlignmentResidual(bool compute);   // computes the accumulated residual up to and including segment recLevel
    mstreal boundOnRemainder(bool compute);           // computes the lower bound expected from segments recLevel+1 and on
    int resToAtomIdx(int resIdx) { return resIdx * atomsPerRes; }
    int atomToResIdx(int atomIdx) { return atomIdx / atomsPerRes; }
    mstreal centToCentTol(int i, int j, bool recomputeResidual = false, bool recomputeBound = false);
    void rebuildProximityGrids();
    void stripSidechains(Structure& S);

  private:
    Structure queryStruct;
    vector<Structure*> targetStructs;
    vector<AtomPointerVector> query;         // just the part of the query that will be sought, split by segment
    vector<AtomPointerVector> targets;       // just the part of the target structure that will be searched over
    vector<targetInfo> targetSource;         // where each target was read from (in case need to re-read it)
    mstreal xlo, ylo, zlo, xhi, yhi, zhi;    // bounding box of the search database
    vector<Transform> tr;                    // transformations from the original frame to the common frames of reference for each target
    bool memSave;                            // save memory by storing only the backbone of targets?

    vector<vector<int> > minGap, maxGap;     // minimum and maximum sequence separations allowed between each pair of segments
    vector<vector<bool> > minGapSet, maxGapSet;
    bool gapConstSet;

    // segmentResiduals[i][j] is the residual of the alignment of segment i, in which
    // its starting residue aligns with the residue index j in the target
    vector<vector<mstreal> > segmentResiduals;

    int recLevel;
    int atomsPerRes;
    searchType type;
    vector<string> searchableAtomTypes;
    int querySize;
    int maxNumMatches, minNumMatches;

    // remOptions[L][i] is a set of alignments for segment i (i >= L) at recursion level
    // L, which are stored sorted by their own residual, through the optList
    // data structure. Note that segments 0 through L-1 have already been place
    // at recursion level L.
    vector<vector<optList> > remOptions;

    // ccTol[L][i][j], j > i, is the acceptable tolerance (the delta) on the
    // center-to-center distance between segments i and j at recursion level L
    vector<vector<vector<mstreal> > > ccTol;

    // alignment indices for segments visited up to the current recursion level
    vector<int> currAlignment;

    // the residual of the above alignment, computed and stored
    mstreal currResidual;

    // the bound on the parts remaining to align, computed and stored
    mstreal currRemBound;

    // center-to-center distances between segments of the query
    vector<vector<mstreal> > centToCentDist;

    // set of solutions, sorted by RMSD
    set<fasstSolution> solutions;

    // ProximitySearch for finding nearby centroids of target segments (there
    // will be one ProximitySearch object per query segment)
    vector<ProximitySearch> ps;

    // various solution constraints
    mstreal rmsdCutRequested, rmsdCut, residualCut;

    // Atom subsets needed at different recursion levels. So queryMasks[i] stores
    // all atoms of the first i+1 segments of the query combined. The same for
    // targetMasks, although (of course), the content of the latter will change
    // depending on the alignment
    vector<AtomPointerVector> queryMasks, targetMasks;

    // every time a new target is added, this flag will be set so we will know
    // to update proximity grids required for the search
    bool updateGrids;

    // grid spacing for ProximitySearch object
    mstreal gridSpacing;

    RMSDCalculator RC;
};

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Implements the FASSA (FAst Structure Search Algorithm). Options:");
  op.addOption("q", "query PDB file.", true);
  op.addOption("d", "a database file with a list of PDB files.", true);
  op.addOption("r", "RMSD cutoff.", true);
  op.addOption("min", "min number of matches.");
  op.addOption("max", "max number of matches.");
  op.setOptions(argc, argv);

  FASST S;
  cout << "Building the database..." << endl;
  S.setQuery(op.getString("q"));
  S.setMemorySaveMode(true);
  S.addTargets(MstUtils::fileToArray(op.getString("d")));
  S.setRMSDCutoff(op.getReal("r"));
  S.setMaxNumMatches(op.getInt("max", -1));
  S.setMinNumMatches(op.getInt("min", -1));
  S.setMaxGap(0, 1, 10); S.setMinGap(0, 1, 0);
  if (false) {
    int N = 10;
    for (int i = 0; i < N; i++) {
      mstreal d = i*20.0/(N-1) + (N-i-1)*5.0/(N-1);
      cout << "Searching with grid spacing " << d << "..." << endl;
      S.setGridSpacing(d);
      S.search();
    }
  } else {
    S.search();
  }
  cout << "found " << S.numMatches() << " matches:" << endl;
  set<fasstSolution> matches = S.getMatches(); int i = 0;
  for (auto it = matches.begin(); it != matches.end(); ++it, ++i) {
    cout << *it << endl;
    Structure match = S.getMatchStructure(*it, true, FASST::matchType::WITHGAPS);
    match.writePDB("/tmp/match" + MstUtils::toString(i) + ".pdb");
  }
}

FASST::FASST() {
  recLevel = 0;
  setRMSDCutoff(1.0);
  setSearchType(searchType::FULLBB);
  querySize = 0;
  updateGrids = false;
  gridSpacing = 15.0;
  memSave = false;
  maxNumMatches = minNumMatches = -1;
  gapConstSet = false;
}

FASST::~FASST() {
  // need to delete atoms only on the lowest level of recursion, because at
  // higher levels we point to the same atoms
  if (targetMasks.size()) targetMasks.back().deletePointers();
  for (int i = 0; i < targetStructs.size(); i++) delete targetStructs[i];
}

void FASST::setCurrentRMSDCutoff(mstreal cut) {
  rmsdCut = cut;
  residualCut = cut*cut*querySize;
}

void FASST::setMaxNumMatches(int _max) {
  maxNumMatches = _max;
  if ((minNumMatches >= 0) && (minNumMatches > maxNumMatches)) MstUtils::error("invalid pairing of min and max number of matches: " + MstUtils::toString(minNumMatches) + " and " + MstUtils::toString(maxNumMatches), "FASST::setMaxNumMatches");
}

void FASST::setMinNumMatches(int _min) {
  minNumMatches = _min;
  if ((maxNumMatches >= 0) && (minNumMatches > maxNumMatches)) MstUtils::error("invalid pairing of min and max number of matches: " + MstUtils::toString(minNumMatches) + " and " + MstUtils::toString(maxNumMatches), "FASST::setMinNumMatches");
}

void FASST::setMinGap(int i, int j, int gapLim) {
  if (gapLim < 0) MstUtils::error("gap constraints must be defined in the positive direction", "FASST::setMinGap");
  minGap[i][j] = gapLim;
  minGapSet[i][j] = true;
  gapConstSet = true;
}

void FASST::setMaxGap(int i, int j, int gapLim) {
  if (gapLim < 0) MstUtils::error("gap constraints must be defined in the positive direction", "FASST::setMinGap");
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

  // center-to-center distances in the qyery
  centToCentDist.resize(query.size(), vector<mstreal>(query.size(), 0));
  for (int i = 0; i < query.size(); i++) {
    for (int j = i+1; j < query.size(); j++) {
      centToCentDist[i][j] = query[i].getGeometricCenter().distancenc(query[j].getGeometricCenter());
      centToCentDist[j][i] = centToCentDist[i][j];
    }
  }

  // set query masks (subsets of query segments involved at each recursion level)
  queryMasks.resize(query.size());
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
  if (memSave) stripSidechains(*targetStruct);
  targetStructs.push_back(targetStruct);
  targets.push_back(AtomPointerVector());
  AtomPointerVector& target = targets.back();
  // we don't care about the chain topology of the target, so append all residues
  for (int i = 0; i < targetStruct->chainSize(); i++) {
    parseChain(targetStruct->getChain(i), target);
  }
  MstUtils::assert(target.size() > 0, "empty target in file '" + pdbFile + "'", "FASST::setTarget");

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

bool FASST::parseChain(const Chain& C, AtomPointerVector& searchable) {
  bool foundAll = true;
  for (int i = 0; i < C.residueSize(); i++) {
    Residue& res = C.getResidue(i);
    switch(type) {
      case searchType::CA:
      case searchType::FULLBB: {
        AtomPointerVector bb;
        for (int k = 0; k < searchableAtomTypes.size(); k++) {
          Atom* a = res.findAtom(searchableAtomTypes[k], false);
          if (a == NULL) { foundAll = false; break; }
          bb.push_back(a);
        }
        if (bb.size() == searchableAtomTypes.size()) searchable.insert(searchable.end(), bb.begin(), bb.end());
        break;
      }
      default:
        MstUtils::error("uknown search type '" + MstUtils::toString(type) + "' specified", "FASST::parseStructure");
    }
  }
  return foundAll;
}

void FASST::setSearchType(searchType _searchType) {
  type = _searchType;
  switch(type) {
    case searchType::CA:
      searchableAtomTypes =  {"CA"};
      break;
    case searchType::FULLBB:
      searchableAtomTypes =  {"N", "CA", "C", "O"};
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
        Atom& A = R[k];
        if (A.isNamed("N") || A.isNamed("CA") || A.isNamed("C") || A.isNamed("O")) continue;
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
  ps.clear(); ps.resize(query.size());
  for (int i = 0; i < query.size(); i++) {
    ps[i] = ProximitySearch(xlo, ylo, zlo, xhi, yhi, zhi, N);
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
    ps[i].dropAllPoints();
    AtomPointerVector& seg = query[i];
    int Na = atomToResIdx(target.size()) - atomToResIdx(seg.size()) + 1; // number of possible alignments
    segmentResiduals[i].resize(Na);
    for (int j = 0; j < Na; j++) {
      // NOTE: can save on this in several ways:
      // 1. the centroid calculation is effectively already done inside RMSDCalculator::bestRMSD
      // 2. updating just one atom involves a simple centroid adjustment, rather than recalculation
      // 3. is there a speedup to be gained from re-calculating RMSD with one atom updated only?
      AtomPointerVector targSeg = target.subvector(resToAtomIdx(j), resToAtomIdx(j) + query[i].size());
      segmentResiduals[i][j] = RC.bestResidual(query[i], targSeg);
      targSeg.getGeometricCenter(xc, yc, zc);
      ps[i].addPoint(xc, yc, zc, j);
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
  currRemBound = boundOnRemainder(true);

  // initialize center-to-center tolerances
  ccTol.clear(); ccTol.resize(query.size());
  for (int i = 0; i < query.size(); i++) {
    ccTol[i].resize(query.size(), vector<mstreal>(query.size(), -1.0)); // unnecessary (unused) entries will stay as -1
  }

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
}

mstreal FASST::boundOnRemainder(bool compute) {
  if (compute) {
    currRemBound = 0;
    for (int i = recLevel + 1; i < query.size(); i++) currRemBound += remOptions[recLevel][i].bestCost();
  }
  return currRemBound;
}

mstreal FASST::centToCentTol(int i, int j, bool recomputeResidual, bool recomputeBound) {
  mstreal remRes = residualCut - currentAlignmentResidual(recomputeResidual) - boundOnRemainder(recomputeBound);
  if (remRes < 0) return -1.0;
  return sqrt((remRes * (query[i].size() + query[j].size())) / (query[i].size() * query[j].size()));
}

void FASST::search() {
  if (isMinNumMatchesSet()) setCurrentRMSDCutoff(999.0);
  else setCurrentRMSDCutoff(rmsdCutRequested);

  solutions.clear();
  auto begin = chrono::high_resolution_clock::now();
  int prepTime = 0;
  for (int currentTarget = 0; currentTarget < targets.size(); currentTarget++) {
    auto beginPrep = chrono::high_resolution_clock::now();
    prepForSearch(currentTarget);
    auto endPrep = chrono::high_resolution_clock::now();
    prepTime += chrono::duration_cast<std::chrono::microseconds>(endPrep-beginPrep).count();
    AtomPointerVector& target = targets[currentTarget];
    vector<int> okLocations, badLocations;
    okLocations.reserve(target.size()); badLocations.reserve(target.size());
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
      // fill up sub-alignment with target atoms
      int N = targetMasks[recLevel].size();
      int n = query[recLevel].size();
      int si = resToAtomIdx(currAlignment[recLevel]);
      for (int i = 0; i < n; i++) {
        targetMasks[recLevel][N - n + i]->setCoor(target[si + i]->getX(), target[si + i]->getY(), target[si + i]->getZ());
      }

      // 2. compute the total residual from the current alignment
      mstreal curBound = currentAlignmentResidual(true) + boundOnRemainder(true);
      if (curBound > residualCut) continue;

      // 3. update update remaining options for subsequent segments based on the
      // newly made choice. The set of options on the next recursion level is a
      // subset of the set of options on the previous level.
      int remSegs = query.size() - (recLevel + 1);
      if (remSegs > 0) {
        recLevel++;
        // copy remaining options from the previous recursion level. This way,
        // we can compute bounds on this level and can do set intersections to
        // further narrow this down
        for (int i = recLevel; i < query.size(); i++) remOptions[recLevel][i].copyIn(remOptions[recLevel-1][i]);
        // if any gap constraints exist, limit options at this recursion level accordingly
        if (gapConstraintsExist()) {
          bool levelExhausted = false;
          for (int j = 0; j < recLevel; j++) {
            for (int i = recLevel; i < query.size(); i++) {
              if (minGapSet[i][j]) remOptions[recLevel][i].constrainLE(currAlignment[j] - minGap[i][j]);
              if (maxGapSet[i][j]) remOptions[recLevel][i].constrainGE(currAlignment[j] - maxGap[i][j]);
              if (minGapSet[j][i]) remOptions[recLevel][i].constrainGE(currAlignment[j] + minGap[j][i]);
              if (maxGapSet[j][i]) remOptions[recLevel][i].constrainLE(currAlignment[j] + maxGap[j][i]);
              if (remOptions[recLevel][i].empty()) {
                recLevel--;
                levelExhausted = true;
                break;
              }
            }
            if (levelExhausted) break;
          }
          if (levelExhausted) continue;
        }
        mstreal dij, de, d, eps = 10E-8;
        for (int c = 0; true; c++) {
          bool updated = false, levelExhausted = false;
          for (int i = recLevel; i < query.size(); i++) {
            optList& remSet = remOptions[recLevel][i];
            for (int j = 0; j < recLevel; j++) {  // j runs over already placed segments
              de = centToCentTol(i, j);
              if (de < 0) {
                recLevel--;
                levelExhausted = true;
                break;
              }
              dij = centToCentDist[i][j];
              mstreal dePrev = ((c == 0) ? ccTol[recLevel-1][j][i] : ccTol[recLevel][j][i]);
              int numLocs = remSet.size();
              CartesianPoint& cj = ps[j].getPoint(currAlignment[j]); // the current centroid of segment j

              // If the set of options for the current segment was arrived at,
              // in part, by limiting center-to-center distances, then we will
              // tighten that list by removing options that are outside of the
              // range allowed at this recursion level. Otherwise, we will do a
              // general proximity search given the current tolerance and will
              // tighten the list that way.
              if (dePrev < 0) {
                okLocations.resize(0);
                ps[i].pointsWithin(cj, max(dij - de, 0.0), dij + de, &okLocations);
                remSet.intersectOptions(okLocations);
              } else {
                badLocations.resize(0);
                ps[i].pointsWithin(cj, max(dij - dePrev, 0.0), dij - de - eps, &badLocations);
                ps[i].pointsWithin(cj, dij + de + eps, dij + dePrev, &badLocations);
                for (int k = 0; k < badLocations.size(); k++) {
                  remSet.removeOption(badLocations[k]);
                }
              }
              ccTol[recLevel][j][i] = de;
              if (numLocs != remSet.size()) {
                // this both updates the bound and checks that there are still
                // feasible solutions left
                if ((remSet.empty()) || (currResidual + boundOnRemainder(true) > residualCut)) {
                  recLevel--;
                  levelExhausted = true;
                  break;
                }
                updated = true;
              }
            }
            if (levelExhausted) break;
          }
          if (levelExhausted) break;
          if (!updated) break;
        }
      } else {
        // if at the lowest recursion level already, then record the solution
        solutions.insert(fasstSolution(currAlignment, sqrt(currResidual/querySize), currentTarget, solutions.size()));
        // if ((maxNumMatches > 0) && (solutions.size() > maxNumMatches)) {
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
  auto end = chrono::high_resolution_clock::now();
  int searchTime = chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
  prepTime = prepTime/1000;
  cout << "prep time " << prepTime << " ms" << std::endl;
  cout << "total time " << searchTime << " ms" << std::endl;
  cout << "prep time was " << (100.0*prepTime/searchTime) << " % of the total" << std::endl;
}

mstreal FASST::currentAlignmentResidual(bool compute) {
  if (compute) {
    if (recLevel == 0) {
      currResidual = segmentResiduals[recLevel][currAlignment[recLevel]];
    } else {
      currResidual = RC.bestResidual(queryMasks[recLevel], targetMasks[recLevel]);
    }
  }
  return currResidual;
}

void FASST::getMatchStructure(const fasstSolution& sol, Structure& match, bool detailed, matchType type) {
  vector<Structure> matches;
  getMatchStructures(vector<fasstSolution>(1, sol), matches, detailed, type);
  match = matches[0];
}

Structure FASST::getMatchStructure(const fasstSolution& sol, bool detailed, matchType type) {
  Structure match; getMatchStructure(sol, match, detailed, type); return match;
}

void FASST::getMatchStructures(const vector<fasstSolution>& sols, vector<Structure>& matches, bool detailed, matchType type) {
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

  // flatten query into a single array of atoms
  AtomPointerVector queryAtoms;
  for (int k = 0; k < query.size(); k++) {
    for (int ai = 0; ai < query[k].size(); ai++) queryAtoms.push_back(query[k][ai]);
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
      // re-read structure
      if (targetSource[idx].type == targetFileType::PDB) {
        dummy.reset();
        dummy.readPDB(targetSource[idx].file, "QUIET");
        tr[idx].apply(dummy);
      } else {
        MstUtils::error("don't know how to read target of this type", "ASST::getMatchStructures");
      }
      targetStruct = &dummy;
    }

    // visit each solution from this target
    vector<int>& solIndices = it->second;
    for (int i = 0; i < solIndices.size(); i++) {
      int solIndex = solIndices[i];
      const fasstSolution& sol = sols[solIndex];
      Structure& match = matches[solIndex];
      vector<int> alignment = sol.getAlignment();
      if (alignment.size() != query.size()) {
        MstUtils::error("solution alignment size inconsistent with number of query segments", "FASST::getMatchStructures");
      }
      // walk over target alignment and pick out corresponding residues
      matchAtoms.resize(0);
      for (int k = 0; k < alignment.size(); k++) {
        // copy matching residues from the original target structure
        int si = resToAtomIdx(alignment[k]);
        int L = query[k].size();
        if (si + L > target.size()) {
          MstUtils::error("solution points to atom outside of target range", "FASST::getMatchStructures");
        }
        for (int ai = si; ai < si + L; ai++) matchAtoms.push_back(target[ai]);
        switch(type) {
          case matchType::REGION: {
            for (int ri = atomToResIdx(si); ri < atomToResIdx(si + L); ri++) {
              Residue* res = target[resToAtomIdx(ri)]->getResidue();
              if (reread) res = &(targetStruct->getResidue(res->getResidueIndex()));
              match.addResidue(res);
            }
            break;
          }
          case matchType::WITHGAPS: {
            // figure out the range of residues to excise from target structure
            vector<bool> toInclude(targetStruct->residueSize());
            for (int i = 0; i < query.size(); i++) {
              // each individual segments should be included
              for (int k = alignment[i]; k < alignment[i] + atomToResIdx(query[i].size()); k++) toInclude[k] = true;
              // then some gaps also
              for (int j = 0; j < query.size(); j++) {
                if (i == j) continue;
                // if gap constrained, fill in between these two segments
                if (gapConstrained(i, j)) {
                  if (alignment[i] > alignment[j]) MstUtils::error("solution not consistent with current gap constraints", "FASST::getMatchStructures");
                  for (int k = alignment[i] + atomToResIdx(query[i].size()); k < alignment[j]; k++) {
                    toInclude[k] = true;
                  }
                }
              }
            }
            for (int ri = 0; ri < toInclude.size(); ri++) {
              if (toInclude[ri]) match.addResidue(&(targetStruct->getResidue(ri)));
            }
            break;
          }
          case matchType::FULL:
            match = *targetStruct;
            break;
          default:
            MstUtils::error("unknown match output type", "FASST::getMatchStructures");
        }
      }
      // align matching region onto query, transforming the match itself
      rc.align(matchAtoms, queryAtoms, match);
    }
  }
}
