#ifndef _MSTFASST_H
#define _MSTFASST_H

#include "msttypes.h"
#include "msttransforms.h"
#include "mstsequence.h"
#include <list>
#include <chrono>

using namespace MST;

class fasstSolution {
  public:
    fasstSolution(const vector<int>& _alignment, mstreal _rmsd, int _target, int _foundOrder, vector<int> segOrder = vector<int>()) {
      alignment = _alignment; rmsd = _rmsd; targetIndex = _target; foundOrder = _foundOrder;
      if (!segOrder.empty()) {
        for (int i = 0; i < _alignment.size(); i++) alignment[segOrder[i]] = _alignment[i];
      }
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
    enum targetFileType { PDB = 1, BINDATABASE, STRUCTURE };

    class targetInfo {
      public:
        targetInfo(const string& _file, targetFileType _type, int _loc, bool _memSave) {
          file = _file; type = _type; loc = _loc; memSave = _memSave;
        }
        targetInfo(const targetInfo& I) {
          file = I.file; type = I.type; loc = I.loc; memSave = I.memSave;
        }
        string file;         // source file
        targetFileType type; // file type (PDB or database)
        int loc;             // location info within file
        bool memSave;        // was the target read with memory save on?
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
        // after this, the only remaining allowed options will be those that were
        // allowed before AND are in the specified list
        void intersectOptions(const vector<int>& opts);
        // keep only options with indices less than or equal to the given one
        void constrainLE(int idx);
        // keep only options with indices greater than or equal to the given one
        void constrainGE(int idx);
        // keep only options with indices in the range [idxLow, idxHigh]
        void constrainRange(int idxLow, int idxHigh);
        // wipe all options
        void removeAllOptions();
        // copy the valid options from a different set of the same options
        void copyIn(const optList& opt);
        void removeOptions(int b, int e);
        void removeOption(int k);
        mstreal bestCost() { return (numIn == 0) ? 0.0 : costs[bestCostRank]; }
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

    /* TODO: add getMatchSequence and getMatchSequences, which return Sequence or vector<Sequence>
     * TODO: pre-center queryMasks at each recursion level (means allocating separate space for each level)
     * TODO: Jianfu and Craig will look for fast ways of NN searches in Hamming distnce space
     * TODO: goal is to define redundancy as follows:
     * 1. define the redundancy of each segment alignment separately
     * 2. if a segment is shorter than L, expand it to L (L = 30)
     * 3. compute RMSD between these L residues. If it is above 1.5 A, then the
     *    two alignments are different, so done. If it is below, then do an
     *    ungapped sequence alignment and judge by the score.
     * 4. some windows will not be expandable (hit a chain terminus), in which
     *    case we will compare the common portion of any two windows. Both the
     *    RMSD cutoff and the sequence identity cutoff have to scale.
     * TODO: need to define chain ends. Do so via the standard information.
     * TODO: enable checking for continuity of gaps */

    ~FASST();
    FASST();
    void setQuery(const string& pdbFile);
    void setQuery(const Structure& Q);
    Structure getQuery() { return queryStruct; }
    void addTarget(const Structure& T);
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
    void setSufficientNumMatches(int _suff);
    void unsetSufficientNumMatches() { suffNumMatches = -1; }
    int getMaxNumMatches() { return maxNumMatches; }
    int getMinNumMatches() { return minNumMatches; }
    int getSufficientNumMatches() { return suffNumMatches; }
    bool isMaxNumMatchesSet() { return (maxNumMatches > 0); }
    bool isMinNumMatchesSet() { return (minNumMatches > 0); }
    bool isSufficientNumMatchesSet() { return (suffNumMatches > 0); }
    set<fasstSolution> getMatches() { return solutions; }
    string toString(const fasstSolution& sol);

    bool gapConstrained(int i, int j) { return (minGapSet[i][j] || maxGapSet[i][j]); }
    void setMinGap(int i, int j, int gapLim); // target topology: [segment i] [gap of at list gapLim long] [segment j]
    void setMaxGap(int i, int j, int gapLim); // target topology: [segment i] [gap of at most gapLim long] [segment j]
    void resetGapConstraints();
    bool gapConstraintsExist() { return gapConstSet; }
    bool validateSearchRequest(); // make sure all user specified requirements are consistent
    void getMatchStructure(const fasstSolution& sol, Structure& match, bool detailed = false, matchType type = matchType::REGION);
    Structure getMatchStructure(const fasstSolution& sol, bool detailed = false, matchType type = matchType::REGION);
    void getMatchStructures(const set<fasstSolution>& sols, vector<Structure>& matches, bool detailed = false, matchType type = matchType::REGION);
    void getMatchStructures(const vector<fasstSolution>& sols, vector<Structure>& matches, bool detailed = false, matchType type = matchType::REGION);
    void writeDatabase(const string& dbFile);
    void readDatabase(const string& dbFile);

    // TODO:
    // results caching
    int cacheSearch();            // save current search results (needs to also save the query)
    void recallFromCache(int i);  // overwrite the current results with those from the cache (solutions, query, ...)
    // TODO: name match structures; match.setName(targetStruct->getName()) on line 768

  protected:
    void processQuery();
    void setCurrentRMSDCutoff(mstreal cut);
    void prepForSearch(int ti);
    bool parseChain(const Chain& S, AtomPointerVector& searchable, Sequence* seq = NULL);
    mstreal currentAlignmentResidual(bool compute);   // computes the accumulated residual up to and including segment recLevel
    mstreal boundOnRemainder(bool compute);           // computes the lower bound expected from segments recLevel+1 and on
    // void updateQueryCentroids();                      // assumes that appropriate transformation matrices were previously set with a call to currentAlignmentResidual(true)
    int resToAtomIdx(int resIdx) { return resIdx * atomsPerRes; }
    int atomToResIdx(int atomIdx) { return atomIdx / atomsPerRes; }
    mstreal centToCentTol(int i);
    mstreal segCentToPrevSegCentTol(int i);
    void rebuildProximityGrids();
    void addTargetStructure(Structure* targetStruct);
    void stripSidechains(Structure& S);
    bool areNumMatchConstraintsConsistent();

  private:
    vector<Structure*> targetStructs;
    vector<AtomPointerVector> targets;       // target structures (just the parts that will be searched over)
    vector<Sequence> targSeqs;               // target sequences (of just the parts that will be searched over)
    vector<targetInfo> targetSource;         // where each target was read from (in case need to re-read it)
    vector<int> targChainBeg, targChainEnd;  // targChainBeg[i] and targChainEnd[i] contain the chain start and end indices for the chain
                                             // that contains the residue with index i (in the overal fused sequence)
    vector<Transform> tr;                    // transformations from the original frame to the common frames of reference for each target
    int currentTarget;                       // the index of the target currently being searched for

    Structure queryStruct;
    vector<AtomPointerVector> queryOrig;     // just the part of the query that will be sought, split by segment
    vector<AtomPointerVector> query;         // same as above, but with segments re-orderd for optimal searching
    vector<int> qSegOrd;                     // qSegOrd[i] is the index (in the original queryOrig) of the i-th segment in query
    mstreal xlo, ylo, zlo, xhi, yhi, zhi;    // bounding box of the search database
    bool memSave;                            // save memory by storing only the backbone of targets?

    vector<vector<int> > minGap, maxGap;     // minimum and maximum sequence separations allowed between each pair of segments
    vector<vector<bool> > minGapSet, maxGapSet;
    bool gapConstSet;

    // segmentResiduals[i][j] is the residual of the alignment of segment i, in which
    // its starting residue aligns with the residue index j in the target
    vector<vector<mstreal> > segmentResiduals;

    int recLevel;                            // segments up to index recLevel are already placed
    int atomsPerRes;
    searchType type;
    vector<vector<string> > searchableAtomTypes;
    int querySize;
    int maxNumMatches, minNumMatches, suffNumMatches;

    // remOptions[L][i] is a set of alignments for segment i (i >= L) at recursion level
    // L, which are stored sorted by their own residual, through the optList
    // data structure. Note that segments 0 through L-1 have already been place
    // at recursion level L.
    vector<vector<optList> > remOptions;

    // // ccTol[L][i][j], j > i, is the acceptable tolerance (the delta) on the
    // // center-to-center distance between segments i and j at recursion level L
    // vector<vector<vector<mstreal> > > ccTol;

    // alignment indices for segments visited up to the current recursion level
    vector<int> currAlignment;

    // the residual of the above alignment, computed and stored
    mstreal currResidual;

    // the same residuals for each recursion level
    vector<mstreal> currResiduals;

    // the centroids of the currently aligned portion, at each recursion level
    vector<CartesianPoint> currCents;

    // the bound on the parts remaining to align, computed and stored
    mstreal currRemBound;

    // the distance between the centroid of each segment and the centroid of the
    // "previous" segment, in the order in which they will be placed
    vector<mstreal> segCentToPrevSegCentDist;

    // distances between the centroid of the already placed sub-query and every
    // sub-sequence segment. I.e., centToCentDist[L][i] is the distance between
    // the centroid of the sub-query placed at recursion level L (i.e., segments
    // 0 through L) and some sub-sequent segment i (i > L)
    vector<vector<mstreal> > centToCentDist;

    // set of solutions, sorted by RMSD
    set<fasstSolution> solutions;

    // ProximitySearch for finding nearby centroids of target segments (there
    // will be one ProximitySearch object per query segment)
    vector<ProximitySearch*> ps;

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

#endif
