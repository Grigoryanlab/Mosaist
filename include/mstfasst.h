#ifndef _MSTFASST_H
#define _MSTFASST_H

#include "msttypes.h"
#include "msttransforms.h"
#include "mstsequence.h"
#include <list>
#include <chrono>

using namespace MST;

class fasstSolutionSet;
class fasstSolution {
  friend class fasstSolutionSet;
  public:
    fasstSolution() { rmsd = 0.0; context = NULL; }
    fasstSolution(const vector<int>& _alignment, mstreal _rmsd, int _target, const Transform& _tr, const vector<int>& _segLengths, vector<int> segOrder = vector<int>());
    fasstSolution(const fasstSolution& _sol);
    ~fasstSolution() { if (context != NULL) delete context; }

    mstreal getRMSD() const { return rmsd; }
    int getTargetIndex() const { return targetIndex; }
    vector<int> getAlignment() const { return alignment; }
    vector<int> getSegLengths() const { return segLengths; }
    int segLength(int i) const { return segLengths[i]; }
    int numSegments() const { return alignment.size(); }
    int operator[](int i) const { return alignment[i]; }
    Transform getTransform() const { return tr; }
    bool seqContextDefined() const { return ((context != NULL) && (!context->segSeq.empty())); }
    void setSeqContext(const vector<Sequence>& _segSeq, const vector<Sequence>& _nSeq, const vector<Sequence>& _cSeq);
    void setStructContext(const vector<AtomPointerVector>& _segStr, const vector<AtomPointerVector>& _nStr, const vector<AtomPointerVector>& _cStr);
    void setRMSD(mstreal r) { rmsd = r; }
    void setTransform(const Transform& _tr) { tr = _tr; }

    static bool foundBefore(const fasstSolution* solA, const fasstSolution* solB) {
      if (solA->targetIndex != solB->targetIndex) return solA->targetIndex < solB->targetIndex;
      if (solA->alignment.size() != solB->alignment.size()) MstUtils::error("cannot compare solutions of different topology", "fasstSolution::foundBefore");
      for (int i = 0; i < solA->alignment.size(); i++) {
        if (solA->alignment[i] != solB->alignment[i]) return solA->alignment[i] < solB->alignment[i];
      }
      return false;
    }
    friend bool operator<(const fasstSolution& si, const fasstSolution& sj) {
      if (si.rmsd != sj.rmsd) return (si.rmsd < sj.rmsd);
      if (si.targetIndex != sj.targetIndex) return (si.targetIndex < sj.targetIndex);
      if (si.alignment.size() != sj.alignment.size()) MstUtils::error("cannot compare solutions of different topology", "fasstSolution::operator<");
      for (int k = 0; k < si.alignment.size(); k++) {
        if (si.alignment[k] != sj.alignment[k]) return (si.alignment[k] < sj.alignment[k]);
      }
      // this will only happen if the two solutioins are identical: have the same
      // RMSD, come from the same target, and signify exactly the same alignment
      return false;
    }

    void write(ostream& _os) const;
    void read(istream& _is);

    friend ostream& operator<<(ostream &_os, const fasstSolution& _sol) {
      _os << std::setprecision(6) << std::fixed << _sol.rmsd << " " << _sol.targetIndex << " [" << MstUtils::vecToString(_sol.alignment, ", ") << "]";
      return _os;
    }

  protected:
    const vector<Sequence>& segmentSeqs() const { return context->segSeq; }
    const vector<Sequence>& nTermContext() const { return context->nSeq; }
    const vector<Sequence>& cTermContext() const { return context->cSeq; }

  private:
    vector<int> alignment, segLengths;
    mstreal rmsd;
    int targetIndex;
    Transform tr; // how the target needs to be transformed for this match to optimally fit onto the query

    class solContext {
      public:
        solContext() {}
        solContext(const solContext& sol) {
          segSeq = sol.segSeq; nSeq = sol.nSeq; cSeq = sol.cSeq;
          segStr = sol.segStr; nStr = sol.nStr; cStr = sol.cStr;
          for (int i = 0; i < sol.segStr.size(); i++) segStr[i] = sol.segStr[i].clone();
          for (int i = 0; i < sol.nStr.size(); i++) nStr[i] = sol.nStr[i].clone();
          for (int i = 0; i < sol.cStr.size(); i++) cStr[i] = sol.cStr[i].clone();
        }
        ~solContext() {
          for (int i = 0; i < segStr.size(); i++) segStr[i].deletePointers();
          for (int i = 0; i < nStr.size(); i++) nStr[i].deletePointers();
          for (int i = 0; i < cStr.size(); i++) cStr[i].deletePointers();
        }
        void write(ostream &_os) const { // write solContext to a binary stream
          MstUtils::writeBin(_os, (int) segSeq.size());
          for (int i = 0; i < segSeq.size(); i++) {
            segSeq[i].write(_os); nSeq[i].write(_os); cSeq[i].write(_os);
          }
          MstUtils::writeBin(_os, (int) segStr.size());
          for (int i = 0; i < segStr.size(); i++) {
            segStr[i].write(_os); nStr[i].write(_os); cStr[i].write(_os);
          }
        }
        void read(istream &_is) { // read solContext from a binary stream
          int len; MstUtils::readBin(_is, len);
          segSeq.resize(len); nSeq.resize(len); cSeq.resize(len);
          for (int i = 0; i < len; i++) {
            segSeq[i].read(_is); nSeq[i].read(_is); cSeq[i].read(_is);
          }
          MstUtils::readBin(_is, len);
          segStr.resize(len); nStr.resize(len); cStr.resize(len);
          for (int i = 0; i < len; i++) {
            segStr[i].read(_is); nStr[i].read(_is); cStr[i].read(_is);
          }
        }

        vector<Sequence> segSeq, nSeq, cSeq; // sequence of each segment and N- and C-terminal contexts
        vector<AtomPointerVector> segStr, nStr, cStr; // structure of each segment and N- and C-terminal contexts
    };
    solContext* context;
};

class fasstSolutionSet {
  public:
    fasstSolutionSet() { updated = false; }
    fasstSolutionSet(const fasstSolution& sol);
    fasstSolutionSet(const fasstSolutionSet& sols);
    fasstSolutionSet& operator=(const fasstSolutionSet& sols);
    bool insert(const fasstSolution& sol, mstreal redundancyCut = 1); // returns whether the insert was performed
    set<fasstSolution>::iterator begin() const { return solsSet.begin(); }
    set<fasstSolution>::reverse_iterator rbegin() const { return solsSet.rbegin(); }
    set<fasstSolution>::iterator end() const { return solsSet.end(); }
    set<fasstSolution>::reverse_iterator rend() const { return solsSet.rend(); }
    set<fasstSolution>::iterator erase(const set<fasstSolution>::iterator it) { return solsSet.erase(it); updated = true; }
    // const fasstSolution& operator[] (int i) const;
    fasstSolution& operator[] (int i);
    int size() const { return solsSet.size(); }
    void clear() { solsSet.clear(); updated = true; }
    mstreal worstRMSD() { return (solsSet.rbegin())->getRMSD(); }
    mstreal bestRMSD() { return (solsSet.begin())->getRMSD(); }
    vector<fasstSolution*> orderByDiscovery();

    void write(ostream &_os) const; // write fasstSolutionSet to a binary stream
    void read(istream &_os);  // read fasstSolutionSet from a binary stream

    friend ostream& operator<<(ostream &_os, const fasstSolutionSet& _sols) {
      for (auto it = _sols.solsSet.begin(); it != _sols.solsSet.end(); ++it) {
        _os << *it;
      }
      return _os;
    }

  protected:
    // decides whether numID identities within an alignment of length numTot
    // meets the identity cutoff cut established for alignment length L0
    bool isWithinSeqID(int L0, mstreal cut, int numTot, int numID);

  private:
    // An important point about std::set is that it never invalidates pointers
    // to its elements (i.e., it does not copy them upon resizing like vector).
    set<fasstSolution> solsSet;
    vector<fasstSolution*> solsVec;
    bool updated;
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

    ~FASST();
    FASST();
    void setQuery(const string& pdbFile, bool autoSplitChains = true);
    void setQuery(const Structure& Q, bool autoSplitChains = true);
    Structure getQuery() const { return queryStruct; }
    AtomPointerVector getQuerySearchedAtoms() const;
    void addTarget(const Structure& T);
    void addTarget(const string& pdbFile);
    void addTargets(const vector<string>& pdbFiles);
    void stripSidechains(Structure& S);

    void addResidueProperties(int ti, const string& propType, const vector<mstreal>& propVals);
    void addResiduePairProperties(int ti, const string& propType, const map<int, map<int, mstreal> >& propVals);
    // void addResidueRelationships(int ti, const string& propType, const map<int, map<int, map<int, set<int> > >& resRels);
    void addResidueRelationship(int ti, const string& propType, int ri, int tj, int rj);
    mstreal isResiduePropertyPopulated(const string& propType);
    bool hasResidueProperties(int ti, const string& propType, int ri);
    mstreal getResidueProperty(int ti, const string& propType, int ri);
    bool hasResiduePairProperties(int ti, const string& propType, int ri);
    mstreal isResiduePairPropertyPopulated(const string& propType);
    map<int, mstreal> getResiduePairProperties(int ti, const string& propType, int ri);

    bool isPropertyDefined(const string& propType);
    bool isPropertyDefined(const string& propType, int ti);
    int numTargets() const { return targetStructs.size(); }
    Structure getTargetCopy(int i) const { return *(targetStructs[i]); }
    Structure* getTarget(int i) { return targetStructs[i]; }
    void setRMSDCutoff(mstreal cut) { rmsdCutRequested = cut; }
    void setSearchType(searchType _searchType);
    void setMemorySaveMode(bool _memSave) { memSave = _memSave; }
    void setGridSpacing(mstreal _spacing) { gridSpacing = _spacing; updateGrids = true; }
    fasstSolutionSet search();
    int numMatches() { return solutions.size(); }
    void setMaxNumMatches(int _max);
    void setMinNumMatches(int _min);
    void unsetMaxNumMatches() { maxNumMatches = -1; }
    void unsetMinNumMatches() { minNumMatches = -1; }
    void setSufficientNumMatches(int _suff);
    void unsetSufficientNumMatches() { suffNumMatches = -1; }
    mstreal getRMSDCutoff() { return rmsdCutRequested; }
    int getMaxNumMatches() { return maxNumMatches; }
    int getMinNumMatches() { return minNumMatches; }
    int getSufficientNumMatches() { return suffNumMatches; }
    bool isMaxNumMatchesSet() { return (maxNumMatches > 0); }
    bool isMinNumMatchesSet() { return (minNumMatches > 0); }
    bool isSufficientNumMatchesSet() { return (suffNumMatches > 0); }
    fasstSolutionSet getMatches() { return solutions; }
    string toString(const fasstSolution& sol);

    bool gapConstrained(int i, int j) { return (minGapSet[i][j] || maxGapSet[i][j]); }
    void setMinGap(int i, int j, int gapLim); // target topology: [segment i] [gap of at list gapLim long] [segment j]
    void setMaxGap(int i, int j, int gapLim); // target topology: [segment i] [gap of at most gapLim long] [segment j]
    void resetGapConstraints();
    bool gapConstraintsExist() { return gapConstSet; }
    bool validateSearchRequest(); // make sure all user specified requirements are consistent
    void writeDatabase(const string& dbFile);
    void readDatabase(const string& dbFile);

    // get various match properties
    void getMatchStructure(const fasstSolution& sol, Structure& match, bool detailed = false, matchType type = matchType::REGION, bool algn = true);
    Structure getMatchStructure(const fasstSolution& sol, bool detailed = false, matchType type = matchType::REGION, bool algn = true);
    void getMatchStructures(fasstSolutionSet& sols, vector<Structure>& matches, bool detailed = false, matchType type = matchType::REGION, bool algn = true);
    vector<Sequence> getMatchSequences(fasstSolutionSet& sols, matchType type = matchType::REGION);
    Sequence getMatchSequence(const fasstSolution& sol, matchType type = matchType::REGION);
    vector<vector<mstreal> > getResidueProperties(fasstSolutionSet& sols, const string& propType, matchType type = matchType::REGION);
    vector<mstreal> getResidueProperties(const fasstSolution& sol, const string& propType, matchType type = matchType::REGION);
    vector<int> getMatchResidueIndices(const fasstSolution& sol, matchType type = matchType::REGION); // figure out the range of residues to excise from target structure

    /* Computes the RMSD of the given match to a query that is (possibly)
     * different from the one it was found with (if it even came from a search). */
    // mstreal matchRMSD(const fasstSolution& sol, const AtomPointerVector& query);
    vector<mstreal> matchRMSDs(fasstSolutionSet& sols, const AtomPointerVector& query, bool update = false);

    /* Normally, the redundancy cutoff is between 0 and 1. But one can set it to
     * values outside of this range, in principle. Setting it to a value above 1
     * will cause no redundancy cutoff to be applied, but will populate solution
     * objects with sequence context information, in case (for example) a filter
     * for redundancy will need to be applied later. */
    void pruneRedundancy(mstreal _redundancyCut = 0.5) { redundancyCut = _redundancyCut; }
    bool isRedundancyCutSet() { return redundancyCut < 1; }
    mstreal getRedundancy() { return redundancyCut; }

  protected:
    void processQuery();
    void setCurrentRMSDCutoff(mstreal cut);
    void prepForSearch(int ti);
    bool parseChain(const Chain& S, AtomPointerVector* searchable = NULL, Sequence* seq = NULL);
    mstreal currentAlignmentResidual(bool compute, bool setTransform = false);   // computes the accumulated residual up to and including segment recLevel
    mstreal boundOnRemainder(bool compute);           // computes the lower bound expected from segments recLevel+1 and on
    Transform currentTransform();                     // tansform for the alignment corresponding to the current residual
    // void updateQueryCentroids();                      // assumes that appropriate transformation matrices were previously set with a call to currentAlignmentResidual(true)
    int resToAtomIdx(int resIdx) { return resIdx * atomsPerRes; }
    int atomToResIdx(int atomIdx) { return atomIdx / atomsPerRes; }
    mstreal centToCentTol(int i);
    mstreal segCentToPrevSegCentTol(int i);
    void rebuildProximityGrids();
    void addTargetStructure(Structure* targetStruct);
    bool areNumMatchConstraintsConsistent();
    void addSequenceContext(fasstSolution& sol); // decorate the solution with sequence context

  private:
    /* targetStructs[i] and targets[i] store the original i-th target structure
     * and just the part of it that will be searched over, respectively. NOTE:
     * Atoms* in targets[i] point to Atoms of targetStructs[i]. So there is only
     * one compy of the target structure and there is no amboguity in mapping
     * targets[i] to targetStructs[i], even if when building tagrets[i] some of
     * the residues in targetStructs[i] had to be ignored for whatever reason
     * (e.g., missing backbone). */
    vector<Structure*> targetStructs;
    vector<AtomPointerVector> targets;

    vector<Sequence> targSeqs;               // target sequences (of just the parts that will be searched over)
    vector<targetInfo> targetSource;         // where each target was read from (in case need to re-read it)
    vector<int> targChainBeg, targChainEnd;  // targChainBeg[i] and targChainEnd[i] contain the chain start and end indices for the chain
                                             // that contains the residue with index i (in the overal concatenated sequence). Residue indices
                                             // are based on just the portion of the structure to be searched over.

    /* Object for holding real-valued residue properties. Specifically,
     * resProperties["env"][ti][ri] is the value of the "env" property for
     * residue ri in target with index ti. */
    map<string, map<int, vector<mstreal> > > resProperties;

    /* Object for holding real-valued residue-pair properties. Specifically,
     * resPairProperties["cont"][ti][ri][rj] is the value of the "cont" property
     * (e.g., contact degree) between residues ri and rj in target with index ti.
     * NOTE: this property can be directional (i.e., pairs are not mirrored). */
    map<string, map<int, map<int, map<int, mstreal> > > > resPairProperties;

    /* Object for holding residue-pair relational graphs. Specifically,
     * resRelProperties["sim"][ti][ri][tj] is the set of all residues in target
     * tj that are related by the property "sim" to the residue ri in target ti.
     * NOTE: this property can be directional (i.e., relationships are not mirrored). */
    map<string, map<int, map<int, map<int, set<int> > > > > resRelProperties;

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

    // redundancy options
    int contextLength;                       // how long of a local window to consider when comparing segment alingments between solutions
    mstreal redundancyCut;                   // maximum sequence identity (as a fraction), defined over contextLength-long windows

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
    fasstSolutionSet solutions;

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
