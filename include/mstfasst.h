#ifndef _MSTFASST_H
#define _MSTFASST_H

#include "msttypes.h"
#include "msttransforms.h"
#include "mstsequence.h"
#include <list>
#include <chrono>
#include <limits.h>

using namespace MST;

class fasstSolutionSet;
class fasstSolutionAddress {
  public:
    fasstSolutionAddress() { targetIndex = -1; }
    fasstSolutionAddress(int ti, const vector<int>& al) { targetIndex = ti; alignment = al; }
    void write(ostream& _os) const;
    void read(istream& _is);

    friend ostream& operator<<(ostream &_os, const fasstSolutionAddress& addr) {
      _os << addr.targetIndex << " [" << MstUtils::vecToString(addr.alignment, ", ") << "]";
      return _os;
    }

    friend bool operator==(const fasstSolutionAddress& ai, const fasstSolutionAddress& aj) {
      return (ai.targetIndex == aj.targetIndex) && (ai.alignment == aj.alignment);
    }

    int targetIndex;
    vector<int> alignment;
};

class fasstSolution {
  friend class fasstSolutionSet;
  public:
    class resAddress {
      public:
        resAddress(int ti = 0, int ri = 0) {
          if ((ti < 0) || (ri < 0)) MstUtils::error("indices cannot be negative!", "fasstSolution::resAddress::resAddress(int, int)");
          if ((ti > USHRT_MAX) || (ri > USHRT_MAX)) MstUtils::error("residue address indices (" + MstUtils::toString(ti) + ", " + MstUtils::toString(ri) + ") out of range for a short, recompile with int!", "fasstSolution::resAddress::resAddress(int, int)");
          targIdx = ti; resIdx = ri;
        }
        unsigned short& targIndex() { return targIdx; }
        unsigned short& resIndex() { return resIdx; }
        friend bool operator<(const resAddress& ai, const resAddress& aj) {
          if (ai.targIdx != aj.targIdx) return (ai.targIdx < aj.targIdx);
          return (ai.resIdx < aj.resIdx);
        }
        friend bool operator>(const resAddress& ai, const resAddress& aj) {
          if (ai.targIdx != aj.targIdx) return (ai.targIdx > aj.targIdx);
          return (ai.resIdx > aj.resIdx);
        }
        friend bool operator==(const resAddress& ai, const resAddress& aj) {
          return (ai.targIdx == aj.targIdx) && (ai.resIdx == aj.resIdx);
        }

      private:
        unsigned short targIdx, resIdx;
    };
    // typedef pair<short,short> resAddress;

    fasstSolution() { rmsd = 0.0; context = NULL; }
    fasstSolution(const vector<int>& _alignment, mstreal _rmsd, int _target, const Transform& _tr, const vector<int>& _segLengths, vector<int> segOrder = vector<int>());
    fasstSolution(const fasstSolution& _sol);
    fasstSolution(const fasstSolutionAddress& addr, const vector<int> segLen);
    ~fasstSolution() { if (context != NULL) delete context; }

    mstreal getRMSD() const { return rmsd; }
    int getTargetIndex() const { return targetIndex; }
    vector<int> getAlignment() const { return alignment; }
    vector<int> getSegLengths() const { return segLengths; }
    int segLength(int i) const { return segLengths[i]; }
    resAddress segCentralResidue(int i) const { return resAddress(targetIndex, alignment[i] + segLengths[i]/2); }
    int numSegments() const { return alignment.size(); }
    int operator[](int i) const { return alignment[i]; }
    Transform getTransform() const { return tr; }
    bool seqContextDefined() const { return ((context != NULL) && (!context->segSeq.empty())); }
    void addSequenceContext(const Sequence& targSeq, int contLen, const vector<int>& targChainBeg, const vector<int>& targChainEnd);
    void addSequenceContext(const Structure& target, int contLen);
    void setSeqContext(const vector<Sequence>& _segSeq, const vector<Sequence>& _nSeq, const vector<Sequence>& _cSeq);
    void setStructContext(const vector<AtomPointerVector>& _segStr, const vector<AtomPointerVector>& _nStr, const vector<AtomPointerVector>& _cStr);
    void setRMSD(mstreal r) { rmsd = r; }
    void setTransform(const Transform& _tr) { tr = _tr; }
    fasstSolutionAddress getAddress() const { return fasstSolutionAddress(targetIndex, alignment); }

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

    friend bool operator==(const fasstSolution& si, const fasstSolution& sj) {
      return !((si < sj) || (sj < si));
    }

    void write(ostream& _os) const;
    void read(istream& _is);

    friend ostream& operator<<(ostream &_os, const fasstSolution& _sol) {
      _os << std::setprecision(6) << std::fixed << _sol.rmsd << " " << _sol.targetIndex << " [" << MstUtils::vecToString(_sol.alignment, ", ") << "]";
      if (_sol.seqContextDefined()) {
        for (int i = 0; i < _sol.numSegments(); i++) {
          _os << "\n\t" << "[" << _sol.nTermContext()[i].toString() << "] - [" << _sol.segmentSeqs()[i].toString() << "] - [" << _sol.cTermContext()[i].toString() << "]";
        }
      }
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
    typedef fasstSolution::resAddress resAddress;

    fasstSolutionSet() { updated = false; }
    fasstSolutionSet(const fasstSolution& sol);
    fasstSolutionSet(const fasstSolutionSet& sols);
    fasstSolutionSet(const vector<fasstSolutionAddress>& addresses, const vector<int>& segLengths);
    fasstSolutionSet& operator=(const fasstSolutionSet& sols);
    bool insert(const fasstSolution& sol, mstreal redundancyCut = 1); // returns whether the insert was performed
    bool insert(const fasstSolution& sol, simpleMap<resAddress, tightvector<resAddress>>& relMap); // returns whether the insert was performed
    void erase(fasstSolution& sol);
    set<fasstSolution>::iterator erase(const set<fasstSolution>::iterator it);
    set<fasstSolution>::iterator begin() const { return solsSet.begin(); }
    set<fasstSolution>::reverse_iterator rbegin() const { return solsSet.rbegin(); }
    set<fasstSolution>::iterator end() const { return solsSet.end(); }
    set<fasstSolution>::reverse_iterator rend() const { return solsSet.rend(); }
    // const fasstSolution& operator[] (int i) const;
    fasstSolution& operator[] (int i);
    int size() const { return solsSet.size(); }
    void init(int numSegs) { clear(); solsByCenRes.resize(numSegs); algnRedBar.resize(numSegs); algnRedBarSource.resize(numSegs); }
    void clear() { solsSet.clear(); solsByCenRes.clear(); updated = true; }
    mstreal worstRMSD() { return (solsSet.rbegin())->getRMSD(); }
    mstreal bestRMSD() { return (solsSet.begin())->getRMSD(); }
    vector<fasstSolution*> orderByDiscovery();
    vector<fasstSolutionAddress> extractAddresses() const;

    void clearTempData() { solsByCenRes.clear(); }
    void resetAlignRedBarrierData(int targetLen) {
      for (int i = 0; i < algnRedBar.size(); i++) {
        algnRedBar[i].clear(); algnRedBar[i].resize(targetLen, INFINITY);
        algnRedBarSource[i].clear();
      }
    }
    bool isAlignRedBarrierDataSet() const { return !algnRedBar.empty() && !algnRedBar[0].empty(); }
    mstreal alignRedBarrier(int segIdx, int ri) const { return algnRedBar[segIdx][ri]; }

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
    // solsByCenRes[i] is a map that of solutions keyed by the central residue
    // of the i-th segment. This is only used during search (for redundancy re-
    // moval) and is deleted before return. Not copied upon assignment.
    vector<map<fasstSolution::resAddress, set<fasstSolution*>>> solsByCenRes;

    // these two member variables are for storing information about existing
    // matches that have redundant segments to potential solutions from the
    // current target
    vector<vector<mstreal>> algnRedBar;
    vector<map<fasstSolution*, set<int>>> algnRedBarSource;

    bool updated;
};

/* A general virtual class for representing per-segment sequence constraints. */
class fasstSeqConst {
  public:
    /* vector alignments is assumed to be pre-allocated. Upon exit, alignments[i]
     * will be true if the alignment of segment segIdx with the target starting
     * at with position i meets all sequence constraints, and false otherwise. */
    virtual void evalConstraint(int segIdx, const Sequence& target, vector<bool>& alignments) = 0;
    virtual bool isSegmentConstrained(int segIdx) = 0;
    virtual ~fasstSeqConst() {}
};

/* A simple implementation of per-segment sequence constraints that allows only
 * for positional constraints to a fixed set of amino acids. */
class fasstSeqConstSimple : public fasstSeqConst {
  public:
    fasstSeqConstSimple(int numSegs) { positions.resize(numSegs); aminoAcids.resize(numSegs); }
    void evalConstraint(int segIdx, const Sequence& target, vector<bool>& alignments);
    bool isSegmentConstrained(int segIdx) { return !positions[segIdx].empty(); }
    bool hasConstraints() const {
      for (int i = 0; i < positions.size(); i++) { if (!positions[i].empty()) return true; }
      return false;
    }

    void addConstraint(int segIdx, int posIdx, const vector<string>& aas) {
      positions[segIdx].push_back(posIdx);
      aminoAcids[segIdx].push_back(set<res_t>());
      for (int i = 0; i < aas.size(); i++) aminoAcids[segIdx].back().insert(SeqTools::aaToIdx(aas[i]));
    }

  private:
    /* positions is of size equal to the number of query segments, positions[i]
     * is the set of indices within the i-th query segment that have constraints
     * defined for them. So if positions[i][j] is some index within the i-th
     * segment that is constrained, then aminoAcids[i][j] is the corresponding
     * set of amino acids allowed at this position. If positions[i] is empty,
     * then no constraints are present for the i-th segment. */
    vector<vector<int> > positions;
    vector<vector<set<res_t> > > aminoAcids;
};

class fasstSearchOptions {
  public:
    fasstSearchOptions() {
      rmsdCutRequested = 0.0;
      maxNumMatches = minNumMatches = suffNumMatches = -1;
      gapConstSet = false;
      diffChainRestSet = false;
      contextLength = 30;
      redundancyCut = 1.0;
      seqConst = NULL;
      verb = false;
    }
    ~fasstSearchOptions() { if (seqConst != NULL) delete(seqConst); }

    /* -- getters -- */
    int getMinNumMatches() const { return minNumMatches; }
    int getMaxNumMatches() const { return maxNumMatches; }
    int getSufficientNumMatches() const { return suffNumMatches; }
    mstreal getRMSDCutoff() const { return rmsdCutRequested; }
    int getMinGap(int i, int j) const { return minGap[i][j]; }
    int getMaxGap(int i, int j) const { return maxGap[i][j]; }
    int getContextLength() const { return contextLength; }
    mstreal getRedundancyCut() const { return redundancyCut; }
    string getRedundancyProperty() const { return redundancyProp; }
    fasstSeqConst* getSequenceConstraints() const { return seqConst; }

    /* -- setters -- */
    void setMinNumMatches(int _min);
    void setMaxNumMatches(int _max);
    void setSufficientNumMatches(int _suff);
    void setRMSDCutoff(mstreal cut) { rmsdCutRequested = cut; }
    void setMinGap(int i, int j, int gapLim); // target topology: [segment i] [gap of at list gapLim long] [segment j]
    void setMaxGap(int i, int j, int gapLim); // target topology: [segment i] [gap of at most gapLim long] [segment j]
    void setChainsDiff(int i, int j);
    void setContextLength(int len) { contextLength = len; }
    void setVerbose(bool _verb) { verb = _verb; }
    /* Normally, the redundancy cutoff is between 0 and 1. But one can set it to
     * values outside of this range, in principle. Setting it to a value above 1
     * will cause no redundancy cutoff to be applied, but will populate solution
     * objects with sequence context information, in case (for example) a filter
     * for redundancy will need to be applied later. */
    void setRedundancyCut(mstreal _redundancyCut = 0.5) { redundancyCut = _redundancyCut; }
    void setRedundancyProperty(const string& _redProp) { redundancyProp = _redProp; }
    template<class T>
    void setSequenceConstraints(const T& c) { if (seqConst != NULL) delete(seqConst); seqConst = new T(c); }

    /* -- unsetters (resetters) -- */
    void unsetMinNumMatches() { minNumMatches = -1; }
    void unsetMaxNumMatches() { maxNumMatches = -1; }
    void unsetSufficientNumMatches() { suffNumMatches = -1; }
    void resetGapConstraints(int numQuerySegs);
    void resetDiffChainConstraints(int numQuerySegs);
    void unsetRedundancyCut() { redundancyCut = 1; }
    void unsetRedundancyProperty() { redundancyProp = ""; }
    void unsetSequenceConstraints() { if (seqConst != NULL) delete(seqConst); seqConst = NULL; }

    /* -- queriers -- */
    bool isMinNumMatchesSet() const { return (minNumMatches > 0); }
    bool isMaxNumMatchesSet() const { return (maxNumMatches > 0); }
    bool isSufficientNumMatchesSet() const { return (suffNumMatches > 0); }
    bool minGapConstrained(int i, int j) const { return minGapSet[i][j]; }
    bool maxGapConstrained(int i, int j) const { return maxGapSet[i][j]; }
    bool diffChainsConstrained(int i, int j) const { return diffChainSet[i][j]; }
    bool gapConstrained(int i, int j) const { return (minGapSet[i][j] || maxGapSet[i][j]); }
    bool gapConstraintsExist() const { return gapConstSet; }
    bool diffChainsConstsExist() const { return diffChainRestSet; }
    bool isRedundancyCutSet() const { return redundancyCut < 1; }
    bool isRedundancyPropertySet() const { return !redundancyProp.empty(); }
    bool sequenceConstraintsSet() const { return seqConst != NULL; }
    bool isVerbose() const { return verb; }

    /* -- validators -- */
    bool validateGapConstraints(int numQuerySegs) const;
    bool validateSearchRequest(int numQuerySegs) const; // make sure all user specified requirements are consistent
    bool areNumMatchConstraintsConsistent() const;

  private:
    mstreal rmsdCutRequested;
    int contextLength;                       // how long of a local window to consider when comparing segment alingments between solutions
    mstreal redundancyCut;                   // maximum sequence identity (as a fraction), defined over contextLength-long windows
    string redundancyProp;                   // residue relational property to use for checking redundancy of sequence windows

    vector<vector<int> > minGap, maxGap;     // minimum and maximum sequence separations allowed between each pair of segments
    vector<vector<bool> > minGapSet, maxGapSet, diffChainSet;
    bool gapConstSet, diffChainRestSet, verb;
    int maxNumMatches, minNumMatches, suffNumMatches;
    fasstSeqConst* seqConst;
};

/* FASST -- Fast Algorithm for Searching STructure */
class FASST {
  public:
    enum matchType { REGION = 1, FULL, WITHGAPS };
    enum searchType { CA = 1, FULLBB };
    enum targetFileType { PDB = 1, BINDATABASE, STRUCTURE };
    typedef fasstSolution::resAddress resAddress;

    class targetInfo {
      public:
        targetInfo(const string& _file, targetFileType _type, streampos _loc, short _memSave) {
          file = _file; type = _type; loc = _loc; memSave = _memSave;
        }
        targetInfo(const targetInfo& I) {
          file = I.file; type = I.type; loc = I.loc; memSave = I.memSave;
        }
        string file;         // source file
        targetFileType type; // file type (PDB or database)
        streampos loc;       // location info within file
        short memSave;       // what memory save setting was the target read with?
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
        // same as above, but instead of a list of okay options, takes a vector
        // indicating whether any of the previously known options are okay.
        void intersectOptions(const vector<bool>& isOK);
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
    int getNumQuerySegments() const { return queryStruct.chainSize(); }
    AtomPointerVector getQuerySearchedAtoms() const;
    void addTarget(const Structure& T, short memSave = 0);
    void addTarget(const string& pdbFile, short memSave = 0);
    void addTargets(const vector<string>& pdbFiles, short memSave = 0);
    void stripSidechains(Structure& S);

    void addResidueStringProperties(int ti, const string& propType, const vector<string>& propVals);
    void addResidueProperties(int ti, const string& propType, const vector<mstreal>& propVals);
    void addResiduePairProperties(int ti, const string& propType, const map<int, map<int, mstreal> >& propVals);
    // void addResidueRelationships(int ti, const string& propType, const map<int, map<int, map<int, set<int> > >& resRels);
    void addResidueRelationship(int ti, const string& propType, int ri, int tj, int rj);
    map<int, vector<resAddress>> getResidueRelationships(int ti, const string& propType);

    bool isResiduePropertyDefined(const string& propType);
    bool isResiduePropertyDefined(const string& propType, int ti);
    bool hasResidueProperty(int ti, const string& propType, int ri);
    bool isResidueStringPropertyDefined(const string& propType);
    bool isResidueStringPropertyDefined(const string& propType, int ti);
    bool hasResidueStringProperty(int ti, const string& propType, int ri);
    mstreal getResidueProperty(int ti, const string& propType, int ri);
    string getResidueStringProperty(int ti, const string& propType, int ri);
    bool isResiduePairBoolPropertyDefined(int ti, const string& propType);
    bool isResiduePairBoolPropertyDefined(int ti, const string& propType, int ri);
    set<int> getResiduePairBoolProperty(int ti, const string& propType, int ri);
    bool hasResiduePairProperties(int ti, const string& propType, int ri);
    mstreal isResiduePairPropertyPopulated(const string& propType);
    map<int, mstreal> getResiduePairProperties(int ti, const string& propType, int ri);
    mstreal isResidueRelationshipPopulated(const string& propType);
    void dropResidueRelationship(const string& propType) { resRelProperties.erase(propType); }

    // access to search options
    fasstSearchOptions& options() { return opts; }
    void setOptions(const fasstSearchOptions& _opts) { opts = _opts; }
    // NOTE: these getters/setters are depricated! Use FASST::options().[get/set whatever]
    void setMinNumMatches(int min) { opts.setMinNumMatches(min); }
    void setMaxNumMatches(int max) { opts.setMaxNumMatches(max); }
    void setSufficientNumMatches(int suff) { opts.setSufficientNumMatches(suff); }
    void setRMSDCutoff(mstreal cut) { opts.setRMSDCutoff(cut); }
    void setMinGap(int i, int j, int gapLim) { opts.setMinGap(i, j, gapLim); }
    void setMaxGap(int i, int j, int gapLim) { opts.setMaxGap(i, j, gapLim); }
    void setRedundancyCut(mstreal cut = 0.5) { opts.setRedundancyCut(cut); }
    void setRedundancyProperty(const string& prop) { opts.setRedundancyProperty(prop); }
    void setVerbose(bool _verb) { opts.setVerbose(_verb); }
    int getMinNumMatches() const { return opts.getMinNumMatches(); }
    int getMaxNumMatches() const { return opts.getMaxNumMatches(); }
    int getSufficientNumMatches() const { return opts.getSufficientNumMatches(); }
    mstreal getRMSDCutoff() const { return opts.getRMSDCutoff(); }
    int getMinGap(int i, int j) const { return opts.getMinGap(i, j); }
    int getMaxGap(int i, int j) const { return opts.getMaxGap(i, j); }
    int getContextLength() const { return opts.getContextLength(); }
    mstreal getRedundancyCut() const { return opts.getRedundancyCut(); }
    string getRedundancyProperty() const { return opts.getRedundancyProperty(); }
    bool isMinNumMatchesSet() const { return opts.isMinNumMatchesSet(); }
    bool isMaxNumMatchesSet() const { return opts.isMaxNumMatchesSet(); }
    bool isSufficientNumMatchesSet() const { return opts.isSufficientNumMatchesSet(); }
    bool minGapConstrained(int i, int j) const { return opts.minGapConstrained(i, j); }
    bool maxGapConstrained(int i, int j) const { return opts.maxGapConstrained(i, j); }
    bool gapConstrained(int i, int j) const { return opts.gapConstrained(i, j); }
    bool gapConstraintsExist() const { return opts.gapConstraintsExist(); }
    bool isRedundancyCutSet() const { return opts.isRedundancyCutSet(); }
    bool isRedundancyPropertySet() const { return opts.isRedundancyPropertySet(); }
    bool isVerbose() const { return opts.isVerbose(); }

    int numTargets() const { return targetStructs.size(); }
    Structure getTargetCopy(int i) const { return *(targetStructs[i]); }
    Structure* getTarget(int i) { return targetStructs[i]; }
    int getTargetResidueSize(int i) const { return (targetStructs[i] == NULL) ? atomToResIdx(targets[i].size()) : targetStructs[i]->residueSize(); }
    string getTargetName(int ti) const { return (targetStructs[ti] == NULL) ? "not-saved" : targetStructs[ti]->getName(); }
    Sequence getTargetSequence(int i) { return targSeqs[i]; }
    void setSearchType(searchType _searchType);
    void setGridSpacing(mstreal _spacing) { gridSpacing = _spacing; updateGrids = true; }
    fasstSolutionSet search();
    int numMatches() { return solutions.size(); }

    fasstSolutionSet getMatches() { return solutions; }
    string toString(const fasstSolution& sol);
    void writeDatabase(const string& dbFile);
    /* The memSave parameter can be used to reduce the memory footprint. 0 is the
     * default and does not do any memory savings. 1 means strip the side-chains
     * (for when we will not typically need this information). 2 means destroy
     * the original target structure upon reading, and only keep backbone coordi-
     * nates. This is quite useful in practice, but it does mean that no reagions
     * in the original target that are skipped over (e.g., due to missing backbone
     * atoms) can be tollerated.*/
    void readDatabase(const string& dbFile, short memSave = 0);

    // get various match properties
    void getMatchStructure(const fasstSolution& sol, Structure& match, bool detailed = false, matchType type = matchType::REGION, bool algn = true);
    Structure getMatchStructure(const fasstSolution& sol, bool detailed = false, matchType type = matchType::REGION, bool algn = true);
    void getMatchStructures(fasstSolutionSet& sols, vector<Structure>& matches, bool detailed = false, matchType type = matchType::REGION, bool algn = true);
    vector<Sequence> getMatchSequences(fasstSolutionSet& sols, matchType type = matchType::REGION);
    Sequence getMatchSequence(const fasstSolution& sol, matchType type = matchType::REGION);
    vector<vector<mstreal> > getResidueProperties(fasstSolutionSet& sols, const string& propType, matchType type = matchType::REGION);
    vector<mstreal> getResidueProperties(const fasstSolution& sol, const string& propType, matchType type = matchType::REGION);
    vector<int> getMatchResidueIndices(const fasstSolution& sol, matchType type = matchType::REGION); // figure out the range of residues to excise from target structure
    void addSequenceContext(fasstSolutionSet& sol); // decorate all solutions in the set with sequence context

    /* Computes the RMSD of the given match to a query that is (possibly)
     * different from the one it was found with (if it even came from a search). */
    vector<mstreal> matchRMSDs(fasstSolutionSet& sols, const AtomPointerVector& query, bool update = false);

    /* Filters out redundancy from the match set, building the local sequence
     * context, if necessary, and returns the resulting non-redundancy solution set. */
    fasstSolutionSet removeRedundancy(const fasstSolutionSet& matches, mstreal cut);

    /* Filters out redundancy from the match set using the pre-computed redundancy
     * property in the FASST database. */
    fasstSolutionSet removeRedundancy(const fasstSolutionSet& matches);

    simpleMap<resAddress, tightvector<resAddress>>& getRedundancyPropertyMap() { return resRelProperties[opts.getRedundancyProperty()]; }

  protected:
    void processQuery();
    void setCurrentRMSDCutoff(mstreal cut, int p = -1); // set the current RMSD for this priority level
    void resetCurrentRMSDCutoff(int p = -1);            // reset current RMSD cutoff back to the value set for the given priority level
    int rmsdPriority() const { return rPrior; }
    mstreal getCurrentRMSDCutoff() const { return rmsdCut; }
    void prepForSearch(int ti);
    bool parseChain(const Chain& S, AtomPointerVector* searchable = NULL, Sequence* seq = NULL);
    mstreal currentAlignmentResidual(bool compute, bool setTransform = false);   // computes the accumulated residual up to and including segment recLevel
    mstreal boundOnRemainder(bool compute);           // computes the lower bound expected from segments recLevel+1 and on
    Transform currentTransform();                     // tansform for the alignment corresponding to the current residual
    // void updateQueryCentroids();                      // assumes that appropriate transformation matrices were previously set with a call to currentAlignmentResidual(true)
    int resToAtomIdx(int resIdx) const { return resIdx * atomsPerRes; }
    int atomToResIdx(int atomIdx) const { return atomIdx / atomsPerRes; }
    mstreal centToCentTol(int i);
    mstreal segCentToPrevSegCentTol(int i);
    void rebuildProximityGrids();
    void addTargetStructure(Structure* targetStruct, short memSave = 0);
    void addSequenceContext(fasstSolution& sol); // decorate the solution with sequence context
    void fillTargetChainInfo(int ti);

  private:
    fasstSearchOptions opts;

    /* targetStructs[i] and targets[i] store the original i-th target structure
     * and just the part of it that will be searched over, respectively. NOTE:
     * Atoms* in targets[i] point to Atoms of targetStructs[i]. So there is only
     * one copy of the target structure and there is no amboguity in mapping
     * targets[i] to targetStructs[i], even if when building tagrets[i] some of
     * the residues in targetStructs[i] had to be ignored for whatever reason
     * (e.g., missing backbone). NOTE: it is also possible for any targetStructs[i]
     * entry to be NULL. This means that the full structure was not retained, that
     * atoms in targets[i] are "disembodied" copies that have been stripped of
     * any extra information except coordinates, and therefore is a one-to-one
     * correspondence between residues in the original full structure and atom
     * indices in targets[i] (i.e., nothing was skipped upon processing). This
     * is done for memory reasons. */
    vector<Structure*> targetStructs;
    vector<AtomPointerVector> targets;

    vector<Sequence> targSeqs;               // target sequences (of just the parts that will be searched over)
    vector<targetInfo> targetSource;         // from where and how each target was read (e.g., in case need to re-read it)
    vector<tightvector<int>> targetChainLen; // chain lengths in each target, listed in the order chains appear in the corresponding Structure
    vector<int> targChainBeg, targChainEnd;  // targChainBeg[i] and targChainEnd[i] contain the chain start and end indices for the chain
                                             // that contains the residue with index i (in the overal concatenated sequence). Residue indices
                                             // are based on just the portion of the structure to be searched over.

    /* Object for holding real-valued residue properties. Specifically,
     * resProperties["env"][ti][ri] is the value of the "env" property for
     * residue ri in target with index ti. */
    map<string, map<int, vector<mstreal> > > resProperties;

    /* Object for holding string residue properties. Specifically,
     * resProperties["stride"][ti][ri] is the string value of the "stride" property for
     * residue ri in target with index ti. */
    map<string, map<int, vector<string> > > resStringProperties;

    /* Object for holding binary residue-pair properties. Specifically,
     * resPairProperties["int"][ti][ri] is the set of all rj that interact
     * with ri in a target with index ti.
     * NOTE: this property can be directional (i.e., pairs are not mirrored). */
    map<string, map<int, map<int, set<int> > > > resPairBoolProperties;

    /* Object for holding real-valued residue-pair properties. Specifically,
     * resPairProperties["cont"][ti][ri][rj] is the value of the "cont" property
     * (e.g., contact degree) between residues ri and rj in target with index ti.
     * NOTE: this property can be directional (i.e., pairs are not mirrored). */
    map<string, map<int, map<int, map<int, mstreal> > > > resPairProperties;

    /* Object for holding residue-pair relational graphs. Specifically,
     * resRelProperties["sim"][(ti, ri)] is the list of all residues in the data-
     * base that are related by the property "sim" to the residue with the address
     * (ti, ri) (i.e., target ti, residue ri).
     * NOTE: this property can be directional (i.e., relationships are not mirrored). */
    map<string, simpleMap<resAddress, tightvector<resAddress>>> resRelProperties;

    vector<Transform> tr;                    // transformations from the original frame to the common frames of reference for each target
    int currentTarget;                       // the index of the target currently being searched for

    Structure queryStruct;
    vector<AtomPointerVector> queryOrig;     // just the part of the query that will be sought, split by segment
    vector<AtomPointerVector> query;         // same as above, but with segments re-orderd for optimal searching
    vector<int> qSegOrd;                     // qSegOrd[i] is the index (in the original queryOrig) of the i-th segment in query
    mstreal xlo, ylo, zlo, xhi, yhi, zhi;    // bounding box of the search database

    // segmentResiduals[i][j] is the residual of the alignment of segment i, in which
    // its starting residue aligns with the residue index j in the target
    vector<vector<mstreal> > segmentResiduals;

    int recLevel;                            // segments up to index recLevel are already placed
    int atomsPerRes;
    searchType type;
    vector<vector<string> > searchableAtomTypes;
    int querySize;

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

    // current RMSD and residual cutoffs (not user-set, but internal)
    mstreal rmsdCut, residualCut;

    // facilitate the storage of a series of temporary RMSD cutoffs, by priority
    vector<mstreal> rmsdCutTemp; // rmsdCutTemp[i] is the RMSD cutoff at priority level i (value is negative if not set)
    mstreal rmsdCutDef;          // RMSD cutoff for the top priority (priority value -1)
    int rPrior;                  // the priority level of the current RMSD (-1 if not currently at a temporary RMSD)

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
