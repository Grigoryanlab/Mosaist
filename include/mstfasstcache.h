#ifndef _MSTFASSTCACHE_H
#define _MSTFASSTCACHE_H

#include "msttypes.h"
#include "mstfasst.h"

// TODO:
// 1. write a test function for the cache.
// 2. write a "build" function for cache to build up a "good" cache.
// 3. consider all permutations of segment elements when looking for a cached search:
//    * would it be better to keep all queries as a separate place FASST database?
//    * if so, we will want to search against it and recover ALL matches (or set maxN to say 10000)
//    * then, if the new query fits within the old query, but has room left, then the RMSD radius and so on will need to be updated
//      (and usually, it would be highly disadvantageous to have room left, because we can assume nothing about those atoms)
// 4. test the impact of doing #3 above
// 5. re-introduce minN capability back inot cFASST

class cFASST : public FASST {
  public:
    cFASST(int max = 1000) : FASST() {
      maxNumResults = max;
      sf = 10.0;
      maxNumPressure = (exp(1.0) - 1.0)*maxNumResults/sf; // so that max factor is 2.0 at the start
      errTolPressure = (exp(0.1) - 1.0)*maxNumResults/sf; // so that error factor is 1.1 at the start
      readPerm = modPerm = true;
      /* Because removal of redundancy is done on-the-fly in FASST (as matches
       * are discovered) and after the fact in cFASST (after fully redundant
       * matches are first fetched from the cache), the precise set of matches
       * may be slightly different in the two cases when redundancy filter is
       * turned on. We can, of course, make sure that cFASST returns the same
       * exact results as FASST when redundancy is on, but this takes a bit of
       * extra computation (not too much, but still). Considering that it is
       * conceptually cleaner to do after-the-fact filtering (as it is indepedent
       * of the order in which matches are discovered and thus indepedent of the
       * order in which the database is listed), it does not seem worth ensuring
       * this exact equivalence in most cases. But in some cases, like when we
       * are testing the implementation of the cache for correctness, this may
       * be desired. Setting this flag to true will have this effect. */
      strictEquiv = false;
    }
    ~cFASST() { clear(); }
    void clear(); // clears cache (removes all solutions)
    fasstSolutionSet search();
    int getMaxNumResults() const { return maxNumResults; }
    void incErrTolPressure(mstreal del = 1.0) { errTolPressure += fabs(del); }
    void incMaxNumPressure(mstreal del = 1.0) { maxNumPressure += fabs(del); }
    void decErrTolPressure(mstreal del = 1.0) { errTolPressure -= fabs(del); if (errTolPressure < 0) errTolPressure = 0; }
    void decMaxNumPressure(mstreal del = 1.0) { maxNumPressure -= fabs(del); if (maxNumPressure < 0) maxNumPressure = 0; }
    void setReadPermission(bool _perm) { readPerm = _perm; }
    void setModifyPermission(bool _perm) { modPerm = _perm; }
    void setStrictEquivalence(bool _eq) { strictEquiv = _eq; }
    mstreal getErrTolPressure() const { return errTolPressure; }
    mstreal getMaxNumPressure() const { return maxNumPressure; }
    mstreal maxNumFactor() const { return log(sf*maxNumPressure/maxNumResults + 1.0) + 1.5; }
    mstreal errTolFactor() const { return log(sf*errTolPressure/maxNumResults + 1.0) + 1.05; }
    bool readPermission() const { return readPerm; }
    bool modifyPermission() const { return modPerm; }
    bool strictEquivalence() const { return strictEquiv; }

    // write/read cache object to/from a binary file
    void write(const string& filename) const;
    void write(ostream &_os) const;
    void read(const string& filename);
    void read(istream &_is);

  protected:
    static vector<int> getStructureTopology(const Structure& S);

  private:
    class cachedResult {
      friend class cFASST;
      public:
        cachedResult() { priority = 1.0; searchRMSDcut = rmsdCut = 0.0; searchMaxNumMatches = 0; }
        cachedResult(const AtomPointerVector& q, const fasstSolutionSet& sols, mstreal cut, int max, vector<int> topo);
        cachedResult(const cachedResult& r);
        ~cachedResult();
        void upPriority(mstreal del = 1.0) { priority += del; }
        void elapsePriority(int Thalf) { priority *= pow(0.5, 1.0/Thalf); }

        AtomPointerVector getQuery() const { return query; }
        vector<fasstSolutionAddress>& getSolutions() { return solSet; }
        mstreal getRMSDCutoff() const { return rmsdCut; }
        mstreal getSearchRMSDCutoff() const { return searchRMSDcut; }
        int getSearchMaxNumMatches() const { return searchMaxNumMatches; }
        mstreal getPriority() const { return priority; }
        vector<int> getTopology() const { return topology; }
        bool isSameTopology(const vector<int>& compTopo) const;
        bool isLimitedByMaxNumMatches() const { return solSet.size() == searchMaxNumMatches; }

        friend bool operator<(const cachedResult& ri, const cachedResult& rj) {
          if (ri.priority != rj.priority) return (ri.priority > rj.priority);
          if (ri.solSet.size() != rj.solSet.size()) return ri.solSet.size() < rj.solSet.size();
          if (ri.query.size() != rj.query.size()) return ri.query.size() < rj.query.size();
          return &ri < &rj;
        }

        friend ostream& operator<<(ostream &_os, const cachedResult& _res) {
          _os << "[" << MstUtils::vecToString(_res.topology) << "] " << _res.searchRMSDcut << "/" << _res.searchMaxNumMatches << " (" << _res.priority << ")" << endl;
          _os << _res.query << endl << MstUtils::vecToString(_res.solSet) << endl;
          return _os;
        }

        void write(ostream &_os) const;
        void read(istream &_is);

      private:
        AtomPointerVector query;
        vector<int> topology;
        vector<fasstSolutionAddress> solSet;
        mstreal priority;
        // search parameters
        mstreal searchRMSDcut, rmsdCut;
        int searchMaxNumMatches;
    };

    struct compResults {
      bool operator() (const cachedResult* lhs, const cachedResult* rhs) const {
        return (*lhs < *rhs);
      }
    };

    int maxNumResults; // max number of searches to cache
    set<cachedResult*, compResults> cache;
    RMSDCalculator rc;
    mstreal errTolPressure, maxNumPressure, sf;
    bool readPerm, modPerm, strictEquiv;
};

typedef cFASST fasstCache; // cFASST seems to be a better name, but the initial one was cFASST

#endif
