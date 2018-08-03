#ifndef _MSTFASSTCACHE_H
#define _MSTFASSTCACHE_H

#include "msttypes.h"
#include "mstfasst.h"

class fasstCache {
  public:
    fasstCache(FASST* s, int max = 1000) {
      S = s; maxNumResults = max;
      sf = 10.0;
      maxNumPressure = (exp(1.0) - 1.0)*maxNumResults/sf; // so that max factor is 2.0 at the start
      errTolPressure = (exp(0.1) - 1.0)*maxNumResults/sf; // so that error factor is 1.1 at the start
    }
    ~fasstCache() { clear(); }
    void clear(); // clears cache (removes all solutions)
    FASST* getFASST() const { return S; }
    fasstSolutionSet search(bool verb = false);
    int getMaxNumResults() const { return maxNumResults; }
    void incErrTolPressure(mstreal del = 1.0) { errTolPressure += fabs(del); }
    void incMaxNumPressure(mstreal del = 1.0) { maxNumPressure += fabs(del); }
    void decErrTolPressure(mstreal del = 1.0) { errTolPressure -= fabs(del); if (errTolPressure < 0) errTolPressure = 0; }
    void decMaxNumPressure(mstreal del = 1.0) { maxNumPressure -= fabs(del); if (maxNumPressure < 0) maxNumPressure = 0; }
    mstreal getErrTolPressure() const { return errTolPressure; }
    mstreal getMaxNumPressure() const { return maxNumPressure; }
    mstreal maxNumFactor() const { return log(sf*maxNumPressure/maxNumResults + 1.0) + 1.5; }
    mstreal errTolFactor() const { return log(sf*errTolPressure/maxNumResults + 1.0) + 1.05; }

    // write/read cache object to/from a binary file
    void write(const string& filename) const;
    void write(ostream &_os) const;
    void read(const string& filename);
    void read(istream &_is);

  protected:
    static vector<int> getStructureTopology(const Structure& S);

  private:
    class fasstCachedResult {
      friend class fasstCache;
      public:
        fasstCachedResult() { priority = 1.0; searchRMSDcut = rmsdCut = 0.0; searchMaxNumMatches = 0; }
        fasstCachedResult(const AtomPointerVector& q, const fasstSolutionSet& sols, mstreal cut, int max, vector<int> topo);
        fasstCachedResult(const fasstCachedResult& r);
        ~fasstCachedResult();
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

        friend bool operator<(const fasstCachedResult& ri, const fasstCachedResult& rj) {
          if (ri.priority != rj.priority) return (ri.priority > rj.priority);
          if (ri.solSet.size() != rj.solSet.size()) return ri.solSet.size() < rj.solSet.size();
          if (ri.query.size() != rj.query.size()) return ri.query.size() < rj.query.size();
          return &ri < &rj;
        }

        friend ostream& operator<<(ostream &_os, const fasstCachedResult& _res) {
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
      bool operator() (const fasstCachedResult* lhs, const fasstCachedResult* rhs) const {
        return (*lhs < *rhs);
      }
    };

    int maxNumResults; // max number of searches to cache
    FASST* S;
    set<fasstCachedResult*, compResults> cache;
    RMSDCalculator rc;
    mstreal errTolPressure, maxNumPressure, sf;
};

#endif
