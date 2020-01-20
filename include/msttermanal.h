#ifndef _MSTTERMANAL_H
#define _MSTTERMANAL_H

#include "msttypes.h"
#include "mstrotlib.h"
#include "mstcondeg.h"
#include "mstfasst.h"
#include "mstsequence.h"

using namespace MST;

class TERMANAL {
  public:
    TERMANAL(FASST* _F = NULL) { F = _F; cdCut = 0.1; pad = 2; pseudoCount = 0.01; matchCount = 50; rmsdCut = 2.0; compatMode = false; compatSearchLimit = 1000; }
    void setFASST(FASST* _F) { F = _F; }

    /* Computes the structure score for the given TERM, following roughly the
     * definition and procedure defined in Zheng, Zhang, and Grigoryan, Structure,
     * 23(5), 2015 (DOI 10.1016/j.str.2015.03.015). This is a slight generalization
     * of the scoring described in the paper, in that the design score component
     * can account for the amino-acid in not just one, but any number of positions
     * (by adding the contribution of each). The set of positions at which to
     * account for amino acids is given as the second optional parameter. */
    mstreal structureScore(const Structure& term, const vector<Residue*>& central = vector<Residue*>(), bool verbose = false);
    pair<mstreal, mstreal> structureScoreParts(const Structure& term, const vector<Residue*>& central = vector<Residue*>(), bool verbose = false);
    mstreal combineScoreParts(const pair<mstreal, mstreal>& parts);
    // Scores an entire structure or a subregion thereof, optionally storing the design and abundance scores in the provided scoreParts vector
    vector<mstreal> scoreStructure(const Structure& S, const vector<Residue*>& subregion, vector<pair<mstreal, mstreal>>* scoreParts = NULL, bool verbose = false);
    vector<mstreal> scoreStructure(const Structure& S, vector<pair<mstreal, mstreal>>* scoreParts = NULL, bool verbose = false) {
      return scoreStructure(S, S.getResidues(), scoreParts, verbose);
    }
    FASST* getFasstDB() { return F; }
    void setFasstDB(FASST* _F) { F = _F; }
    RotamerLibrary& getRotamerLibrary() { return RL; }
    void readRotamerLibrary(const string& rotLibFile) { RL.readRotamerLibrary(rotLibFile); }
    mstreal getCDCut() { return cdCut; }
    void setCDCut(mstreal _cdCut) { cdCut = _cdCut; }
    int getPad() { return pad; }
    void setPad(int _pad) { pad = _pad; }
    mstreal getPseudocount() { return pseudoCount; }
    void setPseudocount(mstreal _pseudoCount) { pseudoCount = _pseudoCount; }
    int getMatchCount() { return matchCount; }
    void setMatchCount(int _matchCount) { matchCount = _matchCount; }
    mstreal getRMSDCut() { return rmsdCut; }
    void setRMSDCut(mstreal _rmsdCut) { rmsdCut = _rmsdCut; }
    bool getCompatMode() { return compatMode; }
    void setCompatMode(bool _compatMode) { compatMode = _compatMode; }
    int getCompatSearchLimit() { return compatSearchLimit; }
    void setCompatSearchLimit(int _compatSearchLimit) { compatSearchLimit = _compatSearchLimit; }

  protected:
    fasstSearchOptions setupSearch(const Structure& S);
    vector<fasstSolution*> getTopMatches(FASST* F, fasstSolutionSet& matches, vector<Atom*>& queryA, vector<mstreal>* topRmsds = NULL);
    mstreal calcSeqLikelihood(vector<fasstSolution*>& matches, const vector<Residue*>& central, bool verbose);
    mstreal calcSeqFreq(Residue* res, vector<Sequence>& matchSeqs);
    mstreal calcStructFreq(vector<fasstSolution*>& matches, vector<mstreal>& rmsds, bool verbose);
    vector<Structure> collectTERMs(ConFind& C, const vector<Residue*>& subregion, vector<Residue*>& centrals, vector<vector<int>>& resOverlaps);
    pair<mstreal, mstreal> smoothScores(vector<pair<mstreal, mstreal>>& structScoreParts, vector<int>& overlapResInds);

  private:
    FASST* F;
    RotamerLibrary RL;
    mstreal cdCut; // CD cutoff for creating TERMs
    int pad; // TERM padding
    mstreal pseudoCount; // pseudocount for the structure score
    int matchCount; // maximum number of matches to consider per TERM
    mstreal rmsdCut; // RMSD cutoff when searching for matches
    // For compatibility with the original TERMANAL method, set the compatibility mode to true
    // Note that whether a FASST object searches by CA or full backbone RMSD must be determined before the object is loaded, so it is up to the user to setup F correctly
    // Note also that the original method uses a different sequence redundancy filter, which is not implemented the way it is here (FASST::setRedundancyCut)
    bool compatMode; // if true, assumes the FASST search is by CA RMSD and then sorts by full backbone RMSD
    int compatSearchLimit; // maximum number of matches to return when searching in compatibility mode
};


#endif
