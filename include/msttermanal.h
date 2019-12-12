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
    TERMANAL(FASST* _F = NULL) { F = _F; cdCut = 0.1; pad = 2; pseudoCount = 0.01; matchCount = 50; rmsdCut = 2.0; }
    void setFASST(FASST* _F) { F = _F; }

    /* Computes the structure score for the given TERM, following roughly the
     * definition and procedure defined in Zheng, Zhang, and Grigoryan, Structure,
     * 23(5), 2015 (DOI 10.1016/j.str.2015.03.015). This is a slight generalization
     * of the scoring described in the paper, in that the design score component
     * can account for the amino-acid in not just one, but any number of positions
     * (by adding the contribution of each). The set of positions at which to
     * account for amino acids is given as the second optional parameter. */
    mstreal structureScore(const Structure& S, const vector<Residue*>& central = vector<Residue*>(), bool verbose = false);
    pair<mstreal, mstreal> structureScoreParts(const Structure& S, const vector<Residue*>& central = vector<Residue*>(), bool verbose = false);
    mstreal combineScoreParts(const pair<mstreal, mstreal>& parts);
    vector<mstreal> scoreStructure(const Structure& S, const vector<Residue*>& subregion, bool verbose = false);
    vector<mstreal> scoreStructure(const Structure& S, bool verbose = false) { return scoreStructure(S, S.getResidues(), verbose); }
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

  protected:
    fasstSearchOptions setupSearch(const Structure& S);
    mstreal calcSeqLikelihood(const vector<Residue*>& central, vector<Sequence>& matchSeqs, bool verbose);
    mstreal calcSeqFreq(Residue* res, vector<Sequence>& matchSeqs);
    mstreal calcStructFreq(fasstSolutionSet& matches, bool verbose);
    vector<Structure> collectTERMs(ConFind& C, const vector<Residue*>& subregion, vector<vector<Residue*>>& termCentrals, vector<vector<int>>& resOverlaps);
    mstreal smoothScore(vector<pair<mstreal, mstreal>>& structScoreParts, int resInd, vector<int>& overlapResInds);

  private:
    FASST* F;
    RotamerLibrary RL;
    mstreal cdCut; // CD cutoff for creating TERMs
    int pad; // TERM padding
    mstreal pseudoCount; // pseudocount for the structure score
    int matchCount; // minimum and maximum number of matches to search for
    mstreal rmsdCut; // RMSD cutoff when searching for matches
};


#endif
