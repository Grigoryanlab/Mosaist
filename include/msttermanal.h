#ifndef _MSTTERMANAL_H
#define _MSTTERMANAL_H

#include "msttypes.h"
#include "mstfasst.h"
#include "mstsequence.h"

using namespace MST;

class TERMANAL {
  public:
    TERMANAL(FASST* _F = NULL) { F = _F; }
    void setFASST(FASST* _F) { F = _F; }

    /* Computes the structure score for the given TERM, following roughly the
     * definition and procedure defined in Zheng, Zhang, and Grigoryan, Structure,
     * 23(5), 2015 (DOI 10.1016/j.str.2015.03.015). This is a slight generalization
     * of the scoring described in the paper, in that the design score component
     * can account for the amino-acid in not just one, but any number of positions
     * (by adding the contribution of each). The set of positions at which to
     * account for amino acids is given as the second optional parameter. */
    map<Residue*, vector<Residue*>> getOverlappingTerms(const Structure& S, const vector<Structure*>& otherS);
    map<Residue*, vector<Residue*>> getOverlappingTerms(const vector<Residue*>& R, const vector<Structure*>& otherS);
    mstreal structureScore(const Structure& S, const vector<Residue*>& central = vector<Residue*>(), const map<Residue*, vector<Residue*>>& overlapSets = map<Residue*, vector<Residue*>>(), bool verbose = false);
    fasstSolutionSet findTopN(const Structure& S, int N);

  protected:
    map<Residue*, Residue*> getOverlaps(const vector<Residue*>& R, const Structure& otherS);
    map<pair<char, int>, Residue*> mapResInds(const vector<Residue*>& R);
    mstreal calcSeqFreq(Residue* res, map<Structure*, fasstSolutionSet>& closestMatches);

  private:
    FASST* F;
};


#endif
