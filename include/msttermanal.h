#ifndef _MSTTERMANAL_H
#define _MSTTERMANAL_H

#include "msttypes.h"
#include "mstfasst.h"
#include "mstsequence.h"
#include "mstmagic.h"
using namespace MST;

class TERMANAL {
  public:
    TERMANAL(FASST* _F = NULL) { F = _F; }
    void setFASST(FASST* _F) { F = _F; }

    /* Computes the structure score for the given TERM, following roughly the
     * definition and procedure defined in Zheng, Zhang, and Grigoryan, Structure,
     * 23(5), 2015 (DOI 10.1016/j.str.2015.03.015), except that no smoothing is
     * performed over neighboring positions (since no larger structural context
     * is given). This is a slight generalization of the scoring described in the
     * paper, in that the design score component can account for the amino-acid
     * in not just one, but any number of positions (by adding the contribution
     * of each). The set of positions at which to account for amino acids is
     * given as the second optional parameter. */
    mstreal structureScore(const Structure& S, const vector<Residue*>& central = vector<Residue*>(), bool verbose = false);

    fasstSolutionSet findTopN(const Structure& S, int N);

  private:
    FASST* F;

};


#endif
