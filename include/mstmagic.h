#ifndef _MSTMAGIC_H
#define _MSTMAGIC_H

#include "msttypes.h"
#include "mstcondeg.h"
#include <algorithm>

using namespace MST;

class TERMUtils {
  public:
    static vector<AtomPointerVector> mostDesignableFragments(Structure& C, vector<Structure*>& TERMs, int n = 10, CartesianPoint* cen = NULL, CartesianPoint* ext = NULL, string outBase = "", bool verb = false);

    /* These functions excise a TERM from a given structure. The first one is
     * defines a TERM as being the central residue (first argument), all residues
     * in contact with it (according to the provided ConFind object and the
     * given contact degree cutoff), and +/- pm residues from each of those. The
     * second one instead takes a list of residues that are the central residue
     * and all of the contacting residues. Both functions then excise the TERM,
     * placing it into the specified Structure object, and also optionally record
     * residue indices (from the original structure) of the residues in the final
     * TERM, in the order they are inserted into the TERM. Disjoint fragments
     * are placed into separate chains within the TERM, with chain topolpgy taken
     * from the original structure that the specified residue(s) belong to. The
     * last function does the same thing as the first, but returns the TERM. */
    static void selectTERM(Residue& cenRes, ConFind& C, Structure& frag, int pm = 2, mstreal cdCut = 0.01, vector<int>* fragResIdx = NULL);
    static void selectTERM(const vector<Residue*>& cenRes, Structure& frag, int pm = 2, vector<int>* fragResIdx = NULL);
    static Structure selectTERM(Residue& cenRes, ConFind& C, int pm = 2, mstreal cdCut = 0.01, vector<int>* fragResIdx = NULL);

};

#endif
