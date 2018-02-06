#ifndef _MSTMAGIC_H
#define _MSTMAGIC_H

#include "msttypes.h"
#include "mstcondeg.h"
#include <algorithm>

using namespace MST;

class TERMUtils {
  public:
    static vector<AtomPointerVector> mostDesignableFragments(Structure& C, vector<Structure*>& TERMs, int n = 10, CartesianPoint* cen = NULL, CartesianPoint* ext = NULL, string outBase = "", bool verb = false);

    /* These functions excise a TERM from a given structure. A TERM is always
     * defined as consisting of some number of central residues, all residues
     * in contact with each (according to the provided ConFind object and the
     * given contact degree cutoff), and +/- pm residues from each of those. But
     * the different functions differ in what they take as input. The last
     * (optional) argument is a pointer to a vector of ints, which is appended
     * with residue indices (from the parent Structure) of the residues that end
     * up comprising the TERM (in the order they are inserted into the TERM).
     * Disjoint fragments are placed into separate chains within the TERM, with
     * chain topolpgy taken from the original structure that the specified
     * residue(s) belong to.*/
    /* case 1: single central residue */
    static void selectTERM(Residue& cenRes, ConFind& C, Structure& frag, int pm = 2, mstreal cdCut = 0.01, vector<int>* fragResIdx = NULL);
    /* case 2: one or more central residues (given as a vector of Residue*) */
    static void selectTERM(const vector<Residue*>& cenRes, ConFind& C, Structure& frag, int pm = 2, mstreal cdCut = 0.01, vector<int>* fragResIdx = NULL);
    /* case 3: same as case 1, but returns a Structure rather than accepting a reference */
    static Structure selectTERM(Residue& cenRes, ConFind& C, int pm = 2, mstreal cdCut = 0.01, vector<int>* fragResIdx = NULL);
    /* case 4: same as case 2, but returns a Structure */
    static Structure selectTERM(const vector<Residue*>& cenRes, ConFind& C, int pm = 2, mstreal cdCut = 0.01, vector<int>* fragResIdx = NULL);
    /* case 5: all central residues + contacting residues are already specified */
    static void selectTERM(const vector<Residue*>& cenRes, Structure& frag, int pm = 2, vector<int>* fragResIdx = NULL);
};

#endif
