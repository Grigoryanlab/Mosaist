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
     * residue(s) belong to. NOTE: the following functions are meant to extract
     * a sensible structure. So, if segments overlap, residues will only be listed
     * once (and the overlapping segments will be part of a single chain). For
     * this reason, the order in which central residues are specified does not
     * influence the order of residues in the selected TERM, and the latter is
     * dictated by the topology of the structure from which the TERM is excised.
     * But this may not always be what one wants, so look at the function below. */
    /* case 1: single central residue */
    static void selectTERM(Residue& cenRes, ConFind& C, Structure& frag, int pm = 2, mstreal cdCut = 0.01, vector<int>* fragResIdx = NULL);
    /* case 2: one or more central residues (given as a vector of Residue*) */
    static void selectTERM(const vector<Residue*>& cenRes, ConFind& C, Structure& frag, int pm = 2, mstreal cdCut = 0.01, vector<int>* fragResIdx = NULL);
    /* case 3: same as case 1, but returns a Structure rather than accepting a reference */
    static Structure selectTERM(Residue& cenRes, ConFind& C, int pm = 2, mstreal cdCut = 0.01, vector<int>* fragResIdx = NULL);
    /* case 4: same as case 2, but returns a Structure */
    static Structure selectTERM(const vector<Residue*>& cenRes, ConFind& C, int pm = 2, mstreal cdCut = 0.01, vector<int>* fragResIdx = NULL);
    /* case 5: all central residues + contacting residues (together referred to as source residues) are already specified */
    static void selectTERM(const vector<Residue*>& cenRes, Structure& frag, int pm = 2, vector<int>* fragResIdx = NULL);

    // /* The following function is a little different, in that it simply cuts out
    //  * one segment at a time, and splices them together. It does not worry about
    //  * any overlap between segments and the order of segments in the resulting
    //  * TERM (i.e., the order of atoms) is according to the order in which the
    //  * central residues are listed in the input. */
    // static void exciseTERM(const vector<Residue*>& cenRes, vector<Atom*>& frag, int pm = 2);

    /* The following function is a little different, in that it simply cuts out
     * one segment at a time, in the order corresponding to the given array out
     * source residues, and places them in separate chains. This is done even
     * if residues from different segments overlap! So that means one can end
     * up with copies of the same residue(s) in different chains. The return
     * value is true if no such overlap occurs (so the resultant TERM is a fine
     * structure) and false if there is overlap. */
    static bool exciseTERM(const vector<Residue*>& cenRes, Structure& frag, int pm = 2);
};

#endif
