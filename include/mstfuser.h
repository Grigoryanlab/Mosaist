#ifndef _MSTFUSER_H
#define _MSTFUSER_H

#include "msttypes.h"
#include "mstoptim.h"

using namespace std;
using namespace MST;

class fusionEvaluator: public optimizerEvaluator {
  public:
    /* vector<vector<Residue*> >& resTopo simultaneously represents the various
     * aligned segments as well as which residues overlap in the alignment. It is
     * assumed that all residues from a give chain move together (i.e., are
     * part of the same aligned fragment). resTopo[i] stores the residues that
     * overlap with residue index i of the eventual fused structure. Thus, the
     * length of resTopo is also the length of the fused structure and resTopo[i]
     * has to have at least one entry for all i. */
    fusionEvaluator(const vector<vector<Residue*> >& resTopo, vector<int> fixedResidues = vector<int>(), bool _verbose = false);

    double eval(const vector<double>& point);

    vector<double> guessPoint() {
      eval(vector<double>());
      noise = 1;
      anchorRes = (fixedResidues.size() > 0) ? fixedResidues[MstUtils::randInt(0, fixedResidues.size() - 1)] : -1;
cout << "anchoring on " << anchorRes << endl;
      return initPoint;
    }
    bool isAnchored() { return anchorRes >= 0; }
    int numDF() {
      int df = 3*numMobileAtoms - 6;
      if (isAnchored()) df += 6;
      return df;
    }
    Structure getStructure() { return fused; }
    void setVerbose(bool _verbose) { verbose = _verbose; }

  protected:
    CartesianPoint bondInstances(int ri, int rj, const string& ai, const string& aj, bool addToCache = false);
    CartesianPoint angleInstances(int ri, int rj, int rk, const string& ai, const string& aj, const string& ak, bool addToCache = false);
    CartesianPoint dihedralInstances(int ri, int rj, int rk, int rl, const string& ai, const string& aj, const string& ak, const string& al, bool addToCache = false);

    // if the last flag is specified as true, will compute angular differences
    double harmonicPenalty(double val, const pair<double, double>& minmax, double K, bool angular = false);

  private:
    Structure fused;
    vector<bool> fixed; // marks whether each residue is to be fixed or not
    vector<int> fixedResidues; // just a list of fixed residue indices. this is
                               // redundant with the above, but helpful to have
    int anchorRes;      // index one of the fixed residues, which will be used
                        // to start placing all other atoms. If there are no
                        // fixed residues, this index is set to -1;
    int numMobileAtoms; // number of non-fixed atoms

    /* Allowed ranges for internal degrees of freedom. Note that the number of
     * these bounds is equal to the number of internal coordinates used in
     * scoring a structure, which is not necessarily the number of degrees of
     * freedom needed to build the structure (i.e., the dimensionality of the
     * optimization problem). This is because some residues may be marked as
     * fixed, such that no degrees of freedom are needed to build them, and yet
     * some degrees of freedom involving those atoms would still neeed to be
     * compared to allowed ranges to score (e.g., the bond between a fixed and
     * a non-fixed atom). */
    vector<pair<double, double> > bounds;

    /* For each residue index in the fused structure, this variable stores the
     * list of all overlapping fragment residues. */
    vector<vector<Residue*> > overlappingResidues;

    /* A list of aligned bits between the fused structure and the overlapping
     * fragments. alignedFrags[i] is a pair of AtomPointerVectors, representing
     * atoms on the fused structure and corresponding atoms on some relevant
     * fragment, respectively. */
    vector<pair<AtomPointerVector, AtomPointerVector> > alignedFrags;

    double kb, ka, kh; // force constants for enforcing bonds, angles, and dihedrals
    vector<double> initPoint;
    bool verbose;
    real noise;
};

class Fuser {
  public:
    static Structure fuse(const vector<vector<Residue*> >& resTopo, const vector<int>& fixed, int Nc = 2, bool verbose = false) {
      fusionEvaluator E(resTopo, fixed, verbose);
      vector<double> bestSolution;
      double bestScore = Optim::fminsearch(E, 1000, bestSolution);
      for (int i = 0; i < Nc-1; i++) {
        vector<double> solution;
        double score = Optim::fminsearch(E, 1000, solution);
        if (score < bestScore) { bestScore = score; bestSolution = solution; }
      }
      cout << "best score = " << bestScore << endl;
      E.eval(bestSolution);
      return E.getStructure();
    }
};


#endif
