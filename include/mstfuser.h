#ifndef _MSTFUSER_H
#define _MSTFUSER_H

#include "msttypes.h"
#include "mstoptim.h"
#include "msttransforms.h"
#include <chrono>

using namespace std;
using namespace MST;

class fusionEvaluator: public optimizerEvaluator {
  public:
    // broken means across a chain break
    enum icType { icBond = 1, icAngle, icDihedral, icBrokenBond, icBrokenAngle, icBrokenDihedral };

    /* vector<vector<Residue*> >& resTopo simultaneously represents the various
     * aligned segments as well as which residues overlap in the alignment. It is
     * assumed that all residues from a give chain move together (i.e., are
     * part of the same aligned fragment). resTopo[i] stores the residues that
     * overlap with residue index i of the eventual fused structure. Thus, the
     * length of resTopo is also the length of the fused structure and resTopo[i]
     * has to have at least one entry for all i. */
    fusionEvaluator(const vector<vector<Residue*> >& resTopo, vector<int> fixedResidues = vector<int>(), bool _verbose = false);

    double eval(const vector<double>& point);

    vector<double> guessPoint() { if (initPoint.empty()) eval(vector<double>()); return initPoint; }
    void setGuessPoint(const vector<double>& _initPoint) { initPoint = _initPoint; }
    void noisifyGuessPoint(mstreal _noise = 1.0) { noise = _noise; initPoint.resize(0); }
    int numResidues() { return overlappingResidues.size(); }
    bool isAnchored() { return fixedResidues.size() > 0; }
    int numDF() {
      int df = 3*numMobileAtoms - 6;
      // TODO: for now, use 3N degrees of freedom when optimizing in XYZ; will change
      // later (need to position the initial conformation so that the first three
      // atoms are in their stardard orientations)
      if (isAnchored()) df += 6;
      return df;
    }
    int getBuildOrigin() { return buildOriginRes; }
    void setBuildOrigin(int _buildOriginRes) { buildOriginRes = _buildOriginRes; }
    Structure getStructure() { return fused; }
    Structure getAlignedStructure();
    void setVerbose(bool _verbose) { verbose = _verbose; }
    int randomizeBuildOrigin() {
      buildOriginRes = (fixedResidues.size() > 0) ? fixedResidues[MstUtils::randInt(0, fixedResidues.size() - 1)] : MstUtils::randInt(0, numResidues() - 1);
      return buildOriginRes;
    }

  class icBound {
    public:
      icBound(icType _type, mstreal _minVal, mstreal _maxVal, string _name = "") { type = _type; minVal = _minVal; maxVal = _maxVal; name = _name; }
      icBound(icType _type, const pair<mstreal, mstreal>& b, string _name = "") { type = _type; minVal = b.first; maxVal = b.second; name = _name; }
      icBound(const icBound& icb) { type = icb.type; minVal = icb.minVal; maxVal = icb.maxVal; name = icb.name; }

      icType type;
      mstreal minVal, maxVal;
      string name;
  };

  protected:
    AtomPointerVector atomInstances(int ri, const string& ai);
    mstreal bondInitValue(int ri, int rj, const string& ai, const string& aj, bool negateStartWithMean = false);
    mstreal angleInitValue(int ri, int rj, int rk, const string& ai, const string& aj, const string& ak);
    mstreal dihedralInitValue(int ri, int rj, int rk, int rl, const string& ai, const string& aj, const string& ak, const string& al);
    CartesianPoint bondInstances(int ri, int rj, const string& ai, const string& aj, bool addToCache = false);
    CartesianPoint angleInstances(int ri, int rj, int rk, const string& ai, const string& aj, const string& ak, bool addToCache = false);
    CartesianPoint dihedralInstances(int ri, int rj, int rk, int rl, const string& ai, const string& aj, const string& ak, const string& al, bool addToCache = false);

    // if the last flag is specified as true, will compute angular differences
    double harmonicPenalty(double val, const icBound& b);

  private:
    Structure fused, guess;
    vector<bool> fixed; // marks whether each residue is to be fixed or not
    vector<int> fixedResidues; // just a list of fixed residue indices. this is
                               // redundant with the above, but helpful to have
    int buildOriginRes;      // index one of the fixed residues, which will be used
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
    vector<icBound> bounds;

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
    bool verbose, startWithMean, optimCartesian;
    mstreal noise;
};

class Fuser {
  public:
    /* Fuses given residues into a single polymer. The parameter resTopo specifies
     * which residues overlap. That is, resTopo.size() is the length of the
     * eventual chain, and resTopo[i] is a list of residues, of length 1 or more
     * that overlap with the i-th residue in the chain, in the N-to-C order,
     * starting with 0. Returns the the fused chain as a Structure object. */
    static Structure fuse(const vector<vector<Residue*> >& resTopo, const vector<int>& fixed = vector<int>(), int Ni = 100, int Nc = 1, bool verbose = false);

    /* This function is a simplified version of Fuser::fuse(), in that it guesses
     * the topology automatically. Argument residues is a flat vector of all the
     * residues from all the different fragments that are to be fused. residues
     * from the different fragments can be intermixed in this vector, BUT when
     * one considers only residues of a given fragment, these must be in N-to-C
     * order. The function guesses the topology by finding residues that are
     * likely overlapping (i.e., have close CA atoms). This should not give
     * incorrect topologies in "normal" circumstances, but in strange cases can. */
    static Structure autofuse(const vector<Residue*>& residues, int flexOnlyNearOverlaps = -1, int Ni = 100, int Nc = 1, bool verbose = false);
};


#endif
