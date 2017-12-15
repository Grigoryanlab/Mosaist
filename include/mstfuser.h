#ifndef _MSTFUSER_H
#define _MSTFUSER_H

#include "msttypes.h"
#include "mstoptim.h"
#include "msttransforms.h"
#include <chrono>

using namespace std;
using namespace MST;

class fusionParams {
  public:
    enum coorInitType { meanCoor = 1, meanIC };
    fusionParams() { // default optimization params
      startType = fusionParams::coorInitType::meanCoor;
      verbose = false;
      optimCartesian = true;
      noise = 0;
      Ni = 100;
      Nc = 1;
      tol = 0.0;
    }
    mstreal getNoise() const { return noise; }
    fusionParams::coorInitType getCoorInitType() const { return startType; }
    bool getOptimCartesian() const { return optimCartesian; }
    bool isVerbose() const { return verbose; }
    int numCycles() const { return Nc; }
    int numIters() const { return Ni; }
    mstreal errTol() const { return tol; }

    void setNoise(mstreal _noise) { noise = _noise; }
    void setVerbose(bool _verbose = true) { verbose = _verbose; }
    void setCoorInitType(fusionParams::coorInitType _startType) { startType = _startType; }
    void setOptimCartesian(bool _optimCartesian) { optimCartesian = _optimCartesian; }
    void setNumCycles(int nc) { Nc = nc; }
    void setNumIters(int ni) { Ni = ni; }
    void setErrTol(mstreal _tol) { tol = _tol; }

  private:
    // start optimization from the averaged Cartesian structure or the structure
    // that results from average internal coordinates?
    fusionParams::coorInitType startType;
    bool verbose, optimCartesian;
    mstreal noise; // noise level for initalizing the starting point
    int Ni, Nc;    // number of iterations per cycle and number of cycles
    mstreal tol;   // error tolerance stopping criterion for optimization
};

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
    fusionEvaluator(const vector<vector<Residue*> >& resTopo, vector<int> fixedResidues = vector<int>(), const fusionParams& _params = fusionParams());

    mstreal eval(const vector<mstreal>& point);
    mstreal eval(const vector<mstreal>& point, Vector& grad);
    vector<mstreal> guessPoint();
    void setGuessPoint(const vector<mstreal>& _initPoint) { initPoint = _initPoint; }
    void noisifyGuessPoint(mstreal _noise = 1.0) { params.setNoise(_noise); initPoint.resize(0); }
    int numResidues() { return overlappingResidues.size(); }
    bool isAnchored() { return fixedResidues.size() > 0; }
    int numDF() {
      int df = 3*numMobileAtoms - 6;
      if (isAnchored()) df += 6;
      return df;
    }
    int getBuildOrigin() { return buildOriginRes; }
    void setBuildOrigin(int _buildOriginRes) { buildOriginRes = _buildOriginRes; }
    Structure getStructure() { return fused; }
    Structure getAlignedStructure();
    void setVerbose(bool _verbose) { params.setVerbose(_verbose); }
    int randomizeBuildOrigin() {
      buildOriginRes = (fixedResidues.size() > 0) ? fixedResidues[MstUtils::randInt(0, fixedResidues.size() - 1)] : MstUtils::randInt(0, numResidues() - 1);
      return buildOriginRes;
    }

  class icBound {
    public:
      icBound(icType _type, mstreal _minVal, mstreal _maxVal, const vector<Atom*>& _atoms, string _name = "") {
        type = _type; minVal = _minVal; maxVal = _maxVal; name = _name; atoms = _atoms;
      }
      icBound(icType _type, const pair<mstreal, mstreal>& b, const vector<Atom*>& _atoms, string _name = "") {
        type = _type; minVal = b.first; maxVal = b.second; name = _name; atoms = _atoms;
      }
      icBound(const icBound& icb) { type = icb.type; minVal = icb.minVal; maxVal = icb.maxVal; name = icb.name; atoms = icb.atoms; }
      friend ostream & operator<<(ostream &_os, const icBound& _b) {
        _os << "BOUND of type '" << _b.type << "', name '" << _b.name << "', limits [" << _b.minVal << ", " << _b.maxVal << "]";
        return _os;
      }
      mstreal getCurrentValue() const {
        switch (type) {
          case icBrokenDihedral:
          case icDihedral:
            return atoms[0]->dihedral(atoms[1], atoms[2], atoms[3]);
          case icBrokenAngle:
          case icAngle:
            return atoms[0]->angle(atoms[1], atoms[2]);
          case icBrokenBond:
          case icBond:
            return atoms[0]->distance(atoms[1]);
          default:
            MstUtils::error("uknown variable type", "icBound::getCurrentValue");
        }
      }

      vector<mstreal> getCurrentGradient() const {
        vector<mstreal> grad;
        switch (type) {
          case icBrokenDihedral:
          case icDihedral:
            grad.resize(12);
            CartesianGeometry::dihedral(atoms[0], atoms[1], atoms[2], atoms[3], grad);
            break;
          case icBrokenAngle:
          case icAngle:
            grad.resize(9);
            CartesianGeometry::angle(atoms[0], atoms[1], atoms[2], grad);
            break;
          case icBrokenBond:
          case icBond:
            grad.resize(6);
            CartesianGeometry::distance(atoms[0], atoms[1], grad);
            break;
          default:
            MstUtils::error("uknown variable type", "icBound::getCurrentGradient");
        }
        return grad;
      }

      icType type;
      mstreal minVal, maxVal;
      vector<Atom*> atoms; // atoms in the fused structure corresponding to this IC
      string name;
  };

  protected:
    AtomPointerVector atomInstances(int ri, const string& ai);
    mstreal bondInitValue(int ri, int rj, const string& ai, const string& aj, bool doNotAverage = false);
    mstreal angleInitValue(int ri, int rj, int rk, const string& ai, const string& aj, const string& ak);
    mstreal dihedralInitValue(int ri, int rj, int rk, int rl, const string& ai, const string& aj, const string& ak, const string& al);
    // for the following three functions, the index always identifies the residue of the final atom
    CartesianPoint bondInstances(int rj, Atom* atomI, Atom* atomJ, bool addToCache = true);
    CartesianPoint angleInstances(int rk, Atom* atomI, Atom* atomJ, Atom* atomK, bool addToCache = true);
    CartesianPoint dihedralInstances(int rl, Atom* atomI, Atom* atomJ, Atom* atomK, Atom* atomL, bool addToCache = true);

    // if the last flag is specified as true, will compute angular differences
    mstreal harmonicPenalty(mstreal val, const icBound& b);

    void resetScore();
    void scoreIC(const icBound& b);
    void scoreRMSD();

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

    mstreal kb, ka, kh; // force constants for enforcing bonds, angles, and dihedrals
    vector<mstreal> initPoint;
    fusionParams params;

    /* Class for storing the gradient of a Cartesian coordinate set with respect
     * to some alternative set of coordinates. E.g., could be used for storing
     * partial derivatives of Cartesian coordinates with respect to BAT
     * coordinates. Can also simply store partial derivatives of Cartesian
     * coordinates with respect to the same Cartesian coordinates, which serves
     * the purpose of defining a mapping between Atoms and corresponding
     * Cartesian coordinate indices. */
    class coordinateGradient {
      public:
        coordinateGradient() {}
        /* store the partial derivative of the dim-th Cartesian coordinate of atom
         * A (i.e., the X-, Y-, or the Z-coordinate) with respect to the alternative
         * coordinate with index coorIdx. */
        void addPartial(Atom* A, int dim, int coorIdx, mstreal partialVal) {
          partials[A][dim][coorIdx] = partialVal;
        }
        void addPartial(Atom& A, int dim, int coorIdx, mstreal partialVal) { addPartial(&A, dim, coorIdx, partialVal); }

        /* store the partial derivative of all Cartesian coordinates of atom A
         * with respect to all alternative coordinates with indices listed in
         * coorIdx (i.e., all combinations). */
        void addPartials(Atom* A, const vector<int>& coorIdx, const vector<mstreal>& partialVals) {
          if (coorIdx.size() * 3 != partialVals.size()) MstUtils::error("inconsistent number of defining variables and partial values provided", "coordinateGradient::addPartials(Atom*, const vector<int>&, const vector<mstreal>&)");
          int k = 0;
          for (int dim = 0; dim < 3; dim++) {
            for (int ci = 0; ci < coorIdx.size(); ci++) {
              addPartial(A, dim, coorIdx[ci], partialVals[k]);
              k++;
            }
          }
        }

        /* This is for the case when the Cartesian coordinates of some target atom
         * T are dependent (e.g., built from) the Cartesian coordinates of some
         * atom D, and the latter are already known within this object to be
         * dependent on some set of alternative coordinates {A}, with computed
         * partials. The function then computes and stores the partial derivatives
         * of the Cartesian coordinates of T with respect to all coordinates in {A}
         * by applying the chain rule. */
        void addRecursivePartials(Atom* T, Atom* D, const vector<mstreal>& partialVals) {
          if (partials.find(D) == partials.end()) return;
          if (partialVals.size() != 9) MstUtils::error("nine partial derivatives are expected for the dependance of two atoms", "coordinateGradient::addPartialss(Atom*, Atom*, const vector<mstreal>&)");
          int k = 0;
          for (int dimTarg = 0; dimTarg < 3; dimTarg++) {
            for (int dimDef = 0; dimDef < 3; dimDef++) {
              map<int, mstreal>& defPartials = getPartials(D, dimDef);
              for (auto it = defPartials.begin(); it != defPartials.end(); ++it) {
                addPartial(T, dimTarg, it->first, it->second * partialVals[k]);
              }
              k++;
            }
          }
        }

        int numDefiningVars(Atom* a, int dim) { return partials[a][dim].size(); }
        vector<int> getDefiningVars(Atom* a, int dim) { return MstUtils::keys(partials[a][dim]); }
        vector<mstreal> getPartialValues(Atom* a, int dim) { return MstUtils::values(partials[a][dim]); }
        map<int, mstreal>& getPartials(Atom* a, int dim) { return partials[a][dim]; }
        mstreal getPartial(Atom* a, int dim, int coorIdx) { return partials[a][dim][coorIdx]; }

        void clear() { partials.clear(); }

      private:
        map<Atom*, map<int, map<int, mstreal> > > partials;
    };

    // stores information on the gradient of Cartesian coordinates of fused atoms
    // with respect to the search coordinate vector (either also Cartesian or BAT)
    coordinateGradient gradOfXYZ;
    mstreal bondPenalty, anglPenalty, dihePenalty, rmsdScore, score; // current scores
    vector<mstreal> gradient; // current search gradient
};

class Fuser {
  public:
    /* Fuses given residues into a single polymer. The parameter resTopo specifies
     * which residues overlap. That is, resTopo.size() is the length of the
     * eventual chain, and resTopo[i] is a list of residues, of length 1 or more
     * that overlap with the i-th residue in the chain, in the N-to-C order,
     * starting with 0. Returns the the fused chain as a Structure object. */
    static Structure fuse(const vector<vector<Residue*> >& resTopo, const vector<int>& fixed = vector<int>(), const fusionParams& params = fusionParams());

    /* This function is a simplified version of Fuser::fuse(), in that it guesses
     * the topology automatically. Argument residues is a flat vector of all the
     * residues from all the different fragments that are to be fused. residues
     * from the different fragments can be intermixed in this vector, BUT when
     * one considers only residues of a given fragment, these must be in N-to-C
     * order. The function guesses the topology by finding residues that are
     * likely overlapping (i.e., have close CA atoms). This should not give
     * incorrect topologies in "normal" circumstances, but in strange cases can. */
    static Structure autofuse(const vector<Residue*>& residues, int flexOnlyNearOverlaps = -1, const fusionParams& params = fusionParams());
};


#endif
