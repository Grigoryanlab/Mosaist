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
    enum minimizerType { gradDescent = 1, conjGrad, NelderMead, langevinDyna };
    fusionParams() { // default optimization params
      startType = fusionParams::coorInitType::meanCoor;
      verbose = false;
      optimCartesian = true;
      minMethod = fusionParams::gradDescent;
      normRMSD = true; fragRedWeighting = adapRedWeighting = false;
      aBeta = 0.1;
      noise = 0;
      Ni = 100;
      Nc = 1;
      tol = 10E-8;
      kb = 10;
      ka =  0.02;
      kh = 0.001;
      krep = kcomp = 0;
      Rcomp = 1000;
      outLogBase = "";
      // dynamics-related stuff
      kT = 1;
      ts = 10E-4;
      saveSteps = -1;
      visc = 10;
      mfact = 1;
    }
    mstreal getNoise() const { return noise; }
    fusionParams::coorInitType getCoorInitType() const { return startType; }
    bool getOptimCartesian() const { return optimCartesian; }
    bool isVerbose() const { return verbose; }
    int numCycles() const { return Nc; }
    int numIters() const { return Ni; }
    mstreal errTol() const { return tol; }
    mstreal getBondFC() { return kb; }
    mstreal getAngleFC() { return ka; }
    mstreal getDihedralFC() { return kh; }
    vector<mstreal> getIntCoorFCs() { return {kb, ka, kh}; }
    mstreal getRepFC() { return krep; }
    mstreal getCompFC() { return kcomp; }
    mstreal getCompRad() { return Rcomp; }
    bool isRepOn() { return (krep != 0); }
    bool isCompOn() { return (kcomp != 0); }
    Structure getStartingStructure() const { return startStruct; }
    bool isStartingStructureGiven() const { return (startStruct.chainSize() != 0); }
    bool normalizeRMSD() const { return normRMSD; }
    bool fragRedundancyWeighting() const { return fragRedWeighting; }
    bool adaptiveWeighting() const { return adapRedWeighting; }
    int getMinimizerType() const { return minMethod; }
    string getLogBase() const { return outLogBase; }
    bool logBaseDefined() const { return !(outLogBase.empty()); }
    mstreal adaptiveBeta() const { return aBeta; }
    mstreal timeStep() const { return ts; }
    mstreal viscosity() const { return visc; }
    mstreal thermalEnergy() const { return kT; }
    mstreal massFactor() const { return mfact; }
    int saveInterval() const { return saveSteps; }

    void setNoise(mstreal _noise) { noise = _noise; }
    void setVerbose(bool _verbose = true) { verbose = _verbose; }
    void setCoorInitType(fusionParams::coorInitType _startType) { startType = _startType; }
    void setOptimCartesian(bool _optimCartesian) { optimCartesian = _optimCartesian; }
    void setNumCycles(int nc) { Nc = nc; }
    void setNumIters(int ni) { Ni = ni; }
    void setErrTol(mstreal _tol) { tol = _tol; }
    void setBondFC(mstreal _k) { kb = _k; }
    void setAngleFC(mstreal _k) { ka = _k; }
    void setIntCoorFCs(const vector<mstreal>& ks) { kb = ks[0]; ka = ks[1]; kh = ks[2]; }
    void setDihedralFC(mstreal _k) { kh = _k; }
    void setRepFC(mstreal _k) { krep = _k; }
    void setCompFC(mstreal _k) { kcomp = _k; }
    void setCompRad(mstreal _r) { Rcomp = _r; }
    void setStartingStructure(const Structure& _S);
    void clearStartingStructure() {startStruct = Structure();}
    void setNormalizeRMSD(bool flag) { normRMSD = flag; }
    void setFragRedundancyWeighting(bool flag) { fragRedWeighting = flag; }
    void setAdaptiveWeighting(bool flag) { adapRedWeighting = flag; }
    void setMinimizerType(minimizerType _type) { minMethod = _type; }
    void setLogBase(const string& _base) { outLogBase = _base; }
    void setAdaptiveBeta(mstreal beta) { aBeta = beta; }
    void setTimeStep(mstreal v) { ts = v; }
    void setViscosity(mstreal v) { visc = v; }
    void setThermalEnergy(mstreal v) { kT = v; }
    void setMassFactor(mstreal v) { mfact = v; }
    void setSaveInterval(int v) { saveSteps = v; }

  private:
    // start optimization from the averaged Cartesian structure or the structure
    // that results from average internal coordinates?
    fusionParams::coorInitType startType;
    bool verbose, optimCartesian, normRMSD, fragRedWeighting, adapRedWeighting;
    mstreal noise; // noise level for initalizing the starting point
    int Ni, Nc;    // number of iterations per cycle and number of cycles
    minimizerType minMethod; // optimization method to use
    mstreal tol;   // error tolerance stopping criterion for optimization
    mstreal kb, ka, kh; // force constants for enforcing bonds, angles, and dihedrals
    mstreal krep, kcomp; // force constants for repulsive and attractive "compactness" interactions
    mstreal Rcomp;       // desired compactness radius
    mstreal aBeta;       // (1/kT) factor used in adapting weighting to compute weights as a function of RMSD
    Structure startStruct;
    string outLogBase;   // base name for various optimization output
    // dynamics-related parameters
    mstreal kT;          // k*T for doing dynamics
    mstreal ts;          // integration time step
    int saveSteps;       // periodicity of recording dynamics snapshot
    mstreal visc;        // viscosity coefficient for Langevin dynamics; viscosity*timeStep should be no less than 10^-5
    mstreal mfact;       // the factor by whichh to multiply atomic masses in Da
};

struct fusionOutput {
  public:
    fusionOutput() { reset(); }
    fusionOutput(mstreal _bondPenalty, mstreal _anglPenalty, mstreal _dihePenalty, mstreal _rmsdScore, mstreal _rmsdTot, mstreal _score) {
      reset();
      bondPenalty = _bondPenalty; anglPenalty = _anglPenalty; dihePenalty = _dihePenalty;
      rmsdScore = _rmsdScore; rmsdTot = _rmsdTot; score = _score;
    }
    fusionOutput(const fusionOutput& sc) {
      bondPenalty = sc.bondPenalty; anglPenalty = sc.anglPenalty; dihePenalty = sc.dihePenalty;
      rmsdScore = sc.rmsdScore; rmsdTot = sc.rmsdTot; score = sc.score;
      fused = sc.fused; trajSnaps = sc.trajSnaps; trajScores = sc.trajScores;
    }
    void reset() { rmsdScore = rmsdTot = bondPenalty = anglPenalty = dihePenalty = 0; trajSnaps = vector<vector<vector<mstreal> > >(3); }
    mstreal getBondScore() const { return bondPenalty; }
    mstreal getAngleScore() const { return anglPenalty; }
    mstreal getDihedralScore() const { return dihePenalty; }
    mstreal getRMSDScore() const { return rmsdScore; }
    mstreal getTotRMSDScore() const { return rmsdTot; }
    mstreal getScore() const { return score; }
    int numSnapshots() const { return trajScores.size(); }
    Structure getSnapshot(int j) const;
    friend ostream & operator<<(ostream &_os, const fusionOutput& _s) {
      _os << _s.score << ": rmsdScore = " << _s.rmsdScore << " (RMSDtot = " << _s.rmsdTot << "), bond penalty = ";
      _os << _s.bondPenalty << ", angle penalty = " << _s.anglPenalty << ", dihedral penalty = " << _s.dihePenalty;
      return _os;
    }

    void setFused(const Structure& f) { fused = f; }
    void addSnapshot(const Structure& snap, mstreal ener);

  private:
    double rmsdScore, rmsdTot, bondPenalty, anglPenalty, dihePenalty, score;
    Structure fused;

    /* dynamics-related outputs */
    vector<vector<vector<mstreal> > > trajSnaps; // trajSnaps[0][i][j] is the x-coordinate of the j-th atom in snapshot i
    vector<mstreal> trajScores;
};

/* This class is a logical representation of the "topology" of a structure to be
 * fused. The structure has a certain length, L, and each position is indexed 0
 * through L-1. There is also a set of fragments that are to overlap this structure
 * at various locations, with each fragment represented as an array of residues
 * that "knows" which position indices it is supposed to overlap with. This also
 * allows one to specify a weight factor for each fragment, specify fixed positions,
 * and so on. */
class fusionTopology {
  friend class fusionEvaluator;
  friend class fusionParams;

  public:
    /* Initializes an empty topology of a given length */
    fusionTopology(int L = 0);

    /* initialize a topology from the "old-style" resTopo array. resTopo is of
     * length equal to the number of residues in the fused structure (i.e., the
     * length of the topology) and resTopo[i] is a list of residues that map to
     * the topology position with index i. This constructor identifies fragments
     * by looking at what Structure object each Residue belongs to. */
    fusionTopology(const vector<vector<Residue*> >& resTopo);

    /* Copy constructor */
    fusionTopology(const fusionTopology& topo);

    /* These two functions allow one to add fragments to an existing topology,
     * specifying which position indices residues of the fragment are to overlap
     * with. This is specified either explicitly, via the second parameter OR the
     * indices are taken from resnum fields of the residues (default). One can
     * also specify a weight factor for the added fragment. */
    void addFragment(Structure& S, const vector<int>& fragResIdx = vector<int>(), mstreal weight = 1.0);
    void addFragment(vector<Residue*>& R, const vector<int>& fragResIdx = vector<int>(), mstreal weight = 1.0);

    /* These function mark specified position(s) as fixed in the topology. */
    void addFixedPositions(vector<int> fixedInds);
    void addFixedPosition(int i) { fixed[i] = true; }

    /* Once a topology object is set up, this function finds pairings between the
     * fragments to be overlapped and regions in the fused structure that they
     * are to be overlapped with. Both are represented as AtomPointerVector, with
     * a one-to-one correspondence between the atoms in each pairing for easy
     * superposition and other analysis. To do this, the function accepts a reference
     * to a Structure object that represents the structure to be fused. */
    void setAlignedFrags(Structure& fused);

    /* Checks whether the given structure is consistent with the topology--i.e.,
     * has the right number of chains of the righht length. */
    bool isConsistent(const Structure& S);
    Structure remapChains(const Structure& S);

    /* Getters and setters */
    int numAlignedFrags() const { return alignedFrags.size(); }
    int numFrags() const { return fragments.size(); }
    int numUniqueOverlaps() const { return fragsByOverlap.size(); }
    int numFragsOverlapping(int i) const { return fragsByOverlap[i].size(); }
    int getFragOverlapping(int i, int k) const { return fragsByOverlap[i][k]; }
    AtomPointerVector& getAlignedFragFused(int fi) { return alignedFrags[fi].first; }
    AtomPointerVector& getAlignedFragRef(int fi) { return alignedFrags[fi].second; }
    pair<AtomPointerVector, AtomPointerVector>& getAlignedFragPair(int fi) { return alignedFrags[fi]; }
    vector<Residue*>& getOverlappingResidues(int pi) { return overlappingResidues[pi]; }
    Residue* getOverlappingResidue(int pi, int ri) { return overlappingResidues[pi][ri]; }
    int numOverlappingResidues(int pi) { return overlappingResidues[pi].size(); }
    vector<int> getFixedPositions();
    int numFixedPositions();
    int numFixedInChain(int ci);
    bool isFixed(int i) { return fixed[i]; }
    bool isFixed(int ci, int ri) { return fixedInChain[ci][ri]; }
    int length() { return overlappingResidues.size(); }
    int numMobileAtoms(int ci) { updateConnectivity(); return numMobAtoms[ci]; }
    int numMobileAtoms();
    vector<mstreal> getFragWeights() const { return fragWeights; }
    vector<int> getChainLengths() { updateConnectivity(); return chainLengths; }
    int numChains() { updateConnectivity(); return chainLengths.size(); }
    bool isVerbose() { return verbose; }
    void setVerbose(bool _verbose) { verbose = _verbose; }

    /* Assignment (copy constructor implemented using this) */
    fusionTopology& operator=(const fusionTopology& topo);

    friend ostream & operator<<(ostream &_os, const fusionTopology& _topo) {
      for (int i = 0; i < _topo.overlappingResidues.size(); i++) {
        cout << i << ": ";
        for (int k = 0; k < _topo.overlappingResidues[i].size(); k++) {
          cout << *(_topo.overlappingResidues[i][k]) << " [" << _topo.overlappingResidues[i][k]->getStructure() << "], ";
        }
        cout << endl;
      }
      return _os;
    }

  protected:
    /* This function figures out the connectivity (i.e., chains) implied by the
     * topology object. It returns the length of each chain, in the order they
     * occur within the topology. The determination is made based on available
     * fragments. If two consecutive residues are not found bonded to each other
     * in any fragment they are both a part of, then a chain break is said to
     * occur between them. */
    void updateConnectivity();

    static vector<string> bba;

  private:
    /* fragments[k] stores the correspondance between the atoms of the k-th
     * fragment and residue indices in the topology. */
    vector<pair<AtomPointerVector, vector<int> > > fragments;

    /* overlap[(set<int>) x] stores the overlap "index" corresponding to the
     * set of residues x. If this index is i, then fragsByOverlap[i] stores the
     * list of fragment indices (i.e., indices into the fragments vector above)
     * that overlap the set of residues resis. */
    map<set<int>, int> overlap;
    vector<vector<int> > fragsByOverlap;

    /* the weight of every fragment */
    vector<mstreal> fragWeights;

    /* overlappingResidues[k] is the list of residues that overlap the k-th
     * position in the topology. */
    vector<vector<Residue*> > overlappingResidues;

    /* fixed[k] indicates whether the k-th position in the topology is fixed. */
    vector<bool> fixed;

    /* alignedFrags[k] stores the correspondence between some portion of on the
     * fused structure (alignedFrags[k].first) and the corresponding atoms of
     * some fragment (alignedFrags[k].second). */
    vector<pair<AtomPointerVector, AtomPointerVector> > alignedFrags;

    // per chain information (computed when needed)
    vector<int> chainLengths;
    vector<int> numMobAtoms;
    vector<vector<bool> > fixedInChain;
    bool updated; // have changes to the topology been made that require an update of connectivity?
    bool verbose; // print connectivity information during updateConnectivity
};

class fusionEvaluator: public optimizerEvaluator {
  public:
    // broken means across a chain break
    enum icType { icBond = 1, icAngle, icDihedral, icBrokenBond, icBrokenAngle, icBrokenDihedral, icDistRep, icDistComp };

    fusionEvaluator(const fusionTopology& _topo, const fusionParams& _params = fusionParams());
    /* vector<vector<Residue*> >& resTopo simultaneously represents the various
     * aligned segments as well as which residues overlap in the alignment. It is
     * assumed that all residues from a give chain move together (i.e., are
     * part of the same aligned fragment). resTopo[i] stores the residues that
     * overlap with residue index i of the eventual fused structure. Thus, the
     * length of resTopo is also the length of the fused structure and resTopo[i]
     * has to have at least one entry for all i. */
    fusionEvaluator(const vector<vector<Residue*> >& resTopo, vector<int> fixedResidues = vector<int>(), const fusionParams& _params = fusionParams());
    void init();

    mstreal eval(const vector<mstreal>& point);
    mstreal eval(const vector<mstreal>& point, Vector& grad);
    vector<mstreal> guessPoint();
    vector<mstreal> getMasses() const { return masses; }
    void setGuessPoint(const vector<mstreal>& _initPoint) { initPoint = _initPoint; }
    void noisifyGuessPoint(mstreal _noise = 1.0) { params.setNoise(_noise); initPoint.resize(0); }
    bool isAnchored() { return (topo.numFixedPositions() > 0); }
    int numDF() {
      /* if we have fixed residues, then there is an absolute reference frame and
       * every atom gets exactly three coordinates. Otherwise, the first three
       * atoms are "special" as we want to remove rigid transformations. */
      if (isAnchored()) return 3*topo.numMobileAtoms();
      return 3*topo.numMobileAtoms() - 6;
    }
    int getBuildOrigin() { return buildOriginRes; }
    void setBuildOrigin(int _buildOriginRes) { buildOriginRes = _buildOriginRes; }
    Structure getStructure() { return fused; }
    Structure getAlignedStructure();
    void setVerbose(bool _verbose) { params.setVerbose(_verbose); }
    void chooseBuildOrigin(bool randomize = false);
    fusionOutput getScores();

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
          case icDistRep:
          case icDistComp:
            return atoms[0]->distance(atoms[1]);
          default:
            MstUtils::error("unknown variable type", "icBound::getCurrentValue");
            return 0; // to make the compiler happy
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
          case icDistRep:
          case icDistComp:
            grad.resize(6);
            CartesianGeometry::distance(atoms[0], atoms[1], grad);
            break;
          default:
            MstUtils::error("unknown variable type", "icBound::getCurrentGradient");
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
    mstreal atomRadius(const Atom& a);

    // if the last flag is specified as true, will compute angular differences
    mstreal harmonicPenalty(mstreal val, const icBound& b);

    void resetScore();
    void scoreIC(const icBound& b);
    void scoreRMSD();

  private:
    Structure fused, guess;

    /* Index of some residue, from which building based on internal coordinates
     * will procede. If there are fixed residues, one of them will be chosen. If
     * not, then some other arbitrary residue will be chosen (e.g., the first one). */
    int buildOriginRes;

    fusionTopology topo;   // stores all information about the topology of the fusion

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

    vector<mstreal> initPoint, masses;
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
    mstreal bondPenalty, anglPenalty, dihePenalty, rmsdScore, rmsdTot, score; // current scores
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
    static Structure fuse(const vector<vector<Residue*> >& resTopo, fusionOutput& scores, const vector<int>& fixed = vector<int>(), const fusionParams& params = fusionParams());
    static Structure fuse(const fusionTopology& topo, fusionOutput& scores, const fusionParams& params = fusionParams());
    static Structure fuse(const fusionTopology& topo, const fusionParams& params = fusionParams());

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
