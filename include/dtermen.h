#ifndef _DTERMEN_H
#define _DTERMEN_H

#include "msttypes.h"
#include "mstfasst.h"
#include "mstsequence.h"
using namespace MST;

class dTERMen {
  public:
    // Implementation approach:
    // 1. DONE: get phi, psi, and environment values into the database (fasstDB does this now)
    // 2. DONE: perform all-by-all local sequence window comparison within the database, using the
    //    uclust-like algorithm (not guaranteed, but pretty good for close homology), and get
    //    two things for every position: its multiplicity (i.e., the number of times closely
    //    homologous local sequences are found in the database) and the local cluster index
    //    (i.e., such that all closely homologous local sequence windows belong to the same
    //    cluster).
    // 3. DONE: make a modification in FASST, so that it can optionally use the "local cluster
    //    index" property to remove redundancy (no need to calculate alighments, so very fast),
    //    if the property is available within the database.
    // 4. TODO: implement on-the-fly backbone and environment potential calculation (upon reading
    //    the FASST database), by accounting for multiplicity of each residue.
    // TODO:
    // 1. get location of FASST database and initialize a built-in FASST object
    // 2. from data within FASST database, build phi, psi, omega, and environment potential lookup tables
    dTERMen(const string& configFile);
    dTERMen();
    void readConfigFile(const string& configFile);

    mstreal getkT() const { return kT; }
    void setAminoAcidMap();
    void printAminoAcidMap();
    int globalAlphabetSize() const { return globAlph.size(); }

    bool isInGlobalAlphabet(const string& aa) const;
    bool isInGlobalAlphabet(res_t aa) const;
    // the following four functions convert to and from the "internal" amino-acid index
    int aaToIndex(const string& aa) const;
    int aaToIndex(res_t aa) const;
    string indexToResName(int idx) const;
    res_t indexToAA(int idx) const;

    struct histType {
      vector<vector<int> > bins;
      vector<mstreal> binEdges;
      vector<mstreal> binMasses;
      vector<mstreal> weights;
    };
    struct oneDimPotType {
      vector<mstreal> binEdges;
      vector<vector<mstreal> > aaEnergies;
    };
    struct twoDimPotType {
      vector<mstreal> xBinEdges;
      vector<mstreal> yBinEdges;
      vector<vector<mstreal> > aaEnergies;
    };
    void buildBackgroundPotentials();
    oneDimPotType buildOneDimPotential(const histType& H, const vector<int>& AA, mstreal pc, vector<vector<mstreal> >& backPot, bool updateBackPot = false);
    oneDimPotType buildOneDimPotential(const histType& H, const vector<int>& AA, mstreal pc); // without specifying a background potential
    twoDimPotType buildTwoDimPotential(const vector<mstreal>& x, const vector<mstreal>& y, const vector<mstreal>& xBinSpec, const vector<mstreal>& yBinSpec, const vector<int>& AA, bool isAngle = false, const vector<mstreal>& priorPot = vector<mstreal>(), const vector<mstreal>& mult = vector<mstreal>());
    mstreal lookupOneDimPotential(const oneDimPotType& P, mstreal x, res_t aa);
    mstreal lookupOneDimPotential(const oneDimPotType& P, mstreal x, mstreal y, res_t aa);
    void printOneDimPotential(const oneDimPotType& P);
    void printTwoDimPotential(const twoDimPotType& P); // TODO

    /* Bins the data in the input vector X according to the binning type and
     * parameters. Outputs a struct with members bins and binEdges. The k-th bin
     * (k being between 0 and n-1, where n is the number of bins) is defined as
     * the interval [binEdges[k]; binEdges[k+1]), except for the last bin, which
     * does include the right limit: [binEdges[k]; binEdges[k+1]]. Upon returning,
     * bins[k] will contain the list of indices of data points that map into the
     * k-th bin. Supported binning types are (values for binSpecType):
     * 1 -- uniform binning. binSpec is expected to be {min value, max value,
     *      and desired number of bins}.
     * 2 -- non-uniform binning with some minimal number of elements per bin and
     *      a minimal bin width. binSpec is expected to be {min value, max value,
     *      minimum number of points per bin, minimal bin width}.
     * Optional argument M (must be the same size as X, if specified) supplies
     * the multiplicity of each point. This is used in non-uniform binning to
     * count data "mass" (i.e., sum of inverses of multiplicities) rather the
     * pure number of points. If isAngle is set to true, will treat input data
     * as angles in degrees, mapping them to the interval [-pi; pi). This makes
     * it so that the bin definitions above do the right thing of counting +/- pi
     * just once.
     * NOTE: any data points that fall outside of the range [min; max] are
     * ignored in building the potential. This may include bad diehdral angles,
     * such as phi/psi angles for terminal residues. */
    histType binData(const vector<mstreal>& X, int binSpecType, const vector<mstreal>& binSpec, const vector<mstreal>& M = vector<mstreal>(), bool isAngle = false);

    void readBackgroundPotentials(const string& file); // TODO
    void writeBackgroundPotentials(const string& file); // TODO

    mstreal selfEnergy(Residue* R, vector<Residue*> C);
    mstreal bbOmegaEner(mstreal omg, res_t aa) { return lookupOneDimPotential(omPot, omg, aa); }
    mstreal envEner(mstreal env, res_t aa) { return lookupOneDimPotential(envPot, env, aa); }

  private:
    FASST F;
    string fasstdbPath, backPotFile;
    mstreal kT;
    oneDimPotType omPot, envPot;
    twoDimPotType ppPot;

    /* We may want to deal with different "universal" alphabets (separate from
     * the design alphabet). We may want to interpret SEC (selenocysteine), for
     * example, as CYS (cysteine) in gathering sequence statistics. Or, we may
     * want to keep these as separate counts. The following several variables
     * support this capability by providing a mapping between all the possible
     * amino acid names and indices defined in SeqTools to a smaller sub-set,
     * defined over a contiguous set of indices starting with 0. In all cases
     * residues are identified by their res_t index, and the classification
     * between residue strings and res_t is left to SeqTools. But, for simplicity,
     * in the comments below, I will refer to residues by their string.
     * aaMap -- stores any mapping between any non-standard amino acid names and
     *          standard ones. E.g., aaMap["SEC"] may contain "CYS", saying that
     *          SEC is interpreted as a CYS, or it may not exist, which says that
     *          we want to explicitly track the counts of SEC, separate from CYS.
     * globAlph -- a list of all counted amino acids, referred to by their most
     *             "standard" name. E.g., if aaMap["SEC"] == "CYS", then globAlph
     *             will contain "CYS" but not "SEC". Further, amino-acid names
     *             in globAlph are stored according to their contiguous index
     *             (internal ones for dTERMen, not the same as SeqTools indices).
     * aaIdx    -- effectively the opposite of globAlph. For every counted amino
     *             acid, aaIdx will contain its corresponding internal index. E.g.,
     *             aaIdx["SEC"] and aaIdx["CYS"] would have the same index, if
     *             if aaMap["SEC"] == "CYS".
     * aaMapType -- specifies the type of mapping between non-standard amino
     *              acids and their standard counterparts. Used by setAminoAcidMap() */
    map<res_t, res_t> aaMap;
    vector<res_t> globAlph;
    map<res_t, int> aaIdx;
    int aaMapType;

};

class EnergyTable {
  public:
    EnergyTable(const string& tabFile);
    void readFromFile(const string& tabFile);
    void addSite(const string& siteName); // TODO
    void addSites(const vector<string>& siteNames); // TODO
    int numSites() const { return siteIndices.size(); }

    // -- some simple evaluation routines
    mstreal selfEnergy(int s, int aa);
    mstreal pairEnergy(int si, int sj, int aai, int aaj);
    mstreal scoreSolution(const vector<int>& seq);
    mstreal scoreSequence(const Sequence& seq);
    mstreal scoreMutation(const vector<int>& seq, int mutSite, int mutAA);
    mstreal scoreMutation(const vector<int>& seq, const vector<int>& mutSites, const vector<int>& mutAAs);
    mstreal meanEnergy() const;
    mstreal energyStdEst(int n = 1000);

    // -- optimization routines
    vector<int> randomSolution() const;
    int randomResidue(int si) const;

    /* Monte Carlo sampling/optimization, with optional annealing. Parameters:
     * Nc   -- number of cycles (each cycle starts with a new random solution)
     * Ni   -- number of iterations per cycle (each iteration makes a move and
     *        either accepts or rejects)
     * kTi  -- initial temperature
     * kTf  -- final temperature. If specified, will anneal from kTi to kTf
     *         in each cycle.
     * annealType -- type of annealing to use. Currently supported are 1 for
     *               linear and 2 for exponential. Linear is the default.
     * rec  -- a pointer to some trajectory recorder container. Only used if the
     *         next argument is also specified (see below).
     * add  -- a function pointer that will be called (if specified) for every
     *         accepted solutions (note, a rejected move results in the previous
     *         sequence being accepted). Specifically, the call will be:
     *         (*add)(rec, seq, ener), where sol is the current solution and ener
     *         is its energy, such that variable rec can keep track of the MC
     *         trajectory in an entirely customizable way.
     * Ne   -- number of pre-equilibration steps to do in each cycle before
     *         beginning to record encountered solutions.
     * Returns the lowest-energy sequence encountered. NOTE: each cycle will
     * involve an initial equilibration phase, during which iterations are
     * are performed, but the trajectory is not recorded. If the number of pre-
     * equilibration steps is not specified, it is defaulted to ceil(0.2*Ni).
     * TODO: add a two-residue mutation step (for pairs that have pair energies)
     * and choose that move some fraction of the time. Could choose the pair to
     * mutate based on the variance of interaction strengths at the pair. */
    vector<int> mc(int Nc, int Ni, mstreal kTi, mstreal kTf = -1, int annealType = 1, void* rec = NULL, void (*add)(void*, const vector<int>&, mstreal) = NULL, int Ne = -1);

    Sequence solutionToSequence(const vector<int>& sol);
    vector<int> sequenceToSolution(const Sequence& seq);
    string getResidueString(int si, int ri); // TODO
    int getResidueIndex(int si, const string& aa) { return aaIndices[si][aa]; }

  private:
    /* siteIndices["A,1"] is the indices corresponding to site "A,1" (or however
     * sites are designated). */
    map<string, int> siteIndices;

    /* aaIndices[0]["ALA"] is the index corresponding to "ALA" at the first
     * encountered site. */
    vector<map<string, int> > aaIndices;

    /* indToAA[0][1] is the amino-acid string (e.g., "ALA") that corresponds to
     * residue index 1 at the first site. */
    vector<map<int, string> > indToAA;

    /* selfE[0][3] is the self energy for the amino acid with index 3 at site 0 */
    vector<vector<mstreal> > selfE;

    /* pairs[i] -- list of site indices that site i has interactions with */
    vector<vector<int> > pairs;

    /* pairsMap[i] -- same as above, but the list is stored as an int->bool map,
     * and every pair entry exists with either true (the pair does have
     * interactions) or false (the pair does not have interactions). */
    vector<map<int, bool> > pairMaps;

    /* pairE[i][pairs[i][j]][aai][aaj] -- pair energy for sites (i -- pairs[i][j])
     * occupied with aai-th and aaj-th amino acid in each site, respectively,
     * assuming that i < pairs[i][j]. NOTE: this assumes energies are stored in
     * one direction only! */
    vector<vector<vector<vector<mstreal > > > > pairE;
};

#endif
