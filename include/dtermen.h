#ifndef _DTERMEN_H
#define _DTERMEN_H

#include "msttypes.h"
#include "mstfasst.h"
#include "mstsequence.h"
#include "mstcondeg.h"
#include "mstmagic.h"
using namespace MST;

class EnergyTable;
class termData;
/* Remaining questions:
 * TODO: when there is no +/- pmSelf or +/- pmPair on one side, extend on the other side (current dTERMen)!
 * TODO: make terminus-specific trivial potentials (probably all of them)
 * TODO: A MAJOR issue is that there is considerable overcounting of the local backbone between self residual
 *       and self correction. Possible solutions are:
 *       - pre-compute self residual for all residues in the database; would take a long time, but could be
 *         shortened considerably via FASST cache. With a pre-computed self residual we can always take it
 *         out of the contributions from the mathces (just like we do for phi/psi).
 *   IDEA: instead of splitting into residual, correct, and pair, simply do like we did in the PLoS One paper.
 *         Gather lots of TERM matches (split, but keep as much of the context as possible) and then optimize
 *         19 self energies for each site and 19x19 pair energies for each pair. Do this one site at a time,
 *         meaning that each split into sub-TERMs happens specifically for one site. Can either keep the method
 *         in the paper or do some other split, but the underlying function to tease out self and pair energies
 *         would be the same. One issue of that method is that when a residue is missing in some sub-TERM, it may
 *         still be present in many of its matches. Though by the nature of our decomposition, it should not be
 *         present in too many of the matches, because we actually try to keep as large of a TERM as possible.
 */

class dTERMen {
  public:
    dTERMen();
    dTERMen(const string& configFile);
    void init();
    void readConfigFile(const string& configFile);
    void setEnergyFunction(const string& ver); // sets the energy function version and alters any necessary parameters
    void setRecordFlag(bool record = true) { recordData = record; }

    /* Builds an energy table for design. Parameters:
     * - a list of mutable positions as vector<Residue*>. All residues must belong
     *   to a single Structure objects.
     * - (optional) a list of allowed amino acids at each position. If not given,
     *   the entire amino-acid alphabet known to the dTERMen object will be
     *   allowed at each site.
     * - (optional) a list of images for crystal symmetry design. If given, the
     *   "central" unit cell must be listed first, followed by all other images
     *    that must have exactly the same number of residues listed in the same
     *    order. The list of variable positions must come from the central unit
     *    cell. Similarly, if the specificity context is specified (see below),
     *    this selection must also be from the central unit cell.
     * - (optional) a pointer to an EnergyTable object that will be filled with
     *   "specificity gap" energies. This is a table that quantifies the extent
     *   to which each amino-acid choice at each variable position is specific to
     *   the fixed context (i.e., the specific decoration present at the fixed
     *   sites--i.e., all those that are not variable).
     * - (optional) a list of residues to consider as the fixed context for the
     *   purposes of specificity-gap energy calculation. This can be used if not
     *   all of the fixed residues are to be considered as the "context" for the
     *   purposes of specificity calculation. For example, this could be the list
     *   of sites a fixed interaction partner for which we are trying to design
     *   a binder. While some of the binder positions may also be fixed, they are
     *   not part of the specificity context.
     */
    EnergyTable buildEnergyTable(const vector<Residue*>& variable, const vector<vector<string>>& allowed = vector<vector<string>>(), const vector<vector<Residue*>>& images = vector<vector<Residue*>>(), EnergyTable* specTable = NULL, const vector<Residue*>& specContext = vector<Residue*>());

    mstreal getkT() const { return kT; }
    FASST* getFASST() { return &F; }
    RotamerLibrary* getRotamerLibrary() { return &RL; }
    vector<res_t> getGlobalAlphabet() { return globAlph; }
    void setAminoAcidMap();
    void printAminoAcidMap();
    int globalAlphabetSize() const { return globAlph.size(); }
    mstreal getHomologyCutoff() const { return homCut; }

    bool isInGlobalAlphabet(const string& aa) const;
    bool isInGlobalAlphabet(res_t aa) const;
    // the following four functions convert to and from the "internal" amino-acid index
    int aaToIndex(const string& aa) const;
    int aaToIndex(res_t aa) const;
    bool aaIndexKnown(int aaIdx) const;
    string indexToResName(int idx) const;
    res_t indexToAA(int idx) const;

    struct oneDimHist {
      vector<vector<int> > bins;
      vector<mstreal> binEdges;
      vector<mstreal> binMasses;
      vector<mstreal> weights;
    };
    struct twoDimHist {
      vector<vector<vector<int> > > bins;
      vector<mstreal> xBinEdges, yBinEdges;
      vector<mstreal> weights;
    };
    struct zeroDimPotType { // NOTE: in the future, these should become classes with overloaded "access" operators for lookup
      vector<mstreal> aaEnergies;
    };
    struct oneDimPotType {
      vector<mstreal> binEdges;
      vector<vector<mstreal> > aaEnergies;
    };
    struct twoDimPotType {
      vector<mstreal> xBinEdges;
      vector<mstreal> yBinEdges;
      vector<vector<vector<mstreal> > > aaEnergies;
    };
    void buildBackgroundPotentials();
    zeroDimPotType buildZeroDimPotential(const vector<int>& AA, vector<vector<mstreal> >& backPot);
    oneDimPotType buildOneDimPotential(const oneDimHist& H, const vector<int>& AA, mstreal pc, vector<vector<mstreal> >& backPot, bool updateBackPot = false);
    oneDimPotType buildOneDimPotential(const oneDimHist& H, const vector<int>& AA, mstreal pc); // without specifying a background potential
    twoDimPotType buildTwoDimPotential(const twoDimHist& H, const vector<int>& AA, mstreal pc, vector<vector<mstreal> >& backPot, bool updateBackPot = false);
    mstreal lookupZeroDimPotential(const zeroDimPotType& P, int aa);
    mstreal lookupOneDimPotential(const oneDimPotType& P, mstreal x, int aa);
    mstreal lookupTwoDimPotential(const twoDimPotType& P, mstreal x, mstreal y, int aa);
    void printZeroDimPotential(const zeroDimPotType& P);
    void printOneDimPotential(const oneDimPotType& P);
    void printTwoDimPotential(const twoDimPotType& P);

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
    oneDimHist binData(const vector<mstreal>& X, int binSpecType, const vector<mstreal>& binSpec, const vector<mstreal>& M = vector<mstreal>(), bool isAngle = false);
    twoDimHist binData(const vector<mstreal>& X, const vector<mstreal>& Y, const vector<mstreal>& xBinSpec, const vector<mstreal>& yBinSpec, const vector<mstreal>& M = vector<mstreal>(), bool isAngle = false);

    void readBackgroundPotentials(const string& file);  // TODO
    void writeBackgroundPotentials(const string& file); // TODO
    void writeRecordedData(const string& file); // writes all recorded TERM data to a file

    mstreal backEner(const string& aa) { return lookupZeroDimPotential(bkPot, aaToIndex(aa)); }
    mstreal bbOmegaEner(mstreal omg, const string& aa) { return lookupOneDimPotential(omPot, omg, aaToIndex(aa)); }
    mstreal bbPhiPsiEner(mstreal phi, mstreal psi, const string& aa) { return lookupTwoDimPotential(ppPot, phi, psi, aaToIndex(aa)); }
    mstreal envEner(mstreal env, const string& aa) { return lookupOneDimPotential(envPot, env, aaToIndex(aa)); }

    /* Compute the dTERMen self energy of a given Residue; identifies relevant
     * contacts using ConFind. The first form computes the value for a specific
     * amino acid at the given position (by default the one actually there), the
     * second form returns values for all amino acids at the position, and the
     * third form works as the second, but accepts a pre-built ConFind object.
     * When multiple energies are needed from the same structure and, passing the
     * same ConFind object will make it so that contacts and freedoms will only
     * need to be computed once. */
    mstreal selfEnergy(Residue* R, const string& aa = "");
    vector<mstreal> selfEnergies(Residue* R, bool verbose = false);
    vector<mstreal> selfEnergies(Residue* R, ConFind& C, bool verbose = false);

    mstreal pairEnergy(Residue* Ri, Residue* Rj, const string& aai = "", const string& aaj = "");
    vector<vector<mstreal>> pairEnergies(Residue* Ri, Residue* Rj, bool verbose = false);
    vector<vector<mstreal>> pairEnergiesNew(Residue* Ri, Residue* Rj, bool verbose = false);
    vector<vector<mstreal>> pairEnergiesNew2(Residue* Ri, Residue* Rj, bool verbose = false);

  protected:
    int findBin(const vector<mstreal>& binEdges, mstreal x); // do a binary search to find the bin into which the value falls
    mstreal backEner(int aai) { return lookupZeroDimPotential(bkPot, aai); }
    mstreal bbOmegaEner(mstreal omg, int aai) { return lookupOneDimPotential(omPot, omg, aai); }
    mstreal bbPhiPsiEner(mstreal phi, mstreal psi, int aai) { return lookupTwoDimPotential(ppPot, phi, psi, aai); }
    mstreal envEner(mstreal env, int aai) { return lookupOneDimPotential(envPot, env, aai); }
    void printSelfComponent(const CartesianPoint& ener, const string& prefix);

    /* Given a list of FASST solutions, computes the residual statistical energy
     * for all amino acids at the position with index cInd, after accounting for
     * all "trivial" background statistical contributions at this position
     * accross all of the matches. */
    CartesianPoint singleBodyStatEnergy(fasstSolutionSet& matches, int cInd, int pc);

    /* For a given match, compute the expectation of any given amino acid at the
     * specified position based on "trivial" background statistical contributions. */
    CartesianPoint backExpectation(const fasstSolution& m, int cInd);

    /* Counts observations at the given positions across all matches. */
    CartesianPoint singleBodyObservations(fasstSolutionSet& matches, int cInd);

    /* Same as above, but for amino-acid observations at two sites. */
    CartesianPoint twoBodyObservations(fasstSolutionSet& matches, int cIndI, int cIndJ);

    /* Computes the number of times each amino acid is expected to be found at
     * the given position, across all matches. NOTE: skips any matches that have
     * an unknown residue in the position. This is because when comparing to the
     * observed counts, such matches would have no "vote", so it would be unfair
     * to given them a vote in computing the expectation. The optional pointer
     * can be specified to collect the expectation in each match. */
    CartesianPoint singleBodyExpectations(fasstSolutionSet& matches, int cInd, vector<CartesianPoint>* breakDown = NULL);

    /* Computes the number of times each amino acid is expected to be found at
     * the given position in each match. Identifies an underlying amino-acid bias
     * energy vector to make sure that the marginals (i.e., the total number of
     * times each amino acid is expected across all matches) approximately equal
     * to the expectations. In practice, the agreement should be essentially
     * perfect, limited only by the number of iterations of the underlying itera-
     * tive procedure. */
    vector<CartesianPoint> singleBodyExpectationsMatchedMarginals(fasstSolutionSet& matches, int cInd);

    mstreal enerToProb(vector<mstreal>& ener);
    mstreal enerToProb(const vector<mstreal>& _ener) { vector<mstreal> ener = _ener; return enerToProb(ener); }

    /* Amino-acid indices provide a convenient way to index into a array, while
     * also encoding the amino-acid type. For amino-acid pairs, this is tricky,
     * because a 2D array is needed. These functions convert between two amino-
     * acid indices and a single index that uniquely identifies the pair. */
    int pairToIdx(int aai, int aaj) const;
    pair<int, int> idxToPair(int idx) const;

    /* Given a list of source residues, and a ConFind object (already initialized
     * with the relevant structure), returns a list of residue pairs that corres-
     * pond to contacts with source residues. These are either intra contacts--,
     * contacts formed between the source residues--or contacts between source
     * residues and residues not belonging to that list. This is controlled by the
     * last integer parameter (0, the default, corresponds to inter, 1 to intra,
     * and two to both). In the case of inter-contacts, the source residue of
     * each contact will be listed first in each pair (i.e., will be stored in
     * the "first" field of each pair). The contacts in the returned list will
     * be ordered as follows: first contact-degree based contacts, ordered by
     * contact degree in descending order, then any additional contacts emerging
     * from interference, ordered by interference. */
    vector<pair<Residue*, Residue*>> getContactsWith(const vector<Residue*>& source, ConFind& C, int type = 0, bool verbose = false);

  private:
    FASST F;
    fasstSearchOptions foptsBase; // base FASST options that will be used with every search
    RotamerLibrary RL;
    string fasstdbPath, backPotFile, rotLibFile, efunVer;
    zeroDimPotType bkPot;
    oneDimPotType omPot, envPot;
    twoDimPotType ppPot;
    mstreal kT, cdCut, intCut, selfResidualPC, selfCorrPC, homCut;
    int pmSelf, pmPair;
    int selfResidualMinN, selfResidualMaxN, selfCorrMinN, selfCorrMaxN, selfCorrMaxCliqueSize, pairMinN, pairMaxN;
    bool recordData;
    vector<termData> data;
    Sequence targetOrigSeq;
    vector<int> variableResidues;
    map<string, vector<mstreal>> targetResidueProperties;

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

/* Helper class that stores data about a single TERM searched in the process of the dTERMen
 * energy-table calculation. */
class termData {
  public:
    termData() {}
    termData(vector<Residue*> _centResidues, int pm) {
      vector<int> cenResFlankingRes(_centResidues.size(),pm);
      define(_centResidues, pm);
    }
    void define(vector<Residue*> _centResidues, int pm) {
      centResidues = _centResidues;
      term.reset();
      fragResIdx.clear();
      centResIndices = TERMUtils::selectTERM(centResidues, term, pm, &fragResIdx);
    }
    void addCentralResidue(Residue* res, int pm) {
      centResidues.push_back(res);
      define(centResidues, pm);
    }

    // void setMatches(const fasstSolutionSet& _matches) { matches = _matches; }
    int setMatches(const fasstSolutionSet& _matches, mstreal homologyCutoff, FASST* F = NULL);
    fasstSolutionSet& getMatches() { return matches; }
    fasstSolution& getMatch(int i) { return matches[i]; }
    int numMatches() const { return matches.size(); }
    const Structure& getTERM() const { return term; }
    const vector<int>& getResidueIndices() const { return fragResIdx; }
    const vector<Residue*>& getCentralResidues() const { return centResidues; }
    const vector<int>& getCentralResidueIndices() const { return centResIndices; }
    int numCentralResidues() const { return centResidues.size(); }

    friend ostream & operator<<(ostream &_os, const termData& _td) {
      _os << "TERM clique with " << _td.centResidues.size() << " central residues, " << _td.matches.size() << " matches:";
      for (int i = 0; i < _td.centResidues.size(); i++) _os << " " << *(_td.centResidues[i]);
      return _os;
    }
    string toString() { stringstream ss; ss << *this; return ss.str(); }

  private:
    /* the central residues of each of the TERM's segments and their corres-
     * ponding indices in the template. */
    vector<Residue*> centResidues;
    /* the flanking residues around each of the central residues */
    vector<int> cenResFlankingRes;
    vector<int> centResIndices;

    /* term is the TERM itself, whose residues correspond to template posi-
     * tions at indices fragResIdx. */
    Structure term;
    vector<int> fragResIdx;

    /* matches resulting from querying the TERM. */
    fasstSolutionSet matches;
};

class EnergyTable {
  public:
    EnergyTable() {}
    EnergyTable(const string& tabFile);

    /* restrictSiteAlphabet() constructs a new energy table, copying only the residue types that are
     specified each position in restricted_siteAlphabets. When the new energy table will be applied
     as a constraint to the original one, constraint_table should be set to true. In this case, all
     residue types will be carried over, but their energies will correspond to the energy of the
     residue type specified in restricted_siteAlphabets. Note that in this mode, no more than one
     residue type can be specified per position. If a position in restricted_siteAlphabets is empty,
     all residue types will be copied.
     */

    EnergyTable restrictSiteAlphabet(const vector<vector<string>>& restricted_siteAlphabets, bool constraint_table = false);
    EnergyTable restrictSiteAlphabet(const Structure& S, bool constraint_table = false);

    void clear(); // resets the table to empty
    void readFromFile(const string& tabFile);
    void writeToFile(const string& tabFile);
    void addSite(const string& siteName);
    void addSites(const vector<string>& siteNames);
    int numSites() const { return siteIndices.size(); }
    vector<string> getSites() {return sites; }
    int siteIndex(const string& siteName) { return siteIndices[siteName]; }
    bool siteExists(const string& siteName) { return siteIndices.find(siteName) != siteIndices.end(); }
    void setSiteAlphabet(const string& siteName, const vector<string>& alpha) { setSiteAlphabet(siteIndex(siteName), alpha); }
    void setSiteAlphabet(int siteIdx, const vector<string>& alpha);
    void renameSiteResidue(int si, int ai, const string& aa);
    vector<string> getSiteAlphabet(int siteIdx) { return aaAlpha[siteIdx]; }
    int addToSiteAlphabet(int siteIdx, const string& aa);
    bool inSiteAlphabet(int siteIdx, const string& aa) { return aaIndices[siteIdx].find(aa) != aaIndices[siteIdx].end(); }
    int indexInSiteAlphabet(int siteIdx, const string& aa) { return aaIndices[siteIdx][aa]; }
    bool empty() const { return selfE.empty() && pairE.empty(); }

    // -- get/set energy-table components
    mstreal selfEnergy(int s, int aa);
    mstreal pairEnergy(int si, int sj, int aai, int aaj);
    void setSelfEnergy(int s, int aa, mstreal ener);
    void setPairEnergy(int si, int sj, int aai, int aaj, mstreal ener);

    // -- some simple evaluation routines
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
     * additionalScore -- a pointer to a function that returns the additional
     *                    score component for any solution. The last two arguments
     *                    to the function are used for scoring mutations. The
     *                    function is free to implement these as it wants, even
     *                    if without any inherent savings.
     * Returns the lowest-energy sequence encountered. NOTE: each cycle will
     * involve an initial equilibration phase, during which iterations are
     * are performed, but the trajectory is not recorded. If the number of pre-
     * equilibration steps is not specified, it is defaulted to ceil(0.2*Ni).
     * TODO: add a two-residue mutation step (for pairs that have pair energies)
     * and choose that move some fraction of the time. Could choose the pair to
     * mutate based on the variance of interaction strengths at the pair. */
    vector<int> mc(int Nc, int Ni, mstreal kTi, mstreal kTf = -1, int annealType = 1, void* rec = NULL, void (*add)(void*, const vector<int>&, mstreal) = NULL, int Ne = -1, void* extra = NULL, mstreal (*additionalScore)(void*, const vector<int>&, EnergyTable&, int mutSite, int mutAA) = NULL);

    Sequence solutionToSequence(const vector<int>& sol);
    vector<int> sequenceToSolution(const Sequence& seq, bool strict = false);
    string getResidueString(int si, int ri); // TODO
    int getResidueIndex(int si, const string& aa) { return aaIndices[si][aa]; }

  private:
    /* siteIndices["A,1"] is the index corresponding to site "A,1" (or however
     * sites are designated). */
    map<string, int> siteIndices;

    /* sites[k] is the name of the site with index k. */
    vector<string> sites;

    /* aaIndices[0]["ALA"] is the index corresponding to "ALA" at the first
     * encountered site. */
    vector<map<string, int> > aaIndices;

    /* aaAlpha[0][1] is the amino-acid string (e.g., "ALA") that corresponds to
     * residue index 1 at the first site. */
    vector<vector<string> > aaAlpha;

    /* selfE[0][3] is the self energy for the amino acid with index 3 at site 0 */
    vector<vector<mstreal> > selfE;

    /* pairMaps[i] -- a map, whose keys represent all site indices with which
     * site i has interactions with. If k is some such site, and i < k, then
     * pairMaps[i][k] is the index into pairE[i] at which interactions between
     * sites i and k are stored. If, i > k, then pairMaps[i][k] is the index
     * into pairE[k], at which interactions between sites k and i are stored. */
    vector<map<int, int> > pairMaps;

    /* pairE[i][pairMaps[i][k]][aai][aaj] -- pair energy for site i occupied with
     * amino acid aai interacting with site k occupied with amino acid aaj,
     * assuming that i < pairMaps[i][k]. */
    vector<vector<vector<vector<mstreal > > > > pairE;
};

#endif
