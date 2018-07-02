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
    res_t aaToIndex(const string& aa) const;

    void buildBackgroundPotentials();
    void buildOneDimPotential(const vector<mstreal>& x, int binSpecType, const vector<mstreal>& binSpec, const vector<res_t>& aa, bool isAngle = false, const vector<mstreal>& priorPot = vector<mstreal>(), const vector<mstreal>& mult = vector<mstreal>());
    void buildTwoDimPotential(const vector<mstreal>& x, const vector<mstreal>& y, const vector<mstreal>& xBinSpec, const vector<mstreal>& yBinSpec, const vector<res_t>& aa, bool isAngle = false, const vector<mstreal>& priorPot = vector<mstreal>(), const vector<mstreal>& mult = vector<mstreal>());

    void readBackgroundPotentials(const string& file); // TODO
    void writeBackgroundPotentials(const string& file); // TODO
    mstreal selfEnergy(Residue* R, vector<Residue*> C);

  private:
    FASST F;
    string fasstdbPath, backPotFile;
    mstreal kT;
    vector<mstreal> ppPot, omPot, envPot;
    /* Stores any mapping between natural amino acid names and "standard" amino
     * acids. The idea is that we may want to interpret SEC (selenocysteine) (for
     * example) as CYS (cysteine) in gathering sequence statistics. But we also
     * may want to keep these as separate counts. This map stores any such
     * correspondances we may want to use, while aaMapType reflects which mapping
     * type is being used. The latter is interpreted by setAminoAcidMap(). */
    map<string, string> aaMap;
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
     * equilibration steps is not specified, it is defaulted to ceil(0.2*Ni). */
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
