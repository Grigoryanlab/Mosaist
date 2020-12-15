#ifndef _MSTCONDEG_H
#define _MSTCONDEG_H

#include "msttypes.h"
#include "mstrotlib.h"
#include <set>

using namespace std;
using namespace MST;

template<typename key, typename T>
using fastmap = map<key, T>;

class interactionList {
  public:
    int size() { return resi.size(); }
    Residue* residueA(int i) { return resi[i]; }
    Residue* residueB(int i) { return resj[i]; }
    Residue* srcResidue(int i) { return resi[i]; }
    Residue* dstResidue(int i) { return resj[i]; }
    vector<Residue*> srcResidues() { return resi; }
    vector<Residue*> destResidues() { return resj; }
    mstreal degree(int i) { return degrees[i]; }
    string info(int i) { return infos[i]; }
    void sortByDegree(); // sorts the contact list by contact degree, highest to lowest
  
  protected:
    vector<Residue*> resi;
    vector<Residue*> resj;
    vector<mstreal> degrees;
    vector<string> infos;
};

class contactList: public interactionList {
  public:
    void addContact(Residue* _resi, Residue* _resj, mstreal _degree, string _info = "", bool directional = false);
  
    mstreal degree(Residue* _resi, Residue* _resj);
    vector<pair<Residue*, Residue*>> getOrderedContacts();
    bool areInContact(Residue* _resi, Residue* _resj);

  protected:
    struct contComp {
      bool operator() (const pair<Residue*, Residue*>& lhs, const pair<Residue*, Residue*>& rhs) const {
        int lhsI = lhs.first->getResidueIndex();
        int rhsI = rhs.first->getResidueIndex();
        if (lhsI == rhsI) {
          lhsI = lhs.second->getResidueIndex();
          rhsI = rhs.second->getResidueIndex();
        }
        return lhsI < rhsI;
      }
    };
    fastmap<Residue*, fastmap<Residue*, int> > inContact;
    set<pair<Residue*, Residue*>, contComp> orderedContacts;
};

class aaConstrainedContactList: public interactionList {
  public:
    void addContact(Residue* _resi, Residue* _resj, mstreal _degree, set<string> _resi_aa = {}, set<string> _resj_aa = {}, string _info = "", bool directional = false) {
      resi.push_back(_resi);
      resj.push_back(_resj);
      degrees.push_back(_degree);
      infos.push_back(_info);
      resi_aa.push_back(_resi_aa);
      resj_aa.push_back(_resj_aa);
    }
  
    set<string> residueA_aa(int i) {return resi_aa[i];}
    set<string> residueB_aa(int i) {return resj_aa[i];}
  protected:
    vector<set<string>> resi_aa;
    vector<set<string>> resj_aa;
};

class ConFind {
  public:
    ConFind(string rotLibFile, const Structure& S, bool strict = false);
    ConFind(RotamerLibrary* _rotLib, const Structure& S, bool strict = false);
    ~ConFind();
    void setFreedomParams(mstreal _loCollProbCut, mstreal _hiCollProbCut, int type) { loCollProbCut = _loCollProbCut; hiCollProbCut = _hiCollProbCut; freedomType = type; }

    // precomputes all necessary info and data structures for computing on this Structure
    void cache(const Structure& S);
    void cache(const vector<Residue*>& residues);
    void cache(Residue* res);

    // find those residues that are close enough to affect the passed residue(s)
    vector<Residue*> getNeighbors(Residue* residue);
    vector<Residue*> getNeighbors(vector<Residue*>& residues);
    bool areNeighbors(Residue* resA, Residue* resB);

    /* this function encodes whether a given atom counts as "side-chain" for the
     * purposes of finding sidechain-to-sidechain contacts. */
    bool countsAsSidechain(Atom& a);
  
    /*
     Contact degree is the potential for the sidechain of two residues to interact. The A_aa and B_aa
     options allow the user to specify the amino acid rotamers considered at residue A and B, respectively.
     If nothing is specified, it is assumed that all are allowed.
     
     Note: specifying amino acids that are not in the library (e.g. GLY, PRO, non-canonical) will
     result in an error.
     */

    mstreal contactDegree(Residue* resA, Residue* resB, bool cacheA = true, bool cacheB = true, bool checkNeighbors = true, set<string> aaAllowedA = {}, set<string> aaAllowedB = {});
    contactList getContacts(Residue* res, mstreal cdcut = 0.0, contactList* list = NULL);
    contactList getContacts(Structure& S, mstreal cdcut = 0.0, contactList* list = NULL);
    contactList getContacts(const vector<Residue*>& residues, mstreal cdcut = 0.0, contactList* list = NULL);
    vector<Residue*> getContactingResidues(Residue* res, mstreal cdcut = 0.0);
  
    /* The function below finds all contacts in the set of provided residues, with constraints applied
     to amino acids at each position. Specifically, for some contact between position i and j, j is
     not constrained (all amino acids allowed) and i is constrained to a single amino acid. We iterate
     over all of the available amino acids (18 since Gly/Pro are excluded) at i, and calculate the CD.
     We also reverse the direction so that all are allowed at i, and a single amino acid is allowed
     at j. This results in up to 18x18 = 324 values per contact. */
    aaConstrainedContactList getConstrainedContacts(const vector<Residue*>& residues, mstreal cdcut = 0.0, aaConstrainedContactList* list = NULL);

    /* Interference is a directional sidechain-to-backbone contact. If A and B
     * are listed as the source and destination residues, respectively, of a
     * contact in the resulting contactList, then some fraction of rotamers at A
     * are clashing with the backbone of B; this fraction is the interference.
     * The first set of functions returns a contact list, where the specified
     * residues act as either intefering residues (i.e., their backbones clash
     * with somebody else's sidechains) or the interfered residues (i.e., their
     * sidechains clash with somebody else's backbone). The second set returns
     * only those contacts, in which the specified residues are being interfered
     * with (and hence, what is effectively being returned are the interfering
     * residues, which is why the functions are named as they are). */
    contactList getInterference(const vector<Residue*>& residues, mstreal incut = 0.0, contactList* list = NULL);
    contactList getInterference(const Structure& S, mstreal incut = 0.0, contactList* list = NULL);
    contactList getInterfering(const vector<Residue*>& residues, mstreal incut = 0.0, contactList* list = NULL);
    contactList getInterfering(const Structure& S, mstreal incut = 0.0, contactList* list = NULL);
  
    /*
     interferenceValue computes the interference of residue A sidechain by residue B backbone.
     By setting aaAllowed, the user can restrict the amino acids that are considered at residue A.
     
     Note: specifying amino acids that are not in the library (e.g. GLY, PRO, non-canonical) will
     result in an error.
     */
    mstreal interferenceValue(Residue* resA, Residue* resB, set<string> aaAllowed = {});

    /* Backbone interaction is a backbone-to-backbone contact, defined as ANY
     * of the backbone atoms (N,Ca,C,O) of two residues being within the cutoff
     * distance. If this criterion is met, the exact distance between the closest
     * pair of backbone atoms from the two sets is reported. Note that by default
     * the residues directly adjacent to the residue of interest are not considered
     * when searching for backbone interactions, this can be adjusted by setting
     * ignoreFlanking. */
    
    mstreal bbInteraction(Residue* resA, Residue* resB);
    contactList getBBInteraction(Residue* res, mstreal dcut = 0.0, int ignoreFlanking = 1, contactList* list = NULL);
    contactList getBBInteraction(Structure& S, mstreal dcut = 0.0, int ignoreFlanking = 1, contactList* list = NULL);
    contactList getBBInteraction(const vector<Residue*>& residues, mstreal dcut = 0.0, int ignoreFlanking = 1, contactList* list = NULL);
    vector<Residue*> getBBInteractingResidues(Residue* res, mstreal dcut = 0.0, int ignoreFlanking = 1);
    
    mstreal getCrowdedness(Residue* res);
    vector<mstreal> getCrowdedness(vector<Residue*>& residues);

    mstreal getFreedom(Residue* res);
    vector<mstreal> getFreedom(vector<Residue*>& residues);
    void clearFreedom() { freedom.clear(); } // useful if one wants to force re-calculation (e.g., with new parameters)
  
    set<string> getAANames() {return aaNames;}

    void openLogFile(string fname, bool append = false);
    void closeLogFile();

  protected:
    mstreal weightOfAvailableRotamers(Residue* res, set<string> available_aa); // computes the total weight of all rotamers available at this position
    mstreal weightOfAvailableAminoAcids(set<string> available_aa); // computes the total weight of all selected amino acids
    void init(const Structure& S);
    void setParams();
    /* given pre-computed collision probabilities, sums up freedom scores. NOTE,
     * does not check whether all the relevant contacting residues have been
     * visited, so must be called only at the right times (that's why protected) */
    mstreal computeFreedom(Residue* res);
    void collProbUpdateOn(Residue* res) { updateCollProb[res] = true; }
    void collProbUpdateOff(Residue* res) { updateCollProb[res] = false; }

  private:
    RotamerLibrary* rotLib;
    bool isRotLibLocal;
    bool strict; //if true, will only consider the amino acid present in the structure in calculations
    AtomPointerVector backbone, ca;
    ProximitySearch *bbNN, *caNN;
    fastmap<Residue*, set<int> > permanentContacts;
    fastmap<Residue*, mstreal> fractionPruned;
    fastmap<Residue*, mstreal> freedom;
    fastmap<Residue*, int> numLibraryRotamers;
    fastmap<Residue*, vector<rotamerID*> > survivingRotamers;
    fastmap<Residue*, fastmap<Residue*, mstreal> > degrees;
    fastmap<Residue*, fastmap<rotamerID*, mstreal> > collProb;
    fastmap<Residue*, fastmap<string, DecoratedProximitySearch<rotamerID*>* > > rotamerHeavySC;
    fastmap<Residue*, fastmap<Residue*, fastmap<string, mstreal> > > interference; // interference[resA][resB][aa] will store how much the backbone of
                                                                 // resB can potentially interfere amino acid aa at resA
    set<string> aaNames;     // amino acids whose rotamers will be considered (all except GLY and PRO)
    mstreal dcut;                  // CA-CA distance cutoff beyond which we do not consider pairwise interactions
    mstreal clashDist, contDist;   // inter-atomic distances for counting main-chain clashes and inter-rotamer contacts, respectively
    fastmap<string, double> aaProp; // amino-acid propensities (in percent)
    bool doNotCountCB;          // if true, CB is not counted as a side-chain atom for counting clashes (except for ALA)
    fstream rotOut;
    /* an internal flag that sets the state of the object with respect to
     * uplading the collision probability mass table. In general, should be
     * false, unless set internally as part of a relevant function (and then
     * unset before returning). */
    fastmap<Residue*, bool> updateCollProb;
    mstreal loCollProbCut, hiCollProbCut; // low and high collision probability cutoffs for computing freedom
    int freedomType;                   // a switch between different formulas for computing freedom
};


#endif
