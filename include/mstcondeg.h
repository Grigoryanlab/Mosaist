#ifndef _MSTCONDEG_H
#define _MSTCONDEG_H

#include "msttypes.h"
#include "mstrotlib.h"
#include <set>

using namespace std;
using namespace MST;

template<typename key, typename T>
using fastmap = map<key, T>;

class contactList {
  public:
    contactList() { }
    contactList(const contactList& other) {
      resi = other.resi;
      resj = other.resj;
      degrees = other.degrees;
      infos = other.infos;
      inContact = other.inContact;
      orderedContacts = other.orderedContacts;
    }
    void addContact(Residue* _resi, Residue* _resj, mstreal _degree, string _info = "", bool directional = false) {
      resi.push_back(_resi);
      resj.push_back(_resj);
      degrees.push_back(_degree);
      infos.push_back(_info);
      inContact[_resi][_resj] = resi.size() - 1;
      inContact[_resj][_resi] = resi.size() - 1;
      if ((!directional) && (_resi->getResidueIndex() > _resj->getResidueIndex())) {
        orderedContacts.insert(pair<Residue*, Residue*>(_resj, _resi));
      } else {
        orderedContacts.insert(pair<Residue*, Residue*>(_resi, _resj));
      }
    }
    int size() { return resi.size(); }
    Residue* residueA(int i) { return resi[i]; }
    Residue* residueB(int i) { return resj[i]; }
    Residue* srcResidue(int i) { return resi[i]; }
    Residue* dstResidue(int i) { return resj[i]; }
    vector<Residue*> srcResidues() { return resi; }
    vector<Residue*> destResidues() { return resj; }
    mstreal degree(int i) { return degrees[i]; }
    mstreal degree(Residue* _resi, Residue* _resj);
    string info(int i) { return infos[i]; }
    void sortByDegree(); // sorts the contact list by contact degree, highest to lowest
    vector<pair<Residue*, Residue*> > getOrderedContacts();
    bool areInContact(Residue* A, Residue* B);


  private:
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

    vector<Residue*> resi;
    vector<Residue*> resj;
    vector<mstreal> degrees;
    vector<string> infos;
    fastmap<Residue*, fastmap<Residue*, int> > inContact;
    set<pair<Residue*, Residue*>, contComp> orderedContacts;
};

class ConFind {
  public:
    ConFind(string rotLibFile, Structure& S);
    ConFind(RotamerLibrary* _rotLib, Structure& S);
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

    mstreal contactDegree(Residue* resA, Residue* resB, bool cacheA = true, bool cacheB = true, bool checkNeighbors = true);
    contactList getContacts(Residue* res, mstreal cdcut = 0.0, contactList* list = NULL);
    contactList getContacts(Structure& S, mstreal cdcut = 0.0, contactList* list = NULL);
    contactList getContacts(const vector<Residue*>& residues, mstreal cdcut = 0.0, contactList* list = NULL);
    vector<Residue*> getContactingResidues(Residue* res, mstreal cdcut = 0.0);
    contactList getInterference(vector<Residue*>& residues, mstreal incut = 0.0, contactList* list = NULL);
    contactList getInterference(Structure& S, mstreal incut = 0.0, contactList* list = NULL);

    mstreal getCrowdedness(Residue* res);
    vector<mstreal> getCrowdedness(vector<Residue*>& residues);

    mstreal getFreedom(Residue* res);
    vector<mstreal> getFreedom(vector<Residue*>& residues);
    void clearFreedom() { freedom.clear(); } // useful if one wants to force re-calculation (e.g., with new parameters)

    void openLogFile(string fname, bool append = false);
    void closeLogFile();

  protected:
    mstreal weightOfAvailableRotamers(Residue* res); // computes the total weight of all rotamers available at this position
    void init(Structure& S);
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
    AtomPointerVector backbone, ca;
    ProximitySearch *bbNN, *caNN;
    fastmap<Residue*, set<int> > permanentContacts;
    fastmap<Residue*, mstreal> fractionPruned;
    fastmap<Residue*, mstreal> freedom;
    fastmap<Residue*, int> numLibraryRotamers;
    fastmap<Residue*, vector<rotamerID*> > survivingRotamers;
    fastmap<Residue*, fastmap<Residue*, mstreal> > degrees;
    fastmap<Residue*, fastmap<rotamerID*, mstreal> > collProb;
    fastmap<Residue*, DecoratedProximitySearch<rotamerID*>* > rotamerHeavySC;
    fastmap<Residue*, fastmap<Residue*, mstreal> > interference; // interferance[resA][resB] will store home much the backbone of
                                                              // resB can potentially interfere with the amino-acid choice at resA

    vector<string> aaNames;     // amino acids whose rotamers will be considered (all except GLY and PRO)
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
