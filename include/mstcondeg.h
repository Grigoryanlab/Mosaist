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
    }
    void addContact(Residue* _resi, Residue* _resj, real _degree, string _info = "") {
      resi.push_back(_resi);
      resj.push_back(_resj);
      degrees.push_back(_degree);
      infos.push_back(_info);
      inContact[_resi][_resj] = resi.size() - 1;
      inContact[_resj][_resi] = resi.size() - 1;
      orderedContacts.insert(pair<Residue*, Residue*>(_resi, _resj));
    }
    int size() { return resi.size(); }
    Residue* residueA(int i) { return resi[i]; }
    Residue* residueB(int i) { return resj[i]; }
    real degree(int i) { return degrees[i]; }
    real degree(Residue* _resi, Residue* _resj);
    string info(int i) { return infos[i]; }
    vector<pair<Residue*, Residue*> > getOrderedContacts();


  private:
    struct contComp {
      bool operator() (const pair<Residue*, Residue*>& lhs, const pair<Residue*, Residue*>& rhs) const {
        int lhsI = lhs.first->getResidueIndex();
        int rhsI = rhs.first->getResidueIndex();
        if (lhsI == rhsI) {
          lhsI = lhs.second->Residue::getResidueIndex();
          rhsI = rhs.second->Residue::getResidueIndex();
        }
        return lhsI < rhsI;
      }
    };

    vector<Residue*> resi;
    vector<Residue*> resj;
    vector<real> degrees;
    vector<string> infos;
    fastmap<Residue*, fastmap<Residue*, int> > inContact;
    set<pair<Residue*, Residue*>, contComp> orderedContacts;
};

class ConFind {
  public:
    ConFind(string rotLibFile, Structure& S);
    ConFind(RotamerLibrary* _rotLib, Structure& S);
    ~ConFind();
    void setFreedomParams(real _loCollProbCut, real _hiCollProbCut) { loCollProbCut = _loCollProbCut; hiCollProbCut = _hiCollProbCut; }

    // precomputes all necessary info and data structures for computing on this Structure
    void cache(Structure& S);
    void cache(vector<Residue*>& residues);
    void cache(Residue* res);

    // find those residues that are close enough to affect the passed residue(s)
    vector<Residue*> getNeighbors(Residue* residue);
    vector<Residue*> getNeighbors(vector<Residue*>& residues);
    bool areNeighbors(Residue* resA, Residue* resB);

    /* this function encodes whether a given atom counts as "side-chain" for the
     * purposes of finding sidechain-to-sidechain contacts. */
    bool countsAsSidechain(Atom& a);

    real contactDegree(Residue* resA, Residue* resB, bool cacheA = true, bool cacheB = true, bool checkNeighbors = true);
    contactList getContacts(Residue* res, real cdcut = 0.0, contactList* list = NULL);
    contactList getContacts(Structure& S, real cdcut = 0.0, contactList* list = NULL);
    contactList getContacts(vector<Residue*>& residues, real cdcut = 0.0, contactList* list = NULL);
    vector<Residue*> getContactingResidues(Residue* res, real cdcut = 0.0);

    real getCrowdedness(Residue* res);
    vector<real> getCrowdedness(vector<Residue*>& residues);

    real getFreedom(Residue* res);
    vector<real> getFreedom(vector<Residue*>& residues);
    void clearFreedom() { freedom.clear(); } // useful if one wants to force re-calculation (e.g., with new parameters)

    void openLogFile(string fname, bool append = false);
    void closeLogFile();

  protected:
    real weightOfAvailableRotamers(Residue* res); // computes the total weight of all rotamers available at this position
    void init(Structure& S);
    void setParams();
    /* given pre-computed collision probabilities, sums up freedom scores. NOTE,
     * does not check whether all the relevant contacting residues have been
     * visited, so must be called only at the right times (that's why protected) */
    real computeFreedom(Residue* res);
    void collProbUpdateOn(Residue* res) { updateCollProb[res] = true; }
    void collProbUpdateOff(Residue* res) { updateCollProb[res] = false; }

  private:
    RotamerLibrary* rotLib;
    bool isRotLibLocal;
    AtomPointerVector backbone, ca;
    ProximitySearch *bbNN, *caNN;
    fastmap<Residue*, set<int> > permanentContacts;
    fastmap<Residue*, real> fractionPruned;
    fastmap<Residue*, real> freedom;
    fastmap<Residue*, int> numLibraryRotamers;
    fastmap<Residue*, vector<rotamerID*> > survivingRotamers;
    fastmap<Residue*, fastmap<Residue*, real> > degrees;
    fastmap<Residue*, fastmap<rotamerID*, real> > collProb;
    fastmap<Residue*, DecoratedProximitySearch<rotamerID*>* > rotamerHeavySC;

    vector<string> aaNames;     // amino acids whose rotamers will be considered (all except GLY and PRO)
    real dcut;                  // CA-AC distance cutoff beyond which we do not consider pairwise interactions
    real clashDist, contDist;   // inter-atomic distances for counting main-chain clashes and inter-rotamer contacts, respectively
    fastmap<string, double> aaProp; // amino-acid propensities (in percent)
    bool doNotCountCB;          // if true, CB is not counted as a side-chain atom for counting clashes (except for ALA)
    fstream rotOut;
    /* an internal flag that sets the state of the object with respect to
     * uplading the collision probability mass table. In general, should be
     * false, unless set internally as part of a relevant function (and then
     * unset before returning). */
    fastmap<Residue*, bool> updateCollProb;
    real loCollProbCut, hiCollProbCut; // low and high collision probability cutoffs for computing freedom
};


#endif
