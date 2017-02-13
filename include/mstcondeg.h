#ifndef _MSTCONDEG_H
#define _MSTCONDEG_H

#include "msttypes.h"
#include "mstrotlib.h"
#include <set>

using namespace std;
using namespace MST;

class cdPartner {

};

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
    }
    int size() { return resi.size(); }
    Residue* residueA(int i) { return resi[i]; }
    Residue* residueB(int i) { return resj[i]; }
    real degree(int i) { return degrees[i]; }
    real degree(Residue* _resi, Residue* _resj);
    string info(int i) { return infos[i]; }

  private:
    vector<Residue*> resi;
    vector<Residue*> resj;
    vector<real> degrees;
    vector<string> infos;
    map<Residue*, map<Residue*, int> > inContact;
};

class ConFind {
  public:
    ConFind(string rotLibFile, Structure& S);
    ConFind(RotamerLibrary* _rotLib, Structure& S);
    ~ConFind();

    // precomputes all necessary info and data structures for computing on this Structure
    bool cache(Structure& S, fstream* rotOut = NULL);
    bool cache(vector<Residue*>& residues, fstream* rotOut = NULL);
    bool cache(Residue* res, fstream* rotOut = NULL);

    // find those residues that are close enough to affect the passed residue(s)
    vector<Residue*> getNeighbors(Residue* residue);
    vector<Residue*> getNeighbors(vector<Residue*>& residues);

    /* this function encodes whether a given atom counts as "side-chain" for the
     * purposes of finding sidechain-to-sidechain contacts. */
    bool countsAsSidechain(Atom& a);

    real contactDegree(Residue* resA, Residue* resB, bool doNotCache = false);
    contactList getContacts(Residue* res, real cdcut = 0.0);
    vector<Residue*> getContactingResidues(Residue* res, real cdcut = 0.0);
    contactList getContacts(Structure& S, real cdcut = 0.0);
    // vector<real> getFreeVolume(Residue* res, real cdcut = 0.0);
    // vector<real> getFreeVolume(vector<Residue*> res, real cdcut = 0.0);
    // vector<real> getFreedom(Residue* res, real cdcut = 0.0);
    // vector<real> getFreedom(vector<Residue*> res, real cdcut = 0.0);
    // ???? TODO: permanent contacts!!!

  protected:
    real weightOfAvailableRotamers(Residue* res); // computes the total weight of all rotamers available at this position
    void init(Structure& S);
    void setParams();

  private:
    RotamerLibrary* rotLib;
    bool isRotLibLocal;
    AtomPointerVector backbone, nonHyd, ca;
    ProximitySearch *bbNN, *nonHydNN, *caNN;
    map<Residue*, set<int> > permanentContacts;
    map<Residue*, double> fractionPruned;
    map<Residue*, double> freeVolume;
    map<Residue*, int> origNumRots;
    map<Residue*, double> freedom;
    map<Residue*, vector<rotamerID*> > survivingRotamers;
    map<Residue*, DecoratedProximitySearch<rotamerID*>* > rotamerHeavySC;

//    DecoratedProximitySearch<rotamerAtomInfo> *rotamerHeavySC; // point cloud of rotamer atoms from ALL rotamers

    vector<string> aaNames;     // amino acids whose rotamers will be considered (all except GLY and PRO)
    real dcut;                  // CA-AC distance cutoff beyond which we do not consider pairwise interactions
    real clashDist, contDist;   // inter-atomic distances for counting main-chain clashes and inter-rotamer contacts, respectively
    map<string, double> aaProp; // amino-acid propensities (in percent)
    bool doNotCountCB;          // if true, CB is not counted as a side-chain atom for counting clashes (except for ALA)
};


#endif
