#ifndef _MSTCONDEG_H
#define _MSTCONDEG_H

#include "msttypes.h"
#include "mstrotlib.h"
#include <set>

using namespace std;
using namespace MST;

class contactList {
  public:
    contactList(Residue* _resi) { resi = _resi; }
    contactList() { resi = NULL; }
    contactList(contactList& other) {
      resi = other.resi;
      resj.assign(other.resj.begin(), other.resj.end());
      degrees.assign(other.degrees.begin(), other.degrees.end());
      infos.assign(other.infos.begin(), other.infos.end());
    }
    void addContact(Residue* _resj, real _degree, string _info = "") {
      resj.push_back(_resj);
      degrees.push_back(_degree);
      infos.push_back(_info);
    }
    int size() { return resj.size(); }
    Residue* source() { return resi; }
    Residue* contResidue(int i) { return resj[i]; }
    real degree(int i) { return degrees[i]; }
    string info(int i) { return infos[i]; }

  private:
    Residue* resi;
    vector<Residue*> resj;
    vector<real> degrees;
    vector<string> infos;
};

class aaRotamers {
  public:
    aaRotamers(string aa, Residue* _pos) { rotamers = new Residue(aa, 1); pos = _pos; }
    ~aaRotamers() { delete rotamers; }
    Residue* position() { return pos; }
    real rotProb(int ri) { return rP[ri]; }
    int rotID(int ri) { return rID[ri]; }
    int atomSize() { return rotamers->atomSize(); }
    Atom& operator[](int i) const { return (*rotamers)[i]; }
    Residue rotamer(int ri) {
      rotamers->makeAlternativeMain(ri);
      Residue ret(*rotamers, false);
      return ret;
    }
    CartesianPoint rotamerAtomCoor(int ri, int ai) {
      return (*rotamers)[ai].getAltCoor(ri);
    }
    string aaName() { return rotamers->getName(); }
    void addRotamer(Residue& res, int _rID, double _rP) {
      // the main coordinate set is some "current" rotamer, so all rotamers are kept in alternative
      // thus, the first time, effectively copy the coordinates twice
      while (1) {
        bool isempty = (rotamers->atomSize() == 0);
        int k = 0;
        for (int i = 0; i < res.atomSize(); i++) {
          Atom& a = res[i];
          if (RotamerLibrary::isHydrogen(a)) continue;
          if (isempty) {
            rotamers->appendAtom(new Atom(a, false)); // the first time into the main set (new atoms)
          } else {
            if ((k >= rotamers->atomSize()) || (!(*rotamers)[k].isNamed(a.getName())))
              MstUtils::error("the new rotamer not consistent with previous ones", "aaRotamers::addRotamer(Residue*)");
            (*rotamers)[k].addAlternative(a.getX(), a.getY(), a.getZ(), 0.0, 1.0); // the into alternatives
            k++;
          }
        }
        if (isempty) { isempty = false; }
        else {
          if (k != rotamers->atomSize()) MstUtils::error("the new rotamer not consistent with previous ones", "aaRotamers::addRotamer(Residue*)");
          break;
        }
      }
      rID.push_back(_rID);
      rP.push_back(_rP);
    }
    int numberOfRotamers() {
      if (rotamers->atomSize() == 0) return 0;
      return (*rotamers)[0].numAlternatives();
    }

  private:
    vector<real> rP;
    vector<int> rID;
    Residue* rotamers;
    Residue* pos;
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

    contactList getContacts(Residue* res, real cdcut = 0.0);
    // vector<contactList> getContacts(vector<Residue*> res, real cdcut = 0.0);
    // vector<Residue*> getContactingResidues(Residue* res, real cdcut = 0.0);
    // vector<Residue*> getContactingResidues(vector<Residue*> res, real cdcut = 0.0);
    // vector<contactList> getContacts(Structure& S, real cdcut = 0.0);
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
    map<Residue*, vector<aaRotamers*> > rotamers;
    map<Residue*, set<int> > permanentContacts;
    map<Residue*, double> fractionPruned;
    map<Residue*, double> freeVolume;
    map<Residue*, int> origNumRots;
    map<Residue*, double> freedom;

    class rotamerAtomInfo {
      public:
        rotamerAtomInfo() { aa = NULL; rotInd = -1; atomInd = -1; } // default constructor needed to be part of vector<>
        rotamerAtomInfo(aaRotamers* _aa, int _rotInd, int _atomInd) { aa = _aa; rotInd = _rotInd; atomInd = _atomInd; }
        rotamerAtomInfo(const rotamerAtomInfo& other) { aa = other.aa; rotInd = other.rotInd; atomInd = other.atomInd; }
        int rotIndex() { return rotInd; }
        int atomIndex() { return atomInd; }
        Residue* position() { return aa->position(); }
        string aaName() { return aa->aaName(); }
        string name() { return (*aa)[atomInd].getName(); }
        CartesianPoint coor() { return (*aa)[atomInd].getAltCoor(rotInd); }
        real rotProb() { return aa->rotProb(rotInd); }
        aaRotamers* getAA() { return aa; }

      private:
        aaRotamers* aa;
        int rotInd;
        int atomInd;
    };

    DecoratedProximitySearch<rotamerAtomInfo> *rotamerHeavySC; // point cloud of rotamer atoms from ALL rotamers

    vector<string> aaNames;     // amino acids whose rotamers will be considered (all except GLY and PRO)
    real dcut;                  // CA-AC distance cutoff beyond which we do not consider pairwise interactions
    real clashDist, contDist;   // inter-atomic distances for counting main-chain clashes and inter-rotamer contacts, respectively
    map<string, double> aaProp; // amino-acid propensities (in percent)
    bool doNotCountCB;          // if true, CB is not counted as a side-chain atom for counting clashes (except for ALA)
};


#endif
