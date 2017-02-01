#ifndef _MSTCONDEG_H
#define _MSTCONDEG_H

#include "msttypes.h"
#include "mstrotlib.h"

class ConFind {
  public:
    ConFind(string rotLibFile, Structure& _S);
    ConFind(RotamerLibrary* _rotLib, Structure& S);
    ~ConFind();

    // precomputes all necessary info and data structures for computing on this Structure
    bool cache(Structure& S);
    bool cache(vector<Residue*>& residues);
    bool cache(Residue* res);

    vector<contact> getContacts(Residue* res, real cdcut = 0.01);
    vector<contact> getContacts(vector<Residue*> res, real cdcut = 0.01);
    vector<Residue*> getContactingResidues(Residue* res, real cdcut = 0.01);
    vector<Residue*> getContactingResidues(vector<Residue*> res, real cdcut = 0.01);
    vector<contact> getAllContacts(Structure& S, real cdcut = 0.01);
    vector<real> getFreeVolume(Residue* res, real cdcut = 0.01);
    vector<real> getFreeVolume(vector<Residue*> res, real cdcut = 0.01);
    vector<real> getFreedom(Residue* res, real cdcut = 0.01);
    vector<real> getFreedom(vector<Residue*> res, real cdcut = 0.01);
    // ???? TODO: permanent contacts!!!

  protected:
    init(Structure& _S);

  private:
    RotamerLibrary* rotLib; // TODO: big problem, don't know if should be deleted or not!!!!!
    AtomPointerVector backbone, nonHyd, ca;
    ProximitySearch bbNN, nonHydNN, caNN;
    map<Residue*, vector<aaRotamers*> > rotamers;
    map<Residue*, set<int> > permanentContacts;
    map<Residue*, double> fractionPruned;
    map<Residue*, double> freeVolume;
    map<Residue*, int> origNumRots;
    map<Residue*, double> freedom;
    ProximitySearch rotamerHeavySC, rotamerHeavyBB; // point cloud of rotamer atoms from ALL rotamers

    /* stores which residues have been fully described (their entire environment),
     * as opposed to as just supporting residues. All this really means is that
     * this residue has been cached AND so have all relevant surrounding residues.*/
    map<Residue*, bool> fullyDescribed;

    // amino acids whose rotamers will be considered (all except GLY and PRO)
    static vector<string> aaNames = {"ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "SER", "THR", "TRP", "TYR", "VAL", "ALA"};
    static real dcut = 25.0; // CA-AC distance cutoff beyond which we do not consider pairwise interactions
    static map<string, double> aaProp; // amino-acid propensities (in percent)
}

class aaRotamers {
  public:
    aaRotamers(string aa) { rotamers = new Residue(aa, 1); }
    double rotProb() { return rP; }
    int rotID() { return rID; }
    void addRotamer(Residue* res, int _rID, double _rP) {
      if (rotamers->atomSize() == 0) {
        rotamers->copyAtoms(*res);
      } else {
        int k = 0;
        for (int i = 0; i < res->atomSize(); i++) {
          Atom& a = (*res)[i];
          if (RotamerLibrary::isHydrogen(a)) continue;
          if ((k >= rotamers->atomSize()) || (!(*rotamers)[k].isNamed(a.atomName())))
            MstUtils::error("the new rotamer not consistent with previous ones", "aaRotamers::addRotamer(Residue*)");
          (*rotamers)[k].addAlternative(a.getX(), a.getY(), a.getZ(), 0.0, 1.0);
          k++;
        }
        if (k != rotamers->atomSize()) MstUtils::error("the new rotamer not consistent with previous ones", "aaRotamers::addRotamer(Residue*)");
        rID.push_back(_rID);
        rP.push_back(_rP);
      }
    }
    int numberOfRotamers() {
      if (rotamers->atomSize() == 0) return 0;
      return (*rotamers)[0].numAlternatives() + 1;
    }

  private:
    vector<real> rP;
    vector<int> rID;
    Residue* rotamers;
};

class contact {
  public:
    Residue* resi;
    Residue* resj;
    real degree;
    string info;
    contact(Residue* _resi, Residue* _resj, real _degree, string _info) {
      resi = _resi;
      resj = _resj;
      degree = _degree;
      info = _info;
    }
    contact() { resi = NULL; resj = NULL; degree = 0; info = ""; }
};


#endif
