#ifndef _MSTCONDEG_H
#define _MSTCONDEG_H

#include "msttypes.h"
#include "mstrotlib.h"

class ConFind {
  public:
    ConFind(string rotLibFile);
    ConFind(string rotLibFile, Structure& S);
    ConFind(RotamerLibrary* _rotLib);
    ConFind(RotamerLibrary* _rotLib, Structure& S);
    ~ConFind(); // destuct all ProximitySearch objects

    // precomputes all necessary info and data structures for computing on this Structure
    void cache(Structure& S, vector<Residue*>* residues = NULL);

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

  private:
    RotamerLibrary* rotLib;
    map<Residue*, vector<rotamer> > rotamers;
    map<Residue*, set<int> > permanentContacts;
    map<Residue*, double> fractionPruned;
    map<Residue*, double> freeVolume;
    map<Residue*, int> origNumRots;
}

class rotamer {
  public:
    rotamer() { atomsSC = atomsBB = NULL; rP = aaP = 1.0; aaN = "XXX"; rID = -1; }
    rotamer(ProximitySearch* _atomsSC, ProximitySearch* _atomsBB, double _aaProp, double _rotProb, string _name, int _rotID) {
      atomsSC = _atomsSC; atomsBB = _atomsBB; rP = _rotProb; aaP = _aaProp; aaN = _name; rID = _rotID;
    }
    ProximitySearch* gridSC() { return atomsSC; }
    ProximitySearch* gridBB() { return atomsBB; }
    double aaProp() { return aaP; }
    double rotProb() { return rP; }
    string aaName() { return aaN; }
    int rotID() { return rID; }

  private:
    ProximitySearch *atomsSC, *atomsBB;
    double rP, aaP;
    string aaN;
    int rID;
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
