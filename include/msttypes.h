#ifndef _MSTTYPES_H
#define _MSTTYPES_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <locale>
#include <stdio.h>
#include <string.h>
#include <iomanip>
#include <vector>
#include <map>
#include <set>
#include <limits>
#include <algorithm>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <chrono>
#include <execinfo.h>
#include <signal.h>
#include <random>
#undef assert

using namespace std;

namespace MST {

// forward declarations
class Chain;
class Residue;
class Atom;
class Structure;
class CartesianPoint;

typedef double mstreal;
typedef Structure System;                // for interchangability with MSL

class Structure {
  friend class Chain;

  public:
    Structure();
    Structure(string pdbFile, string options = "");
    Structure(istream& is, string options = "");
    Structure(const Structure& S);
    Structure(Chain& C);
    Structure(Residue& R);
    Structure(const vector<Atom*>& atoms);
    Structure(const vector<Residue*>& residues);
    ~Structure();

    /* Produces a hybrid structure by coopying atoms from first input structure
     * into residues of the second one, preserving the chain topology of the latter. */
    static Structure combine(const Structure& atomsStruct, const Structure& topoStruct, bool renameResidues = true);

    void readPDB(const string& pdbFile, string options = "");
    void readPDB(istream& is, string options = "");
    void writePDB(const string& pdbFile, string options = "") const;
    void writePDB(ostream& ofs, string options = "") const;
    void writeData(const string& dataFile) const;
    void writeData(ostream& ofs) const;
    void readData(const string& dataFile);
    void readData(istream& ifs);
    void reset();
    Structure& operator=(const Structure& A);
    int chainSize() const { return chains.size(); }
    int residueSize() const { return numResidues; }
    int positionSize() const { return residueSize(); }  // for interchangability with MSL
    int atomSize() const { return numAtoms; }
    Chain* getChainByID(string id) { return (chainsByID.find(id) != chainsByID.end()) ? chainsByID[id] : NULL; }
    Chain* getChainBySegID(string id) { return (chainsBySegID.find(id) != chainsBySegID.end()) ? chainsBySegID[id] : NULL; }
    Chain& getChain(int i) const { return (*this)[i]; }
    Residue& getResidue(int i) const;
    Chain& operator[](int i) const { return *(chains[i]); }
    vector<Atom*> getAtoms() const;
    vector<Residue*> getResidues() const;
    void setName(const string& _name) { name = _name; }
    string getName() const { return name; }
    void renumber(int startResNum=1, int startAtomIndex=1); // make residue numbering consequitive in each chain and atom index consequitive throughout
    // looks at the length of the peptide bond between adjacent residues to figure out where chains break
    void reassignChainsByConnectivity(Structure& dest, mstreal maxPeptideBond = 2.0);
    Structure reassignChainsByConnectivity(mstreal maxPeptideBond = 2.0);

    /* ----- functions that grow/shrink structure ----- */
    /* returns false if the chain name collides with existing chain names and no suitable single-letter
     * chain name was found as replacement OR if renaiming was not allowed. This could still mean that
     * a multi-character name is picked that is unique, but that's not technically correct for output,
     * so false will still be returned. Note that in the latter case, the segment ID will be renamed as
     * well to be the same multi-character name. So, although in the output chain names will repeat,
     * segment names will still be unique. If it fails to find an even multi-character name, errors. */
    bool appendChain(Chain* C, bool allowRename = true);
    Chain* appendChain(string cid, bool allowRename = true);
    void deleteChain(Chain* chain);
    void deleteShortChains(int sizeThreshold = 5);

    /* makes a copy of the atom, then decides where it is supposed to go and inserts it
     * into the Structure, creating a new Chain and/or Residue as needed. */
    void addAtom(Atom* A);
    void addAtom(Atom& A) { addAtom(&A); }
    void addAtoms(vector<Atom*> atoms) { addAtoms(&atoms); }
    void addAtoms(vector<Atom*>* atoms);

    /* makes a copy of the residue, then decides where it is supposed to go and
     * inserts it into the Structure, creating a new Chain as needed. */
    Residue* addResidue(Residue* res);

    /* ----- functions that grow/shrink structure ----- */

    int getResidueIndex(Residue* res);

    /* == and != operators are needed to convert vector<Structure> into python
     * lists via boost.python. This is because python lists are quite a bit more
     * powerful than C++ vectors, enabling, for example, contains queries. */
    bool operator==(const Structure& other) { return (this == &other); }
    bool operator!=(const Structure& other) { return (this != &other); }

  protected:
    void incrementNumAtoms(int delta = 1) { numAtoms += delta; }
    void incrementNumResidues(int delta = 1) { numResidues += delta; }
    void deletePointers();
    void copy(const Structure& S);

  private:
    vector<Chain*> chains;
    string name;
    int numResidues, numAtoms;
    // NOTE: thse two maps are maintained for convenience and will not guarantee the lack of collisions. That is,
    // if more than one chain use the same ID or segment ID, these maps will only store the last one added.
    map<string, Chain*> chainsByID;
    map<string, Chain*> chainsBySegID;
};

class Chain {
  friend class Residue;
  friend class Structure;

  public:
    Chain();
    Chain(const Chain& C);
    Chain(const string& chainID, const string& segID);
    Chain(const string& chainID, const string& segID, const vector<Residue*>& residues);
    ~Chain();

    int residueSize() const { return residues.size(); }
    int positionSize() const { return residueSize(); }  // for interchangability with MSL
    int atomSize() const { return numAtoms; }
    Residue& operator[](int i) const { return *(residues[i]); }
    Residue& getResidue(int i) const { return (*this)[i]; }
    vector<Residue*> getResidues() { return residues; }
    vector<Atom*> getAtoms();
    string getID() const { return cid; }
    string getSegID() const { return sid; }
    Structure* getParent() const { return parent; }
    Structure* getStructure() const { return getParent(); }
    int getResidueIndex(const Residue* res); // this in the in-chain residue index!
    int getIndex() const; // index of this chain within the larger structure

    /* convenience functoins, not efficient (linear search). If you need to do this a lot,
     * call getResidues() and construct your own data structure (e.g., a map<>) for fast lookups. */
    Residue* findResidue(string resname, int resnum);
    Residue* findResidue(string resname, int resnum, char icode);

    void setID(string _cid) { cid = _cid; }
    void setSegID(string _sid) { sid = _sid; }

    /* ----- functions that grow/shrink structure ----- */
    void appendResidue(Residue* R);
    void insertResidue(Residue* R, int index); // insert the Residue in such a way that it ends up being at index i
    void appendResidueCopies(const vector<Residue*>& residues);
    Residue* insertResidueCopy(Residue* R, int index = -1); // same, but copies the residue first
    Residue* insertResidueCopy(Residue& R, int index = -1); // same, but copies the residue first
    /* ----- functions that grow/shrink structure ----- */

  protected:
    void setParent(Structure* p) { parent = p; } // will not itself update residue/atom counts in parent
    void incrementNumAtoms(int delta = 1);

  private:
    vector<Residue*> residues;
    map<Residue*, int> residueIndexInChain; // to enable quick look-ups of up/down-stream residues
    Structure* parent;
    int numAtoms;
    string cid, sid;
};

class Residue {
  friend class Structure;
  friend class Chain;
  friend class Atom;

  public:
    Residue();
    Residue(const Residue& R, bool copyAlt = true);
    Residue(string _resname, int _resnum, char _icode = ' ');
    ~Residue();

    int atomSize() const { return atoms.size(); }
    vector<Atom*> getAtoms() { return atoms; }
    Atom& operator[](int i) const { return *(atoms[i]); }
    Atom& getAtom(int i) const { return *(atoms[i]); }
    Chain* getChain() const { return parent; }
    string getChainID(bool strict = true);
    string getName() const { return resname; }
    int getNum() const { return resnum; }
    char getIcode() const { return icode; }
    bool isNamed(const string& _name) const { return (resname.compare(_name) == 0); }
    bool isNamed(const char* _name) const { return (strcmp(resname.c_str(), _name) == 0); }
    Atom* findAtom(string _name, bool strict = true) const; // returns NULL if not found and if strict is false
    bool atomExists(string _name) { return (findAtom(_name, false) != NULL); } // mostly for interchangeability with MSL, better to use findAtom and check for NULL
    Chain* getParent() const { return parent; }
    Structure* getStructure() const { return (parent == NULL) ? NULL : parent->getParent(); }

    void setName(const char* _name) { resname = (string) _name; }
    void setName(const string& _name) { resname = _name; }
    void setIcode(char _icode) { icode = _icode; }
    void setNum(int num) { resnum = num; }
    void copyAtoms(Residue& R, bool copyAlt = true);
    void copyAtoms(const vector<Atom*>& _atoms, bool copyAlt = true);
    void compactify() { atoms.shrink_to_fit(); } // saves memory by adjusting capacity to match length

    /* for all atoms in this residue, overwrite the main coordinate set with the
     * coordinate set from the alternative with the specified index. */
    void makeAlternativeMain(int altInd);

    /* ----- functions that grow/shrink structure ----- */
    void appendAtom(Atom* A);
    void appendAtoms(const vector<Atom*>& A);
    void deleteAtoms();
    void deleteAtom(int ind);

    /* replaces the residue's atom vector with the given vector of atoms. By default,
     * all old atoms are deleted (i.e., removed from the residue's atom vector and
     * destructed). However, if the second argument is passed, will only delete the
     * atoms that were at the specified indices in the old atom array. Note that
     * this function is flexible enough to do things like erase a set of one or more
     * atoms, insert a set of one or more atoms, both, replace the entire
     * set of atoms with a new set, destroying the old ones, or even simply change
     * the permutation of the existing atoms. The order of atoms in the new vector
     * will be: any old ones that survived, in their initial order, followed by any
     * newly added atoms, in the specified order. */
    void replaceAtoms(const vector<Atom*>& newAtoms, vector<int>* oldAtoms = NULL);
    // same as above, but old Atoms are identified via Atom pointers
    void replaceAtoms(const vector<Atom*>& newAtoms, const vector<Atom*>& oldAtoms);
    /* ----- end functions that grow/shrink structure -- */

    Residue* previousResidue();
    Residue* nextResidue();
    Residue* iPlusDelta(int off);
    mstreal getPhi(bool strict = true);
    mstreal getPsi(bool strict = true);
    mstreal getOmega(bool strict = true);
    static bool isBadDihedral(mstreal dihe) { return dihe == Residue::badDihedral; }
    static bool areBonded(const Residue& resN, const Residue& resC, mstreal maxPeptideBond = 2.0);
    static bool areBonded(Residue* resN, Residue* resC, mstreal maxPeptideBond = 2.0) { return Residue::areBonded(*resN, *resC, maxPeptideBond); }

    int getResidueIndex() const;
    int getResidueIndexInChain() const;

    static const mstreal badDihedral; // value that signals a dihedral angle that could not be computed for some reason

    friend ostream & operator<<(ostream &_os, const Residue& _res) {
      if (_res.getParent() != NULL) {
        _os << _res.getParent()->getID() << ",";
      }
      _os << _res.getNum() << " " << _res.getName();
      return _os;
    }
    // so that sets of residues, for example, are sorted in the right order
    friend bool operator<(const Residue& r1, const Residue& r2) {
      return r1.getResidueIndex() < r2.getResidueIndex();
    }

  protected:
    void setParent(Chain* _parent) { parent = _parent; } // will not itself update residue/atom counts in parents

  private:
    vector<Atom*> atoms;
    Chain* parent;
    int resnum;
    string resname;
    char icode;
};

class Atom {
  friend class Structure;
  friend class Chain;
  friend class Residue;
  friend class AtomPointerVector;
  friend class AtomContainer;

  public:
    Atom();
    Atom(const Atom& A, bool copyAlt = true);
    Atom(const Atom* A, bool copyAlt = true) : Atom(*A, copyAlt) {}
    Atom(int _index, const string& _name, mstreal _x, mstreal _y, mstreal _z, mstreal _B, mstreal _occ, bool _het, char _alt = ' ', Residue* _parent = NULL);
    ~Atom();

    mstreal getX() const { return x; }
    mstreal getY() const{ return y; }
    mstreal getZ() const{ return z; }
    bool hasInfo() const { return info != NULL; }
    mstreal& operator[](int i);
    const mstreal& operator[](int i) const;
    CartesianPoint getCoor() const;
    CartesianPoint getAltCoor(int altInd) const;
    mstreal getAltB(int altInd) const { return info->getAltB(altInd); }
    mstreal getAltOcc(int altInd) const { return info->getAltOcc(altInd); }
    char getAltLocID(int altInd) const  { return info->getAltLocID(altInd); }
    mstreal getB() const { return info->B; }
    mstreal getOcc() const { return info->occ; }
    string getName() const { return string(info->name); }
    char* getNameC() { return info->name; }
    bool isHetero() const { return info->het; }
    int getIndex() const { return info->index; }
    char getAlt() const { return info->alt; }
    void setAlt(char a) const { info->alt = a; }
    bool isNamed(const char* _name) const { return (strcmp(info->name, _name) == 0); }
    bool isNamed(const string& _name) const { return isNamed(_name.c_str()); }
    bool hasAlternatives() const { return (info->alternatives != NULL); }
    int numAlternatives() const { return (info->alternatives == NULL) ? 0 : info->alternatives->size(); }
    Residue* getParent() const { return info->parent; }
    Residue* getResidue() const { return info->parent; }
    Chain* getChain() const { return (info->parent == NULL) ? NULL : info->parent->getParent(); }
    Structure* getStructure() { Chain* chain = getChain(); return (chain == NULL) ? NULL : chain->getParent(); }
    mstreal getMass() const { return Atom::getMass(info->name); }
    static mstreal getMass(const char* name);

    void setName(const char* _name) { info->setName(_name); }
    void setName(const string& _name) { info->setName(_name); }
    void setX(mstreal _x) { x = _x; }
    void setY(mstreal _y) { y = _y; }
    void setZ(mstreal _z) { z = _z; }
    void setCoor(mstreal _x, mstreal _y, mstreal _z) { x = _x; y = _y; z = _z; }
    void setCoor(const CartesianPoint& xyz);
    void setCoor(const Atom& a);
    void setCoor(const Atom* a) { setCoor(*a); }
    void setAltCoor(int ai, mstreal _x, mstreal _y, mstreal _z) { info->setAltCoor(ai, _x, _y, _z); }
    void setOcc(mstreal _occ) { info->occ = _occ; }
    void setB(mstreal _B) { info->B = _B; }
    void seetHet(bool _het) { info->het = _het; }
    void setIndex(int _index) { info->index = _index; }

    /* make the alternative with the specified index the main one, making the current
     * main position the alternative with the specified index. Calling this twice with
     * the same index will return things back to the way they were originally. */
    void swapWithAlternative(int altInd);

    /* overwrite the main coordinate set with the coordinate set from the alternative
     * with the specified index. */
    void makeAlternativeMain(int altInd);

    void addAlternative(mstreal _x, mstreal _y, mstreal _z, mstreal _B, mstreal _occ, char _alt = ' ') { info->addAlternative(_x, _y, _z, _B, _occ, _alt); }
    void addAlternative(const Atom& a) { addAlternative(a.getX(), a.getY(), a.getZ(), a.getB(), a.getOcc()); }
    void removeLastAlternative() { info->removeLastAlternative(); }
    void removeAlternative(int i) { info->removeAlternative(i); }
    void clearAlternatives() { info->clearAlternatives(); }

    string pdbLine() { return pdbLine((info->parent == NULL) ? 1 : info->parent->getNum(), getIndex()); }
    string pdbLine(int resIndex, int atomIndex);

    mstreal distance(const Atom& another) const;
    mstreal distance(const Atom* another) const { return distance(*another); }
    mstreal distance2(const Atom& another) const;
    mstreal distance2(const Atom* another) const { return distance2(*another); }
    mstreal angle(const Atom& A, const Atom& B, bool radians = false) const;
    mstreal angle(const Atom* A, const Atom* B, bool radians = false) const;
    mstreal dihedral(const Atom& A, const Atom& B, const Atom& C, bool radians = false) const;
    mstreal dihedral(const Atom* A, const Atom* B, const Atom* C, bool radians = false) const;

    /* Sets the coordinates of the atom based on internal coordinates relative
     * to three other atoms: thA, anA, diA, A (this atom). The internal
     * coordinates are distance diA-A (di), angle anA-diA-A (an), and dihedral
     * angle thA - anA - diA - A (th). */
    bool build(const Atom& diA, const Atom& anA, const Atom& thA, mstreal di, mstreal an, mstreal th, bool radians = false);
    bool build(const Atom* diA, const Atom* anA, const Atom* thA, mstreal di, mstreal an, mstreal th, bool radians = false) {
      return build(*diA, *anA, *thA, di, an, th, radians);
    }

    void write(ostream& _os) const; // write Atom to a binary stream
    void read(istream& _is);  // read Atom from a binary stream
    friend ostream & operator<<(ostream &_os, const Atom& _atom);

    /* Strips the Atom of all of its additional information, besides the 3D coor-
     * dinate. This can be useful when memory needs to be preserved and all the
     * other info is not necessary. NOTE: after calling this function, attempts
     * to access any properties other than coordinates will lead to undefined
     * behavior, including segmentation faults. For efficiency, checking for
     * whether the information exists or has been stripepd is not performed. */
    void stripInfo();

  protected:
    void setParent(Residue* _parent) { info->parent = _parent; } // will not itself update residue/atom counts in parents

  private:
    mstreal x, y, z;
    class atomInfo {
      public:
        // data structure for storing information about alternative atom locations
        class altInfo {
          public:
            altInfo() { x = y = z = occ = B = 0; alt = ' '; }
            altInfo(const altInfo& A) { x = A.x; y = A.y; z = A.z; B = A.B; occ = A.occ; alt = A.alt; }
            altInfo(mstreal _x, mstreal _y, mstreal _z, mstreal _occ, mstreal _B, char _alt) { x = _x; y = _y; z = _z; B = _B; occ = _occ; alt = _alt; }
            mstreal x, y, z, occ, B;
            char alt;
        };

        atomInfo();
        atomInfo(const atomInfo& other, bool copyAlt = true);
        atomInfo(int _index, const string& _name, mstreal _B, mstreal _occ, bool _het, char _alt = ' ', Residue* _parent = NULL);
        ~atomInfo();

        void setName(const char* _name);
        void setName(const string& _name) { setName(_name.c_str()); }
        CartesianPoint getAltCoor(int altInd) const;
        mstreal getAltB(int altInd) const;
        mstreal getAltOcc(int altInd) const;
        char getAltLocID(int altInd) const;
        void setAltCoor(int ai, mstreal _x, mstreal _y, mstreal _z);
        void addAlternative(mstreal _x, mstreal _y, mstreal _z, mstreal _B, mstreal _occ, char _alt);
        void removeLastAlternative();
        void removeAlternative(int i);
        void clearAlternatives();

        mstreal occ, B;
        char *name, alt;
        Residue* parent;
        bool het;
        int index;
        vector<altInfo>* alternatives; /* since this is a pointer, and will be NULL for most atoms, it's fine
                                        * to use vector here in terms of memory, but very convenient for use */
    };
    atomInfo* info;
};
ostream & operator<<(ostream &_os, const Atom& _atom); // this just to silence a silly compiler warning

/* The following several classes look and feel like MSL classes, BUT (importantly) their
 * use is absolutely optional, and none of the basic MST datastructures use them. On
 * the other hand, they can be created from those basic types for a similar use as
 * in MSL when needed (and only when needed). For example, CartesianPoint knows how
 * to construct itself from atom, both via a constructor and assignment operator.
 * Similarly, AtomPointerVector knows how to construct itself from vector<Atom*>. */
class CartesianPoint : public vector<mstreal> {
  /* this class it no limited to 3D vectors, though some of the functions will only
   * work with 3D vectors. The intention is to make it general, such that if a 3D
   * vector is required (or another dimension), appropriate assertions are made. */
  public:
    // inherit a bunch of useful constructors from vector
    CartesianPoint() : vector<mstreal>() { }
    CartesianPoint(size_t sz) : vector<mstreal>(sz) { }
    CartesianPoint(size_t sz, mstreal val) : vector<mstreal>(sz, val) { }
    CartesianPoint(const CartesianPoint& other) : vector<mstreal>(other) { }
    CartesianPoint(const vector<mstreal>& other) : vector<mstreal>(other) { }
    CartesianPoint(mstreal x, mstreal y, mstreal z) : vector<mstreal>(3, 0) { (*this)[0] = x; (*this)[1] = y; (*this)[2] = z; }

    CartesianPoint(const Atom& A);
    CartesianPoint(const Atom* A) : CartesianPoint(*A) {}
    CartesianPoint& operator+=(const CartesianPoint& rhs);
    CartesianPoint& operator-=(const CartesianPoint& rhs);
    CartesianPoint& operator/=(const mstreal& s);
    CartesianPoint& operator*=(const mstreal& s);
    CartesianPoint& operator+=(const mstreal& d);
    const CartesianPoint operator+(const CartesianPoint &other) const;
    const CartesianPoint operator-(const CartesianPoint &other) const;
    const CartesianPoint operator*(const mstreal& s) const;
    const CartesianPoint operator+(const mstreal& d) const;
    const CartesianPoint operator/(const mstreal& s) const;
    const CartesianPoint operator-() const;
    // CartesianPoint& operator=(const Atom& A);
    const double operator*(const CartesianPoint& other) const { return this->dot(other); }

    mstreal norm() const;
    mstreal norm2() const;
    mstreal mean() const;
    mstreal stdev() const;
    mstreal var() const;
    mstreal median() const;
    mstreal sum() const;
    CartesianPoint cross(const CartesianPoint& other) const;
    mstreal dot(const CartesianPoint& other) const;
    CartesianPoint getUnit() const { double L = norm(); return (*this/L); };
    CartesianPoint elemProd(const CartesianPoint& other) const; // element-wise product

    // a few special access operations
    mstreal getX() const { return (*this)[0]; }
    mstreal getY() const { return (*this)[1]; }
    mstreal getZ() const { return (*this)[2]; }

    void convertToSphericalCoordinates(mstreal& radius, mstreal& polarAngle, mstreal& azimuthalAngle); // physics (ISO 80000-2:2019 convention), returns radians
    void setPositionBySphericalCoordinates(mstreal radius, mstreal polarAngle, mstreal azimuthalAngle); // physics (ISO 80000-2:2019 convention), angles in radians

    mstreal distance(const CartesianPoint& another) const;
    mstreal distance(const CartesianPoint* another) const { return distance(*another); }
    mstreal distancenc(const CartesianPoint& another) const; // no size check (for speed)
    mstreal distancenc(const CartesianPoint* another) const { return distancenc(*another); } // no size check (for speed)
    mstreal distance2(const CartesianPoint& another) const;
    mstreal distance2(const CartesianPoint* another) const { return distance2(*another); }
    mstreal distance2nc(const CartesianPoint& another) const; // no size check (for speed)
    mstreal distance2nc(const CartesianPoint* another) const { return distance2nc(*another); } // no size check (for speed)
    mstreal cosineAngle(const CartesianPoint& another) const;

    friend ostream & operator<<(ostream &_os, const CartesianPoint& _p) {
      for (int i = 0; i < _p.size(); i++) {
        _os << _p[i];
        if (i != _p.size() - 1) _os << " ";
      }
      return _os;
    }
};

class CartesianGeometry {
  public:
    static mstreal dihedral(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4, bool radians = false);
    static mstreal dihedral(const CartesianPoint * _p1, const CartesianPoint * _p2, const CartesianPoint * _p3, const CartesianPoint * _p4, bool radians = false);
    static mstreal angle(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, bool radians = false);
    static mstreal angle(const CartesianPoint * _p1, const CartesianPoint * _p2, const CartesianPoint * _p3, bool radians = false);
    static mstreal angleDiff(mstreal A, mstreal B, bool radians = false);
    static mstreal angleDiffCCW(mstreal A, mstreal B, bool radians = false);

    /* Given a list of angles defined on the full circle, it is actually not clear
     * where the min and max are (since angles are defined on a circle). Exampe:
     * we are given -179 and +179. Is the range essentially the entire circle,
     * from -179 to +179? Or is it just a two-degree arc from 179 to 181 (aka
     * -179)? Both are correct. This function finds the min/max bounds that
     * minimize the arc length corresponding to the range. This is often
     * convenient. The functoin will try each angle value as the min, find its
     * corresponding max, and pick the pair that minimizes the arc length. */
    static pair<mstreal, mstreal> angleRange(const vector<mstreal>& angles, bool radians = false);

    /* Just like in the above, the average of a list of angles is also not
     * necessarily uniquely defined. This functon defines it by first finding
     * the angular range with the smallest arc length, and then finding the mean
     * with respect to that range (i.e., via the average deviation from the
     * minimal value). */
    static mstreal angleMean(const vector<mstreal>& angles, bool radians = false);

    /* gradient vector must be of length 6, and will be filled with partial
     * derivatives d(distance)/dc, where c runs over x, y, and z of the first
     * atom, then second atom. */
    template <class T>
    static mstreal distance(const CartesianPoint& atom1, const CartesianPoint& atom2, T& grad);

    /* gradient vector must be of length 9, and will be filled with partial
     * derivatives d(distance)/dc, where c runs over x, y, and z of the first
     * atom, then second atom, then third atom. */
    template <class T>
    static mstreal angle(const CartesianPoint& atom1, const CartesianPoint& atom2, const CartesianPoint& atom3, T& grad, bool radians = false);

    /* gradient vector must be of length 12, and will be filled with partial
     * derivatives d(distance)/dc, where c runs over x, y, and z of the first
     * atom, then second atom, then third atom, then fourth atom. */
    template <class T>
    static mstreal dihedral(const CartesianPoint& atom1, const CartesianPoint& atom2, const CartesianPoint& atom3, const CartesianPoint& atom4, T& grad, bool radians = false);

    /* Tests implementation of analytical gradients of bond, angle, and dihedral
     * using finite difference for comparison. */
    static bool testPrimitiveGradients(bool radians = false);
};

class AtomPointerVector : public vector<Atom*> {
  public:
    // inherit a bunch of useful constructors from vector
    AtomPointerVector() : vector<Atom*>() { }
    AtomPointerVector(size_t sz, Atom* val = NULL) : vector<Atom*>(sz, val) { }
    AtomPointerVector(const AtomPointerVector& other) : vector<Atom*>(other) { }
    AtomPointerVector(const vector<Atom*>& other) : vector<Atom*>(other) { }
    void copyCoordinates(const AtomPointerVector& other);

    using vector<Atom*>::push_back;    // base push_back of vector class
    void push_back(const Residue& R);  // overloaded push_back for Residues
    void push_back(const Residue* R) { push_back(*R); }
    void push_back(const vector<Atom*>& atoms);  // overloaded push_back for Residues

    CartesianPoint getGeometricCenter();
    void getGeometricCenter(mstreal& xc, mstreal& yc, mstreal& zc);
    void center();
    mstreal radiusOfGyration();
    // computes the smallest radius of a sphere centered at the centroid of the
    // atom vector, which encloses the set of all atoms.
    mstreal boundingSphereRadiusCent();
    void deletePointers();

    AtomPointerVector clone() const;
    void clone(AtomPointerVector& into) const;
    AtomPointerVector subvector(int beg, int end); // returns the range [beg, end)

    AtomPointerVector& operator+=(const AtomPointerVector& rhs);
    AtomPointerVector& operator-=(const AtomPointerVector& rhs);
    AtomPointerVector& operator/=(const mstreal& s);
    AtomPointerVector& operator*=(const mstreal& s);

    void write(ostream& _os) const; // write AtomPointerVector to a binary stream (parent information will be lost)
    void read(istream& _is);  // read AtomPointerVector from a binary stream
    friend ostream & operator<<(ostream &_os, const AtomPointerVector& _atoms);
};
ostream & operator<<(ostream &_os, const AtomPointerVector& _atoms); // this just to silence a silly compiler warning

class expressionTree {
  public:
    enum selProperty { RESID = 1, RESNAME, ICODE, CHAIN, SEGID, NAME, ALL };      // selectable properties
    enum logicalOp { AND = 1, OR, NOT, BYRES, BYCHAIN, IS, AROUND, RANGE };  // logical operators

    expressionTree(logicalOp _op = logicalOp::IS) { op = _op; }
    ~expressionTree() {
      for (int i = 0; i < children.size(); i++) delete children[i];
    }

    void setLogicalOperator(logicalOp _op) { op = _op; }
    void setProperty(selProperty _type) { type = _type; }
    void setNum(int _num) { num = _num; }
    void setNums(const vector<int> _nums) { nums = _nums; }
    void setNumSet(const set<int> _numset) { numset = _numset; }
    void addToNumSet(const vector<int> _newnums) {
      for (int newnum : _newnums) numset.insert(newnum);
    }
    void setVal(mstreal _val) { val = _val; }
    void setString(string _str) { str = _str; }
    void addChild(expressionTree* subtree) { children.push_back(subtree); }
    logicalOp getLogicalOperator() { return op; }
    selProperty getProperty() { return type; }
    int getNum() { return num; }
    vector<int> getNums() { return nums; }
    set<int> getNumSet() { return numset; }
    bool hasNums() { return !nums.empty(); }
    bool hasNumSet() { return !numset.empty(); }
    bool inNumSet(int i) const { return (numset.find(i) != numset.end()); }
    int getNumByIdx(int i) { return nums[i]; }
    mstreal getVal() { return val; }
    string getString() { return str; }
    int numChildren() { return children.size(); }
    expressionTree* getChild(int i) { return children[i]; }

  private:
    logicalOp op;
    selProperty type;
    int num; vector<int> nums; set<int> numset;
    mstreal val;
    string str;
    vector<expressionTree*> children;
};

class selector {
  public:
    selector(const Structure& S);
    AtomPointerVector select(string selStr);
    vector<Residue*> selectRes(string selStr);
    void select(expressionTree* tree, AtomPointerVector& sel);
    expressionTree* buildExpressionTree(string selStr);
    AtomPointerVector around(AtomPointerVector& selAtoms, mstreal dcut);
    AtomPointerVector byRes(AtomPointerVector& selAtoms);
    AtomPointerVector byChain(AtomPointerVector& selAtoms);
    AtomPointerVector invert(AtomPointerVector& selAtoms);
    AtomPointerVector intersect(AtomPointerVector& selA, AtomPointerVector& selB);
    AtomPointerVector combine(AtomPointerVector& selA, AtomPointerVector& selB);

  private:
    string getNextSelectionToken(string& selStr);
    vector<Atom*> atoms;
    vector<Residue*> atomResidues, residues;
    vector<Chain*> atomChains;
};

class RMSDCalculator {
 public:
    RMSDCalculator() {}
    ~RMSDCalculator() {}

    // getters
    mstreal lastRMSD() { return sqrt(_res/_n); }
    vector<mstreal> lastTranslation();
    vector<vector<mstreal> > lastRotation();

    // calculate optimal superposition and the resulting RMSD, applying transformation to given atoms
    bool align(const vector<Atom*> &_align, const vector<Atom*> &_ref, vector<Atom*>& _moveable);
    bool align(const vector<Atom*> &_align, const vector<Atom*> &_ref, Structure& _moveable);

    // simply apply the previously obtained optimal rotatoin/translation
    void applyLastTransformation(vector<Atom*>& _moveable);

    // quickly calculate RMSD upon optimal superposition without generating the rotation matrix
    mstreal bestRMSD(const vector<Atom*> &_align, const vector<Atom*> &_ref, bool setTransRot = false, bool* _suc = NULL);
    // calculates many RMSDs from a vector of atom pointer vectors; the first int determines where this starts the comparison along the vector of atom pointer vectors, and the second determines how far from the end it goes
    vector<mstreal> bestRMSD(const vector<vector<Atom*>> &_align, const vector<Atom*> &_ref, int, int);
    mstreal bestResidual(const vector<Atom*> &_align, const vector<Atom*> &_ref, bool setTransRot = false, bool* _suc = NULL);

    // in-place RMSD (no transformations)
    static mstreal rmsd(const vector<Atom*>& A, const vector<Atom*>& B);
    static mstreal rmsd(const Structure& A, const Structure& B);

    // -- Craig's cutoff function
    /* L stores the number of residues in each chain, assumed to be contiguous.
     * Residues in different chains are assumed to be statistically independent. */
    static mstreal rmsdCutoff(const vector<int>& L, mstreal rmsdMax = 1.1, mstreal L0 = 15);
    /* Takes every chain in S to be independent. */
    static mstreal rmsdCutoff(const Structure& S, mstreal rmsdMax = 1.1, mstreal L0 = 15);
    /* I stores residue indices of residues in each chain, so that non-contiguous
     * residues within a chain can be treated. Residues in different chains are
     * still assumed to be independent. */
    static mstreal rmsdCutoff(const vector<vector<int> >& I, mstreal rmsdMax = 1.1, mstreal L0 = 15);
    /* RMSD is computed for the fragment comprising residues from the given
     * Structure S, whose indides are given in J. Like in the previous function,
     * for residues located within the same chain of S, their relative locations
     * are accounted for when computing degrees of freedom. */
    static mstreal rmsdCutoff(const vector<int>& J, const Structure& S, mstreal rmsdMax = 1.1, mstreal L0 = 15);

    template <class T>
    mstreal qcpRMSD(const T& A, const T& B, bool setTransform = false, bool setResiduals = false);

    template <class T>
    mstreal qcpRMSDGrad(const T& A, const T& B, vector<mstreal>& grad);

    /* Tests QCP implementation for RMSD and RMSD gradient calculation. */
    static bool testQCP(bool testGrad = false);

 protected:
    // implemetation of Kabsch algoritm for optimal superposition
    bool Kabsch(const vector<Atom*> &_align, const vector<Atom*> &_ref, int mode);

 private:
    mstreal _res;
    int _n;
    vector<vector<mstreal> > residuals;
    mstreal t[3];    // translation vector
    mstreal u[3][3]; // rotation matrix

};

class ProximitySearch {
  public:
    ProximitySearch() { xlo = ylo = zlo = xhi = yhi = zhi = xbw = ybw = zbw = 0.0; N = 0;}
    ProximitySearch(const ProximitySearch& ps);
    ProximitySearch(mstreal _xlo, mstreal _ylo, mstreal _zlo, mstreal _xhi, mstreal _yhi, mstreal _zhi, int _N = 20);
    ProximitySearch(const AtomPointerVector& _atoms, int _N, bool _addAtoms = true, vector<int>* tags = NULL, mstreal pad = 0);
    ProximitySearch(const AtomPointerVector& _atoms, mstreal _characteristicDistance, bool _addAtoms = true, vector<int>* tags = NULL, mstreal pad = 0);
    ProximitySearch(const vector<CartesianPoint>& _points, mstreal _characteristicDistance, bool _addPoints = true, vector<int>* tags = NULL, mstreal pad = 0);
    ~ProximitySearch();
    ProximitySearch& operator=(const ProximitySearch& ps);

    mstreal getXLow() { return xlo; }
    mstreal getYLow() { return ylo; }
    mstreal getZLow() { return zlo; }
    mstreal getXHigh() { return xhi; }
    mstreal getYHigh() { return yhi; }
    mstreal getZHigh() { return zhi; }
    int pointSize() { return pointList.size(); }
    CartesianPoint& getPoint(int i) { return *(pointList[i]); }
    int getPointTag(int i) { return pointTags[i]; }
    mstreal distance(int i, int j) { return pointList[i]->distance(pointList[j]); }

    void reinitBuckets(int _N);
    void addPoint(const CartesianPoint& _p, int tag);
    void addPoint(mstreal xc, mstreal yc, mstreal zc, int tag);
    void addAtoms(AtomPointerVector& apv, vector<int>* tags = NULL);
    void dropAllPoints();
    bool isPointWithinGrid(CartesianPoint _p);
    void pointBucket(CartesianPoint* p, int* i, int* j, int* k) { pointBucket(p->getX(), p->getY(), p->getZ(), i, j, k); }
    void pointBucket(const CartesianPoint& p, int* i, int* j, int* k) { pointBucket(p.getX(), p.getY(), p.getZ(), i, j, k); }
    void pointBucket(mstreal px, mstreal py, mstreal pz, int* i, int* j, int* k);
    void limitIndex(int *ind);
    mstreal limitX(mstreal x);
    mstreal limitY(mstreal y);
    mstreal limitZ(mstreal z);
    mstreal gridSpacingX() { return (xhi - xlo)/N; }
    mstreal gridSpacingY() { return (yhi - ylo)/N; }
    mstreal gridSpacingZ() { return (zhi - zlo)/N; }
    static void calculateExtent(const AtomPointerVector& _atoms, mstreal& _xlo, mstreal& _ylo, mstreal& _zlo, mstreal& _xhi, mstreal& _yhi, mstreal& _zhi);
    static void calculateExtent(const Structure& S, mstreal& _xlo, mstreal& _ylo, mstreal& _zlo, mstreal& _xhi, mstreal& _yhi, mstreal& _zhi);
    static void calculateExtent(const vector<CartesianPoint>& points, mstreal& _xlo, mstreal& _ylo, mstreal& _zlo, mstreal& _xhi, mstreal& _yhi, mstreal& _zhi);

    bool pointsWithin(const CartesianPoint& c, mstreal dmin, mstreal dmax, vector<int>* list = NULL, bool byTag = false);
    vector<int> getPointsWithin(const CartesianPoint& c, mstreal dmin, mstreal dmax, bool byTag = false);
    int numPointsWithin(const CartesianPoint& c, mstreal dmin, mstreal dmax) {
      vector<int> closeOnes; pointsWithin(c, dmin, dmax, &closeOnes); return closeOnes.size();
    }

    // Returns true if the grid of the current ProximitySearch object overlaps
    // that of the ProximitySearch specified by more than the padding given
    bool overlaps(ProximitySearch& other, mstreal pad = 0);

  protected:
    void setBinWidths();
    void calculateExtent(const AtomPointerVector& _atoms) { ProximitySearch::calculateExtent(_atoms, xlo, ylo, zlo, xhi, yhi, zhi); }
    void calculateExtent(const vector<CartesianPoint>& _points) { ProximitySearch::calculateExtent(_points, xlo, ylo, zlo, xhi, yhi, zhi); }

  private:
    int N; // dimension of bucket list is N x N x N

    mstreal xlo, ylo, zlo, xhi, yhi, zhi, xbw, ybw, zbw; // extents of coordinates

    /* Each bucket is a vector of point indices (zero-initiated). These
     * indices are into two vectors: a vector of 3D points (CartesianPoint
     * pointers) and a vector of point tags. In this way, it is easy to go
     * from a bucket into points as well as from a point index into its point
     * or tag. One can also go from a point to its bucket location, by doing
     * a simple computations on the coordinates of the point via poitBucket(). */
    vector<vector<vector<vector<int> > > > buckets;
    vector<CartesianPoint*> pointList;
    vector<int> pointTags;
    vector<vector<int>*> fullBuckets; // keeps track of buckets
};

template<class T>
class DecoratedProximitySearch : public ProximitySearch {
  public:
    DecoratedProximitySearch(mstreal _xlo, mstreal _ylo, mstreal _zlo, mstreal _xhi, mstreal _yhi, mstreal _zhi, int _N = 20) :
      ProximitySearch(_xlo, _ylo, _zlo, _xhi, _yhi, _zhi, _N) {}
    DecoratedProximitySearch(AtomPointerVector& _atoms, int _N, vector<T>& _tags, mstreal pad = 0) :
      ProximitySearch(_atoms, _N, true, NULL, pad) {
      tags = _tags;
      checkSizes();
    }
    DecoratedProximitySearch(AtomPointerVector& _atoms, int _N, mstreal pad = 0) : ProximitySearch(_atoms, _N, false, NULL, pad) {}
    DecoratedProximitySearch(AtomPointerVector& _atoms, mstreal _characteristicDistance, vector<T>& _tags, mstreal pad = 0) :
      ProximitySearch(_atoms, _characteristicDistance, true, NULL, pad) {
      tags = _tags;
      checkSizes();
    }
    DecoratedProximitySearch(AtomPointerVector& _atoms, mstreal _characteristicDistance, mstreal pad = 0) :
      ProximitySearch(_atoms, _characteristicDistance, false, NULL, pad) { }

    T getPointTag(int i) { return tags[this->ProximitySearch::getPointTag(i)]; }
    void addPoint(CartesianPoint _p, T tag) {
      this->ProximitySearch::addPoint(_p, tags.size());
      tags.push_back(tag);
    }

    vector<T> getPointsWithin(const CartesianPoint& c, mstreal dmin, mstreal dmax) {
      vector<int> inds;
      this->ProximitySearch::pointsWithin(c, dmin, dmax, &inds, true);
      vector<T> ret(inds.size());
      for (int i = 0; i < inds.size(); i++) ret[i] = tags[inds[i]];
      return ret;
    }

    void getPointsWithin(const CartesianPoint& c, mstreal dmin, mstreal dmax, vector<T>* list) {
      vector<int> inds;
      this->ProximitySearch::pointsWithin(c, dmin, dmax, &inds, true);
      list->resize(inds.size());
      for (int i = 0; i < inds.size(); i++) (*list)[i] = tags[inds[i]];
    }

    vector<int> getPointsWithinIndices(const CartesianPoint& c, mstreal dmin, mstreal dmax) {
      return this->ProximitySearch::getPointsWithin(c, dmin, dmax, true);
    }

  protected:
    bool checkSizes();

  private:
    vector<T> tags;
};


class Clusterer {
  public:
    Clusterer(bool _flag = true) { optimAlign = _flag; }
    void optimizeAlignments(bool _flag) { optimAlign = _flag; }
    bool getOptimizeAlignments() { return optimAlign; }

    /* Will greedy cluster the given set of units (must all have the same number of atoms),
     * using the given RMSD cutoff, while making sure that no more than ~Nmax x Nmax RMSD
     * computations are done per iteration. So, if the number of units is below Nmax, a
     * normal greedy clustering will be performed. However, if above Nmax, then will try
     * be further greedy and find best centroids without ever doing all-by-all comparisons.
     * Output a vector of clusters, in decreasing size, each of which is a vector of
     * indices of elements mapping to that cluster. The centroid is always listed first
     * in each cluster. */
    vector<vector<int> > greedyCluster(const vector<vector<Atom*> >& units, mstreal rmsdCut, int Nmax = 10000, mstreal coverage = 1.0, int maxClusters = -1, bool verbose = true);

    /* Perform k-means clustering of a point cloud in arbitrary dimension, using
     * Euclidean distance as the metric. Returns a list of clusters of size k,
     * where each cluster is represented by a vector point indices. */
    vector<vector<int> > kmeans(const vector<CartesianPoint>& points, int k, int Ntrials = 1, int Niter = 10);

  protected:
    // these functions are protected because they assume that the cache of pre-
    // computed RMSDs is in a good state (so don't want external calls)
    vector<vector<int> > greedyClusterBruteForce(const vector<vector<Atom*> >& units, set<int> remIndices, mstreal rmsdCut, int nClusts = -1);
    vector<int> elementsWithin(const vector<vector<Atom*> >& units, set<int>& remIndices, const vector<Atom*>& fromUnit, mstreal rmsdCut);
    set<int> randomSubsample(set<int>& indices, int N);

  private:
    // won't cache for now (need careful memory management)
    map<int, map<int, mstreal > > coputedRMSDs;
    bool optimAlign;
    RMSDCalculator rCalc;
};

/* A simple map data structure based on storing a sorted list of keys (sorted by
 * the less-than operator). This is much more memory efficient than std::map and
 * the time complexity for lookup is still logarithmic. However, time complexity
 * for insertion (when inserting one by one) is still linear with size, while if
 * all the keys are specified initially, then it is amortized logarithmic too. */
template<class keyType, class valType>
class simpleMap {
  public:
    simpleMap() {}
    simpleMap(const vector<keyType>& _keys, const vector<valType>& _vals);
    size_t size() const { return keys.size(); }

    // access and modification
    valType at(const keyType& key);
    valType& operator[](const keyType& key); // will create entry if key missing, setting the value to the default value
    valType& value(int idx) { return vals[idx]; } // get value by index
    keyType key(int idx) const { return keys[idx]; } // get key by index
    int find(const keyType& key); // returns a negative index if key is not found
    void insert(const keyType& key, const valType& val);
    void erase(const keyType& key);
    void clear() { keys.clear(); vals.clear(); }

    /* Same meaning as for corresponding tightvector<T> functions, applied to keys. */
    int getUpperBound(const keyType& k) { return lower_bound(keys.begin(), keys.end(), k) - keys.begin(); }
    int getLowerBound(const keyType& k) { return upper_bound(keys.begin(), keys.end(), k) - keys.begin(); }

  protected:
    void insert(const keyType& key, const valType& val, int i);

  private:
    vector<keyType> keys; // sorted keys
    vector<valType> vals; // values, stored in corresponding order
};

/* A memory-efficient vector class (leaves no extra space). */
template<class T>
class tightvector {
  public:
    tightvector() { vec = NULL; L = 0; }
    tightvector(int _L, const T& val = T());
    tightvector(const tightvector& other);
    tightvector(const vector<T>& other);
    ~tightvector();

    size_t size() const { return L; }
    bool empty() const { return (L == 0); }
    void resize(int len, const T& val = T());
    void push_back(const T& val);

    /* Inserts value into the vector at the given index by shifting all values
     * after this index to the right by one and overwriting the value at the
     * index with the given value. */
    void insert(const T& val, int idx);

    /* Insert value into vector, assuming the vector is currently sorted (based
     * on the T::operator<), and in a way that leaves the vector sorted. */
    int insertSorted(const T& val);

    /* Linear search for value. Returns the first element (left to right) which
     * matches the value, or a negative index if value not found. */
    int find(const T& val);

    /* Binary search that assumes the vector to be sorted according to T::operator<.
     * The index of the first candidate element with the sought value is returned
     * or a negative index if value not found. */
    int findSorted(const T& val);

    /* Both functions assume the vector to be sorted according to T::operator<.
     * getUpperBound -- finds the first element (left to right) that evaluates
     * to > val. If no such element exists (i.e., all elements are
     * not greater than val), the size of the current vector is returned.
     * getLowerBound -- finds the first element (left to right) that does NOT
     * evaluate to < val. If no such element exists (i.e., all elements are
     * less than val), the size of the current vector is returned. */
    int getUpperBound(const T& val);
    int getLowerBound(const T& val);

    void clear();

    // NOTE: subscript operators do not check for out-of-bounds indices!
    T& operator[](int idx) { return vec[idx]; }
    T operator[](int idx) const { return vec[idx]; }

    tightvector<T>& operator=(const tightvector<T>& other);
    tightvector<T>& operator=(const vector<T>& other);

    // cast into regular vector
    operator vector<T>() const;

  private:
    T* vec;
    int L;
};
}

/* A simple and accurate timer class built on top of chrono::high_resolution_clock */
class MstTimer {
  public:
    enum timeUnits { sec = 0, msec, usec, nsec };

    MstTimer() {
      running = false;
      begin = chrono::high_resolution_clock::now();
      end = begin;
      elapsed = end - begin;
    }
    bool isRunning() { return running; }
    bool start(bool resume = false) {
      if (!running) {
        running = true;
        if (!resume) elapsed = begin - begin; // simulate a zero
        begin = chrono::high_resolution_clock::now();
        return true;
      }
      return false;
    }
    bool stop() {
      if (running) {
        running = false;
        end = chrono::high_resolution_clock::now();
        elapsed += (end - begin);
        return true;
      }
      return false;
    }
    int getDuration(timeUnits units = timeUnits::sec) {
      chrono::high_resolution_clock::duration dt;
      if (!running) {
        dt = elapsed;
      } else {
        dt = chrono::high_resolution_clock::now() - begin;
      }
      switch (units) {
        case timeUnits::sec:
          return chrono::duration_cast<std::chrono::seconds>(dt).count();
        case timeUnits::msec:
          return chrono::duration_cast<std::chrono::milliseconds>(dt).count();
        case timeUnits::usec:
          return chrono::duration_cast<std::chrono::microseconds>(dt).count();
        case timeUnits::nsec:
          return chrono::duration_cast<std::chrono::nanoseconds>(dt).count();
      }
      return 0; // ERROR: should never occur
    }

  private:
    bool running;
    chrono::high_resolution_clock::time_point begin, end;
    chrono::high_resolution_clock::duration elapsed;
};

/* Utilities class, with a bunch of useful static functions, is defined outside of the MST namespace because:
 * 1) it really represents a different beast, not an MST type
 * 2) some of its functions (like assert) are likely to clash with function names in other project
 */
class MstUtils {
  public:
    template <class F>
    static void openFile(F& fs, string filename, ios_base::openmode mode = ios_base::in, string from = "");
    static void fileToArray(const string& _filename, vector<string>& lines); // reads lines from the file and appends them to the given vector
    static vector<string> fileToArray(const string& _filename) { vector<string> lines; fileToArray(_filename, lines); return lines; }
    static FILE* openFileC (const char* filename, const char* mode, string from = "");
    static vector<string> trim(const vector<string>& strings, string delimiters = " \t\n\v\f\r");
    static string trim(const string& str, string delimiters = " \t\n\v\f\r");
    static void warn(const string& message, string from = "");
    static void error(const string& message, string from = "", int code = -1);
    static void assertCond(bool condition, string message = "error: assertion failed", string from = "", int exitCode = -1);
    static string uc(const string& str);                        // returns an upper-case copy of the input string
    static string lc(const string& str);                        // returns an lower-case copy of the input string
    static bool stringsEqual(const string& A, const string& B, bool caseInsensitive = true);
    static string wrapText(const string& message, int width, int leftSkip = 0, int startingOffset = 0);
    static char* copyStringC(const char* str);
    static int toInt(const string& num, bool strict = true);
    static bool isInt(const string& num);
    static MST::mstreal toReal(const string& num, bool strict = true);
    static bool isReal(const string& num);
    static MST::mstreal mod(MST::mstreal num, MST::mstreal den);
    static MST::mstreal sign(MST::mstreal val) { return (val > 0) ? 1.0 : ((val < 0) ? -1.0 : 0.0); }
    static string nextToken(string& str, string delimiters = " ", bool skipTrailingDelims = true);
    static string nextQuoteAwareToken(string& str, string quoteMarks = "'\"", string delimiters = " ", bool skipTrailingDelims = true);
    static string escape(string str, string special);
    static string unescape(string str);
    static vector<string> split(const string& str, string delimiters = " ", bool skipTrailingDelims = true);
    static string removeComment(const string& str, string commStart = "#");
    static vector<MST::mstreal> splitToReal(const string& str, string delimiters = " ", bool skipTrailingDelims = true, bool strict = true);
    static vector<int> splitToInt(const string& str, string delimiters = " ", bool skipTrailingDelims = true, bool strict = true);
    static string join(const string& delim, const vector<string>& words);
    static string readNullTerminatedString(fstream& ifs);
    static string getDate();
    static vector<pair<int, int> > splitTasks(int numTasks, int numJobs);
    static void setSignalHandlers();
    static void errorHandler(int sig);

    // random numbers
    static unsigned seedRandEngine(unsigned seed = std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count()) { mt.seed(seed); return seed; }
    static mt19937& randEngine() { return mt; }

    // returns a random number in the range [lower, upper]
    static int randInt(int lower, int upper);
    // returns a random number in the range [0, upper) (convenient for generating random array subscripts)
    static int randInt(int upper) { return randInt(0, upper - 1); }
    // random number in the unit range 0 and 1
    static MST::mstreal randUnit(MST::mstreal mi = 0, MST::mstreal ma = 1.0);
    // normally-distributed random number with mean mu and standard deviation sig
    static MST::mstreal randNormal(MST::mstreal mu = 0.0, MST::mstreal sig = 1.0);

    template <class T>
    static T pow2(const T& A) { return A*A; }
    template <class T>
    static string toString(const T& obj) { return toString(&obj); }
    template <class T>
    static string toString(const T* obj);
    template <class T>
    static vector<int> sortIndices(const vector<T>& vec, bool descending = false);
    template <class T>
    static vector<T> subVector(const vector<T>& vec, const vector<int>& inds);
    template <class T1, class T2>
    static vector<T1> keys(const map<T1, T2>& _map);
    template <class T>
    static vector<T> keys(const set<T>& _set);
    template <class T1, class T2>
    static vector<T2> values(map<T1, T2>& _map);
    template <class T>
    static string vecToString(const vector<T>& vec, string del = " ");
    template <class T>
    static string vecPtrToString(const vector<T*>& vec, string del = " ");
    template <class T>
    static T min(const T& a, const T& b);
    template <class T>
    static T max(const T& a, const T& b);
    template <class T>
    static T min(const vector<T>& vec, int beg = -1, int end = -1, int* minIndex = NULL);
    template <class T>
    static T max(const vector<T>& vec, int beg = -1, int end = -1, int* maxIndex = NULL);
    template <class T>
    static bool closeEnough(const T& a, const T& b, const T& epsilon = std::numeric_limits<T>::epsilon());
    template <class T>
    static set<T> contents(const vector<T>& vec);
    template <class T>
    static map<T, int> indexMap(const vector<T>& vec);

    /* Returns a vector with all of the elements in vector A, in their original
     * relative order, that do not appear in vector B. */
    template <class T>
    static vector<T> setdiff(const vector<T>& A, const vector<T>& B);

    /* Returns a vector with all of the elements in vector A, in their original
     * relative order, followed by any elements in B that are not in A, in their
     * original order. */
    template <class T>
    static vector<T> setunion(const vector<T>& A, const vector<T>& B);

    /* Returns a vector with all of the elements that are both in vector A and B,
    * in the order that they appear in A. */
    template <class T>
    static vector<T> setintersect(const vector<T>& A, const vector<T>& B);

    template <class T>
    static vector<T> range(const T& from, const T& to, const T& step = 1);

    // randomly shuffle the array in place using the Fisher-Yates shuffle
    template <class T>
    static void shuffle(vector<T>& vec);

    // read to/from binary file
    template <class T>
    static void writeBin(ostream& ofs, const T& val); // binary output for primitive types
    static void writeBin(ostream& ofs, const string& str) { ofs << str << '\0'; } // special overload for strings
    template <class T, class U>
    static void writeBin(ostream& ofs, const pair<T, U>& p) { writeBin(ofs, p.first); writeBin(ofs, p.second); }
    template <class T>
    static void writeBin(ostream& ofs, const vector<T>& vec);
    template <class T>
    static void writeBin(ostream& ofs, const set<T>& s) { writeBin(ofs, vector<T>(s.begin(), s.end())); }
    template <class K, class V>
    static void writeBin(ostream& ofs, const map<K, V>& m);
    static void writeBin(ostream& ofs, const MST::Structure& S) { S.writeData(ofs); }
    template <class T>
    static void readBin(istream& ifs, T& val);
    static void readBin(istream& ifs, string& str) { getline(ifs, str, '\0'); };
    template <class T, class U>
    static void readBin(istream& ifs, pair<T, U>& p) { readBin(ifs, p.first); readBin(ifs, p.second); }
    template <class T>
    static void readBin(istream& ifs, vector<T>& vec);
    template <class T>
    static void readBin(istream& ifs, set<T>& s);
    template <class K, class V>
    static void readBin(istream& ifs, map<K, V>& m);
    static void readBin(istream& ifs, MST::Structure& S) { S.readData(ifs); }

    private:
      static mt19937 mt; // A Mersenne Twister pseudo-random generator of 32-bit numbers with a state size of 19937 bits.
};

template <class F>
void MstUtils::openFile(F& fs, string filename, ios_base::openmode mode, string from) {
  fs.open(filename.c_str(), mode);
  if (!fs.is_open()) {
    if (!from.empty()) from += " -> ";
    MstUtils::error("could not open file '" + filename + "'", from + "MstUtils::openFile");
  }
}


template <class T>
string MstUtils::toString(const T* obj) {
  stringstream ss;
  ss << *obj;
  return ss.str();
}

template <class T>
vector<int> MstUtils::sortIndices(const vector<T>& vec, bool descending) {
  vector<int> sortedIndices(vec.size(), 0);
  for (int i = 0; i < vec.size(); i++) sortedIndices[i] = i;
  if (descending) {
    sort(sortedIndices.begin(), sortedIndices.end(), [&vec](size_t i1, size_t i2) {return vec[i2] < vec[i1];});
  } else {
    sort(sortedIndices.begin(), sortedIndices.end(), [&vec](size_t i1, size_t i2) {return vec[i1] < vec[i2];});
  }
  return sortedIndices;
}

template <class T>
vector<T> MstUtils::subVector(const vector<T>& vec, const vector<int>& inds) {
  vector<T> sub(inds.size());
  for (int i = 0; i < inds.size(); i++) sub[i] = vec[inds[i]];
  return sub;
}

template <class T1, class T2>
vector<T1> MstUtils::keys(const map<T1, T2>& _map) {
  vector<T1> K(_map.size());
  int k = 0;
  for (auto it = _map.begin(); it != _map.end(); ++it, ++k) {
    K[k] = it->first;
  }
  return K;
}

template <class T>
vector<T> MstUtils::keys(const set<T>& _set) {
  vector<T> K(_set.size());
  int k = 0;
  for (auto it = _set.begin(); it != _set.end(); ++it, ++k) {
    K[k] = *it;
  }
  return K;
}

template <class T1, class T2>
vector<T2> MstUtils::values(map<T1, T2>& _map) {
  vector<T2> V(_map.size());
  int k = 0;
  for (typename map<T1, T2>::iterator it = _map.begin(); it != _map.end(); ++it, ++k) {
    V[k] = it->second;
  }
  return V;
}

template <class T>
string MstUtils::vecToString(const vector<T>& vec, string del) {
  string str;
  for (int i = 0; i < vec.size(); i++) {
    str += MstUtils::toString(vec[i]);
    if (i != vec.size() - 1) str += del;
  }
  return str;
}

template <class T>
string MstUtils::vecPtrToString(const vector<T*>& vec, string del) {
  string str;
  for (int i = 0; i < vec.size(); i++) {
    str += MstUtils::toString(*(vec[i]));
    if (i != vec.size() - 1) str += del;
  }
  return str;
}

template <class T>
T MstUtils::min(const T& a, const T& b) {
  if (a < b) return a;
  return b;
}

template <class T>
T MstUtils::max(const T& a, const T& b) {
  if (a > b) return a;
  return b;
}

template <class T>
T MstUtils::min(const vector<T>& vec, int beg, int end, int* minIndex) {
  MstUtils::assertCond(vec.size() > 0, "empty vector passed!", "MstUtils::min(vector<T>&)");
  if (beg < 0) beg = 0;
  if (end < 0) end = vec.size()-1;
  MstUtils::assertCond(beg <= end, "beg = " + MstUtils::toString(beg) + " and end = " + MstUtils::toString(end), "MstUtils::min(vector<T>&)");
  MstUtils::assertCond(end < vec.size(), "end index out of bounds, " + MstUtils::toString(end), "MstUtils::min(vector<T>&)");

  T mv = vec[beg];
  if (minIndex != NULL) *minIndex = beg;
  for (int i = beg; i <= end; i++) {
    if (vec[i] < mv) {
      mv = vec[i];
      if (minIndex != NULL) *minIndex = i;
    }
  }
  return mv;
}

template <class T>
T MstUtils::max(const vector<T>& vec, int beg, int end, int* maxIndex) {
  MstUtils::assertCond(vec.size() > 0, "empty vector passed!", "MstUtils::max(vector<T>&)");
  if (beg < 0) beg = 0;
  if (end < 0) end = vec.size()-1;
  MstUtils::assertCond(beg <= end, "beg = " + MstUtils::toString(beg) + " and end = " + MstUtils::toString(end), "MstUtils::max(vector<T>&)");
  MstUtils::assertCond(end < vec.size(), "end index out of bounds, " + MstUtils::toString(end), "MstUtils::max(vector<T>&)");

  T mv = vec[beg];
  if (maxIndex != NULL) *maxIndex = beg;
  for (int i = beg; i <= end; i++) {
    if (vec[i] > mv) {
      mv = vec[i];
      if (maxIndex != NULL) *maxIndex = i;
    }
  }
  return mv;
}

template <class T>
bool MstUtils::closeEnough(const T& a, const T& b, const T& epsilon) {
  return (a - b > -epsilon) && (a - b < epsilon);
}

template <class T>
void MstUtils::shuffle(vector<T>& vec) {
  int i, j; T tmp;
  for (i = vec.size() - 1; i >= 0; i--) {
    j = MstUtils::randInt(0, i);
    if (i == j) continue;
    tmp = vec[i]; vec[i] = vec[j]; vec[j] = tmp;
  }
}

template <class T>
set<T> MstUtils::contents(const vector<T>& vec) {
  set<T> cont;
  for (int i = 0; i < vec.size(); i++) cont.insert(vec[i]);
  return cont;
}

template <class T>
map<T, int> MstUtils::indexMap(const vector<T>& vec) {
  map<T, int> indices;
  for (int i = 0; i < vec.size(); i++) indices[vec[i]] = i;
  return indices;
}

template <class T>
void MstUtils::writeBin(ostream& ofs, const T& val) {
  ofs.write((char*) (&val), sizeof(T));
}

template <class T>
void MstUtils::writeBin(ostream& ofs, const vector<T>& vec) {
  MstUtils::writeBin(ofs, (int) vec.size());
  for (int i = 0; i < vec.size(); i++) MstUtils::writeBin(ofs, vec[i]);
}

template <class K, class V>
void MstUtils::writeBin(ostream& ofs, const map<K, V>& m) {
  MstUtils::writeBin(ofs, (int) m.size());
  for (auto i = m.begin(); i != m.end(); ++i) {
    MstUtils::writeBin(ofs, i->first);
    MstUtils::writeBin(ofs, i->second);
  }
}

template <class T>
void MstUtils::readBin(istream& ifs, T& val) {
  ifs.read((char*) (&val), sizeof(T));
}

template <class T>
void MstUtils::readBin(istream& ifs, vector<T>& vec) {
  int len; MstUtils::readBin(ifs, len);
  vec.resize(len);
  for (int i = 0; i < vec.size(); i++) MstUtils::readBin(ifs, vec[i]);
}

template <class T>
void MstUtils::readBin(istream& ifs, set<T>& s) {
  vector<T> vec; MstUtils::readBin(ifs, vec);
  s.clear();
  s.insert(vec.begin(), vec.end());
}

template <class K, class V>
void MstUtils::readBin(istream& ifs, map<K, V>& m) {
  int len; MstUtils::readBin(ifs, len);
  m.clear();
  for (int i = 0; i < len; i++) {
    K key; MstUtils::readBin(ifs, key);
    V val; MstUtils::readBin(ifs, val);
    m[key] = val;
  }
}

template <class T>
vector<T> MstUtils::setdiff(const vector<T>& A, const vector<T>& B) {
  vector<T> diff;
  set<T> setB = MstUtils::contents(B);
  for (int i = 0; i < A.size(); i++) {
    if (setB.find(A[i]) == setB.end()) diff.push_back(A[i]);
  }
  return diff;
}

template <class T>
vector<T> MstUtils::setunion(const vector<T>& A, const vector<T>& B) {
  vector<T> U = A;
  set<T> setA = MstUtils::contents(A);
  for (int i = 0; i < B.size(); i++) {
    if (setA.find(B[i]) == setA.end()) U.push_back(B[i]);
  }
  return U;
}


template <class T>
vector<T> MstUtils::setintersect(const vector<T>& A, const vector<T>& B) {
  vector<T> inter;
  set<T> setB = MstUtils::contents(B);
  for (int i = 0; i < A.size(); i++) {
    if (setB.find(A[i]) != setB.end()) inter.push_back(A[i]);
  }
  return inter;
}

template <class T>
vector<T> MstUtils::range(const T& from, const T& to, const T& step) {
  vector<T> r;
  for (T v = from; v < to; v += step) r.push_back(v);
  return r;
}

using namespace MST;

/* --------- simpleMap --------- */
template<class keyType, class valType>
simpleMap<keyType, valType>::simpleMap(const vector<keyType>& _keys, const vector<valType>& _vals) {
  if (_keys.size() != _vals.size()) MstUtils::error("different number of keys and values specified", "simpleMap::simpleMap");
  keys.resize(_keys.size());
  vals.resize(_vals.size());
  vector<int> sInds = MstUtils::sortIndices(_keys);
  for (int i = 0; i < _keys.size(); i++) {
    keys[i] = _keys[sInds[i]];
    vals[i] = _vals[sInds[i]];
    if ((i > 0) && (keys[i] == keys[i-1])) MstUtils::error("some of the specified keys are identical!", "simpleMap::simpleMap");
  }
};

template<class keyType, class valType>
valType simpleMap<keyType, valType>::at(const keyType& key) {
  auto it = lower_bound(keys.begin(), keys.end(), key); // first element which does not compare less than key
  if ((it == keys.end()) || (*it > key)) MstUtils::error("key not in map", "simpleMap::at");
  return vals[it - keys.begin()];
}

template<class keyType, class valType>
valType& simpleMap<keyType, valType>::operator[](const keyType& key) {
  auto it = lower_bound(keys.begin(), keys.end(), key); // first element which does not compare less than key
  int idx = it - keys.begin();
  if ((it == keys.end()) || (*it > key)) insert(key, valType(), idx);
  return vals[idx];
}

template<class keyType, class valType>
int simpleMap<keyType, valType>::find(const keyType& key) {
  auto it = lower_bound(keys.begin(), keys.end(), key); // first element which does not compare less than key
  if ((it == keys.end()) || (*it > key)) return -1;
  return (int) (it - keys.begin());
}

template<class keyType, class valType>
void simpleMap<keyType, valType>::insert(const keyType& key, const valType& val, int i) {
  keys.insert(keys.begin() + i, key);
  vals.insert(vals.begin() + i, val);
}

template<class keyType, class valType>
void simpleMap<keyType, valType>::insert(const keyType& key, const valType& val) {
  auto it = lower_bound(keys.begin(), keys.end(), key); // first element which does not compare less than key
  int idx = it - keys.begin();
  if ((it == keys.end()) || (*it > key)) insert(key, val, idx);
  else vals[idx] = val;
}

template<class keyType, class valType>
void simpleMap<keyType, valType>::erase(const keyType& key) {
  auto it = lower_bound(keys.begin(), keys.end(), key); // first element which does not compare less than key
  if ((it == keys.end()) || (*it > key)) return;
  vals.erase(vals.begin() + (it - keys.begin()));
  keys.erase(it);
}


/* --------- tightvector --------- */
template<class T>
tightvector<T>::tightvector(int _L, const T& val) {
  L = _L;
  vec = new T[L];
  for (int i = 0; i < L; i++) vec[i] = val;
}

template<class T>
tightvector<T>::tightvector(const tightvector& other) {
  if (other.empty()) { L = 0; vec = NULL; }
  else {
    L = other.size();
    vec = new T[L];
    for (int i = 0; i < L; i++) vec[i] = other[i];
  }
}

template<class T>
tightvector<T>::tightvector(const vector<T>& other) {
  if (other.empty()) { L = 0; vec = NULL; }
  else {
    L = other.size();
    vec = new T[L];
    for (int i = 0; i < L; i++) vec[i] = other[i];
  }
}

template<class T>
tightvector<T>::~tightvector() {
  if (vec != NULL) delete[] (vec);
}

template<class T>
void tightvector<T>::resize(int len, const T& val) {
  T* vecNew = (len == 0) ? NULL : new T[len];
  int commonLen = MstUtils::min(len, L);
  for (int i = 0; i < commonLen; i++) vecNew[i] = vec[i];
  for (int i = L; i < len; i++) vecNew[i] = val;
  L = len;
  if (vec != NULL) delete(vec);
  vec = vecNew;
}

template<class T>
void tightvector<T>::clear() {
  if (vec != NULL) delete(vec);
  L = 0;
}

template<class T>
tightvector<T>& tightvector<T>::operator=(const tightvector<T>& other) {
  resize(other.size());
  for (int i = 0; i < size(); i++) vec[i] = other[i];
  return *this;
}

template<class T>
tightvector<T>& tightvector<T>::operator=(const vector<T>& other) {
  resize(other.size());
  for (int i = 0; i < size(); i++) vec[i] = other[i];
  return *this;
}

template<class T>
void tightvector<T>::push_back(const T& val) {
  resize(size() + 1);
  vec[size() - 1] = val;
}

template<class T>
int tightvector<T>::insertSorted(const T& val) {
  int idx = getUpperBound(val);
  insert(val, idx);
  return idx;
}

template<class T>
int tightvector<T>::find(const T& val) {
  for (int i = 0; i < size(); i++) {
    if (vec[i] == val) return i;
  }
  return -1;
}

template<class T>
int tightvector<T>::findSorted(const T& val) {
  int lo = 0, hi = size() - 1;
  while (lo < hi) {
    int mid = (lo + hi)/2;
    if (vec[mid] < val) lo = mid + 1;
    else if (val < vec[mid]) hi = mid - 1;
    else return mid;
  }
  if (lo > hi) return -1;
  return lo;
}

template<class T>
void tightvector<T>::insert(const T& val, int idx) {
  resize(size() + 1);
  for (int i = size() - 1; i >= idx + 1; i--) vec[i] = vec[i-1];
  vec[idx] = val;
}

template<class T>
int tightvector<T>::getUpperBound(const T& val) {
  if (size() == 0) return size();
  int lo = 0, hi = size() - 1;
  while (lo < hi) {
    int mid = (lo + hi)/2;
    if (val < vec[mid]) hi = mid;
    else lo = mid + 1;
  }
  if (val < vec[lo]) return lo;
  return size();
}

template<class T>
int tightvector<T>::getLowerBound(const T& val) {
  if (size() == 0) return size();
  int lo = 0, hi = size() - 1;
  while (lo < hi) {
    int mid = (lo + hi)/2;
    if (vec[mid] < val) lo = mid + 1;
    else hi = mid;
  }
  if (!(vec[lo] < val)) return lo;
  return size();
}

template<class T>
tightvector<T>::operator vector<T>() const {
  vector<T> ret(size());
  for (int i = 0; i < size(); i++) ret[i] = (*this)[i];
  return ret;
}

/* --------- DecoratedProximitySearch --------- */
template<class T>
bool DecoratedProximitySearch<T>::checkSizes() {
  if (tags.size() != pointSize())
    MstUtils::error("number of points and custom tags does not agree", "DecoratedProximitySearch<T>::checkSizes()");
  return true;
}

#endif
