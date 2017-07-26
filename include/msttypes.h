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
#undef assert

using namespace std;

namespace MST {

// forward declarations
class Chain;
class Residue;
class Atom;
class Structure;
class CartesianPoint;

typedef double real;
typedef Structure System;                // for interchangability with MSL

class Structure {
  friend class Chain;

  public:
    Structure();
    Structure(string pdbFile, string options = "");
    Structure(const Structure& S);
    Structure(Chain& C);
    Structure(Residue& R);
    Structure(vector<Atom*>& atoms);
    ~Structure();

    void readPDB(string pdbFile, string options = "");
    void writePDB(string pdbFile, string options = "");
    void writePDB(fstream& ofs, string options = "");
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
    vector<Residue*> getResidues();
    void setName(string _name) { name = _name; }
    string getName() { return name; }
    void renumber(); // make residue numbering consequitive in each chain and atom index consequitive throughout

    // looks at the length of the peptide bond between adjacent residues to figure out where chains break
    void reassignChainsByConnectivity(Structure& dest, real maxPeptideBond = 2.0);
    Structure reassignChainsByConnectivity(real maxPeptideBond = 2.0);

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

    /* makes a copy of the atom, then decides where it is supposed to go and inserts it
     * into the Structure, creating a new Chain and/or Residue as needed. */
    void addAtom(Atom* A);
    void addAtom(Atom& A) { addAtom(&A); }
    void addAtoms(vector<Atom*> atoms) { addAtoms(&atoms); }
    void addAtoms(vector<Atom*>* atoms);

    /* makes a copy of the residue, then decides where it is supposed to go and
     * inserts it into the Structure, creating a new Chain as needed. */
    void addResidue(Residue* res);
    /* ----- functions that grow/shrink structure ----- */

    int getResidueIndex(Residue* res);

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
    Chain(string chainID, string segID);
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
    Structure* getParent() { return parent; }
    Structure* getStructure() { return getParent(); }
    int getResidueIndex(Residue* res);

    /* convenience functoins, not efficient (linear search). If you need to do this a lot,
     * call getResidues() and construct your own data structure (e.g., a map<>) for fast lookups. */
    Residue* findResidue(string resname, int resnum);
    Residue* findResidue(string resname, int resnum, char icode);

    void setID(string _cid) { cid = _cid; }
    void setSegID(string _sid) { sid = _sid; }

    /* ----- functions that grow/shrink structure ----- */
    void appendResidue(Residue* R);
    void insertResidue(Residue* R, int index); // insert the Residue in such a way that it ends up being at index i
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
    Atom& getAtom(int i) { return *(atoms[i]); }
    Chain* getChain() { return parent; }
    string getChainID(bool strict = true);
    string getName() { return resname; }
    int getNum() { return resnum; }
    char getIcode() { return icode; }
    bool isNamed(string& _name) { return (resname.compare(_name) == 0); }
    bool isNamed(const char* _name) { return (strcmp(resname.c_str(), _name) == 0); }
    Atom* findAtom(string _name, bool strict = true); // returns NULL if not found and if strict is false
    bool atomExists(string _name) { return (findAtom(_name, false) != NULL); } // mostly for interchangeability with MSL, better to use findAtom and check for NULL
    Chain* getParent() { return parent; }
    Structure* getStructure() { return (parent == NULL) ? NULL : parent->getParent(); }

    void setName(const char* _name) { resname = (string) _name; }
    void setName(string _name) { resname = _name; }
    void setIcode(char _icode) { icode = _icode; }
    void setNum(int num) { resnum = num; }
    void copyAtoms(Residue& R, bool copyAlt = true);

    /* for all atoms in this residue, overwrite the main coordinate set with the
     * coordinate set from the alternative with the specified index. */
    void makeAlternativeMain(int altInd);

    /* ----- functions that grow/shrink structure ----- */
    void appendAtom(Atom* A);
    void appendAtoms(vector<Atom*>& A);
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
    void replaceAtoms(vector<Atom*>& newAtoms, vector<int>* oldAtoms = NULL);
    /* ----- functions that grow/shrink structure ----- */

    Residue* previousResidue();
    Residue* nextResidue();
    Residue* iPlusDelta(int off);
    real getPhi(bool strict = true);
    real getPsi(bool strict = true);
    real getOmega(bool strict = true);

    int getResidueIndex();

    static const real badDihedral; // value that signals a dihedral angle that could not be computed for some reason

    friend ostream & operator<<(ostream &_os, Residue& _res) {
      if (_res.getParent() != NULL) {
        _os << _res.getParent()->getID() << ",";
      }
      _os << _res.getNum() << " " << _res.getName();
      return _os;
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

  public:
    Atom();
    Atom(const Atom& A, bool copyAlt = true);
    Atom(int _index, string _name, real _x, real _y, real _z, real _B, real _occ, bool _het, char _alt = ' ', Residue* _parent = NULL);
    ~Atom();

    real getX() const { return x; }
    real getY() const{ return y; }
    real getZ() const{ return z; }
    real& operator[](int i);
    CartesianPoint getCoor() const;
    CartesianPoint getAltCoor(int altInd) const;
    real getB() { return B; }
    real getOcc() { return occ; }
    string getName() const { return string(name); }
    char* getNameC() { return name; }
    bool isHetero() const { return het; }
    int getIndex() const { return index; }
    char getAlt() const { return alt; }
    bool isNamed(const char* _name) { return (strcmp(name, _name) == 0); }
    bool isNamed(string _name) { return isNamed(_name.c_str()); }
    int numAlternatives() { return (alternatives == NULL) ? 0 : alternatives->size(); }
    Residue* getParent() { return parent; }
    Residue* getResidue() { return parent; }
    Chain* getChain() { return (parent == NULL) ? NULL : parent->getParent(); }
    Structure* getStructure() { Chain* chain = getChain(); return (chain == NULL) ? NULL : chain->getParent(); }

    void setName(const char* _name);
    void setName(string _name);
    void setX(real _x) { x = _x; }
    void setY(real _y) { y = _y; }
    void setZ(real _z) { z = _z; }
    void setCoor(real _x, real _y, real _z) { x = _x; y = _y; z = _z; }
    void setCoor(const CartesianPoint& xyz);
    void setOcc(real _occ) { occ = _occ; }
    void setB(real _B) { B = _B; }
    void seetHet(bool _het) { het = _het; }
    void setIndex(int _index) { index = _index; }

    /* make the alternative with the specified index the main one, making the current
     * main position the alternative with the specified index. Calling this twice with
     * the same index will return things back to the way they were originally. */
    void swapWithAlternative(int altInd);

    /* overwrite the main coordinate set with the coordinate set from the alternative
     * with the specified index. */
    void makeAlternativeMain(int altInd);

    void addAlternative(real _x, real _y, real _z, real _B, real _occ, char _alt = ' ');

    string pdbLine() { return pdbLine((this->parent == NULL) ? 1 : this->parent->getNum(), index); }
    string pdbLine(int resIndex, int atomIndex);

    real distance(const Atom& another) const;
    real distance(const Atom* another) const { return distance(*another); }
    real distance2(const Atom& another) const;
    real distance2(const Atom* another) const { return distance2(*another); }
    real angle(const Atom& A, const Atom& B, bool radians = false) const;
    real angle(const Atom* A, const Atom* B, bool radians = false) const;
    real dihedral(const Atom& A, const Atom& B, const Atom& C, bool radians = false) const;
    real dihedral(const Atom* A, const Atom* B, const Atom* C, bool radians = false) const;

    /* Sets the coordinates of the atom based on internal coordinates relative
     * to three other atoms: thA, anA, diA, A (this atom). The internal
     * coordinates are distance diA-A (di), angle anA-diA-A (an), and dihedral
     * angle thA - anA - diA - A (th). */
    bool build(const Atom& diA, const Atom& anA, const Atom& thA, real di, real an, real th, bool radians = false);
    bool build(const Atom* diA, const Atom* anA, const Atom* thA, real di, real an, real th, bool radians = false) {
      return build(*diA, *anA, *thA, di, an, th, radians);
    }

    friend ostream & operator<<(ostream &_os, const Atom& _atom) {
      _os << _atom.getName() << _atom.getAlt() << " " << _atom.index << " " << (_atom.isHetero() ? "HETERO" : "");
      _os << _atom.x << " " << _atom.y << " " << _atom.z << " : " << _atom.occ << " " << _atom.B;
      return _os;
    }

  protected:
    void setParent(Residue* _parent) { parent = _parent; } // will not itself update residue/atom counts in parents

  private:
    real x, y, z, occ, B;
    char *name, alt;
    Residue* parent;
    bool het;
    int index;

    // data structure for storing information about alternative atom locations
    class altInfo {
      public:
        altInfo() { x = y = z = occ = B = 0; alt = ' '; }
        altInfo(const altInfo& A) { x = A.x; y = A.y; z = A.z; B = A.B; occ = A.occ; alt = A.alt; }
        altInfo(real _x, real _y, real _z, real _occ, real _B, char _alt) { x = _x; y = _y; z = _z; B = _B; occ = _occ; alt = _alt; }
        real x, y, z, occ, B;
        char alt;
    };
    vector<altInfo>* alternatives; /* since this is a pointer, and will be NULL for most atoms, it's fine
                                    * to use vector here in terms of memory, but very convenient for use */
};

/* The following several classes look and feel like MSL classes, BUT (importantly) their
 * use is absolutely optional, and none of the basic MST datastructures use them. On
 * the other hand, they can be created from those basic types for a similar use as
 * in MSL when needed (and only when needed). For example, CartesianPoint knows how
 * to construct itself from atom, both via a constructor and assignment operator.
 * Similarly, AtomPointerVector knows how to construct itself from vector<Atom*>. */
class CartesianPoint : public vector<real> {
  /* this class it no limited to 3D vectors, though some of the functions will only
   * work with 3D vectors. The intention is to make it general, such that if a 3D
   * vector is required (or another dimension), appropriate assertions are made. */
  public:
    // inherit a bunch of useful constructors from vector
    CartesianPoint() : vector<real>() { }
    CartesianPoint(size_t sz) : vector<real>(sz) { }
    CartesianPoint(size_t sz, real val) : vector<real>(sz, val) { }
    CartesianPoint(const CartesianPoint& other) : vector<real>(other) { }
    CartesianPoint(const vector<real>& other) : vector<real>(other) { }
    CartesianPoint(real x, real y, real z) : vector<real>(3, 0) { (*this)[0] = x; (*this)[1] = y; (*this)[2] = z; }

    CartesianPoint(const Atom& A);
    CartesianPoint(const Atom* A) : CartesianPoint(*A) {}
    CartesianPoint& operator+=(const CartesianPoint& rhs);
    CartesianPoint& operator-=(const CartesianPoint& rhs);
    CartesianPoint& operator/=(const real& s);
    CartesianPoint& operator*=(const real& s);
    const CartesianPoint operator+(const CartesianPoint &other) const;
    const CartesianPoint operator-(const CartesianPoint &other) const;
    const CartesianPoint operator*(const real& s) const;
    const CartesianPoint operator/(const real& s) const;
    const CartesianPoint operator-() const;
    CartesianPoint& operator=(const Atom& A);
    const double operator*(const CartesianPoint& other) const { return this->dot(other); }

    real norm() const;
    real mean() const;
    real median() const;
    real sum() const;
    CartesianPoint cross(CartesianPoint other) const;
    real dot(CartesianPoint other) const;
    CartesianPoint getUnit() const { double L = norm(); return (*this/L); };

    // a few special access operations
    real getX() const { return (*this)[0]; }
    real getY() const { return (*this)[1]; }
    real getZ() const { return (*this)[2]; }

    real distance(const CartesianPoint& another) const;
    real distance(const CartesianPoint* another) const { return distance(*another); }
    real distance2(const CartesianPoint& another) const;
    real distance2(const CartesianPoint* another) const { return distance2(*another); }

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
    static real dihedral(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4, bool radians = false);
    static real dihedral(const CartesianPoint * _p1, const CartesianPoint * _p2, const CartesianPoint * _p3, const CartesianPoint * _p4, bool radians = false);
    static real angle(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, bool radians = false);
    static real angle(const CartesianPoint * _p1, const CartesianPoint * _p2, const CartesianPoint * _p3, bool radians = false);
    static real angleDiff(real A, real B, bool radians = false);
    static real angleDiffCCW(real A, real B, bool radians = false);

    /* Given a list of angles defined on the full circle, it is actually not clear
     * where the min and max are (since angles are defined on a circle). Exampe:
     * we are given -179 and +179. Is the range essentially the entire circle,
     * from -179 to +179? Or is it just a two-degree arc from 179 to 181 (aka
     * -179)? Both are correct. This function finds the min/max bounds that
     * minimize the arc length corresponding to the range. This is often
     * convenient. The functoin will try each angle value as the min, find its
     * corresponding max, and pick the pair that minimizes the arc length. */
    static pair<real, real> angleRange(const vector<real>& angles, bool radians = false);

    /* Just like in the above, the average of a list of angles is also not
     * necessarily uniquely defined. This functon defines it by first finding
     * the angular range with the smallest arc length, and then finding the mean
     * with respect to that range (i.e., via the average deviation from the
     * minimal value). */
    static real angleMean(const vector<real>& angles, bool radians = false);
};

class AtomPointerVector : public vector<Atom*> {
  public:
    // inherit a bunch of useful constructors from vector
    AtomPointerVector() : vector<Atom*>() { }
    AtomPointerVector(size_t sz, Atom* val) : vector<Atom*>(sz, val) { }
    AtomPointerVector(const AtomPointerVector& other) : vector<Atom*>(other) { }
    AtomPointerVector(const vector<Atom*>& other) : vector<Atom*>(other) { }

    CartesianPoint getGeometricCenter();
    void center();
    real radiusOfGyration();
    // computes the smallest radius of a sphere centered at the centroid of the
    // atom vector, which encloses the set of all atoms.
    real boundingSphereRadiusCent();
    void deletePointers();

    AtomPointerVector clone();
    void clone(AtomPointerVector& into);
    AtomPointerVector subvector(int beg, int end);

    friend ostream & operator<<(ostream &_os, const AtomPointerVector& _atoms) {
      for (int i = 0; i < _atoms.size(); i++) {
        _os << _atoms[i]->pdbLine() << endl;
      }
      return _os;
    }
};

class expressionTree {
  public:
    enum selProperty { RESID = 1, RESNAME, CHAIN, SEGID, NAME, AROUND }; // selectable properties
    enum logicalOp { AND = 1, OR, NOT, BYRES, BYCHAIN, IS };             // logical operators

    expressionTree(logicalOp _op = logicalOp::IS) { op = _op; }
    ~expressionTree() {
      for (int i = 0; i < children.size(); i++) delete children[i];
    }

    void setLogicalOperator(logicalOp _op) { op = _op; }
    void setProperty(selProperty _type) { type = _type; }
    void setNum(int _num) { num = _num; }
    void setString(string _str) { str = _str; }
    void addChild(expressionTree* subtree) { children.push_back(subtree); }
    logicalOp getLogicalOperator() { return op; }
    selProperty getProperty() { return type; }
    int getNum() { return num; }
    string getString() { return str; }
    int numChildren() { return children.size(); }
    expressionTree* getChild(int i) { return children[i]; }

  private:
    logicalOp op;
    selProperty type;
    int num;
    string str;
    vector<expressionTree*> children;
};

class selector {
  public:
    selector(Structure& S);
    AtomPointerVector select(string selStr);
    vector<Residue*> selectRes(string selStr);
    void select(expressionTree* tree, AtomPointerVector& sel);
    expressionTree* buildExpressionTree(string selStr);
    AtomPointerVector byRes(AtomPointerVector& selAtoms);
    AtomPointerVector byChain(AtomPointerVector& selAtoms);
    AtomPointerVector invert(AtomPointerVector& selAtoms);
    AtomPointerVector intersect(AtomPointerVector& selA, AtomPointerVector& selB);
    AtomPointerVector combine(AtomPointerVector& selA, AtomPointerVector& selB);

  private:
    string getNextSelectionToken(string& selStr);
    vector<Atom*> atoms;
    vector<Residue*> residues;
    vector<Chain*> chains;
};

class RMSDCalculator {
 public:
    RMSDCalculator() {}
    ~RMSDCalculator() {}

    // getters
    real lastRMSD() { return _rmsd; }
    vector<real> lastTranslation();
    vector<vector<real> > lastRotation();

    // calculate optimal superposition and the resulting RMSD, applying transformation to given atoms
    bool align(const vector<Atom*> &_align, const vector<Atom*> &_ref, vector<Atom*>& _moveable);
    bool align(const vector<Atom*> &_align, const vector<Atom*> &_ref, Structure& _moveable);

    // quickly calculate RMSD upon optimal superposition without generating the rotation matrix
    real bestRMSD(const vector<Atom*> &_align, const vector<Atom*> &_ref, bool setTransRot = false, bool* _suc = NULL);

    // in-place RMSD (no transformations)
    static real rmsd(const vector<Atom*>& A, const vector<Atom*>& B);
    static real rmsd(const Structure& A, const Structure& B);

    // Craig's cutoff function
    static real rmsdCutoff(const vector<int>& L, real rmsdMax = 1.1, real L0 = 15);
    static real rmsdCutoff(const Structure& S, real rmsdMax = 1.1, real L0 = 15);

 protected:
    // implemetation of Kabsch algoritm for optimal superposition
    bool Kabsch(const vector<Atom*> &_align, const vector<Atom*> &_ref, int mode);

 private:
    real _rmsd;
    real t[3];    // translation vector
    real u[3][3]; // rotation matrix

};

class ProximitySearch {
  public:
    ProximitySearch(real _xlo, real _ylo, real _zlo, real _xhi, real _yhi, real _zhi, int _N = 20);
    ProximitySearch(const AtomPointerVector& _atoms, int _N, bool _addAtoms = true, vector<int>* tags = NULL, real pad = 0);
    ProximitySearch(const AtomPointerVector& _atoms, real _characteristicDistance, bool _addAtoms = true, vector<int>* tags = NULL, real pad = 0);
    ~ProximitySearch();

    real getXLow() { return xlo; }
    real getYLow() { return ylo; }
    real getZLow() { return zlo; }
    real getXHigh() { return xhi; }
    real getYHigh() { return yhi; }
    real getZHigh() { return zhi; }
    int pointSize() { return pointList.size(); }
    CartesianPoint& getPoint(int i) { return *(pointList[i]); }
    int getPointTag(int i) { return pointTags[i]; }
    real distance(int i, int j) { return pointList[i]->distance(pointList[j]); }

    void reinitBuckets(int _N);
    void addPoint(CartesianPoint _p, int tag);
    void addAtoms(AtomPointerVector& apv, vector<int>* tags = NULL);
    bool isPointWithinGrid(CartesianPoint _p);
    void pointBucket(CartesianPoint* p, int* i, int* j, int* k) { pointBucket(p->getX(), p->getY(), p->getZ(), i, j, k); }
    void pointBucket(CartesianPoint p, int* i, int* j, int* k) { pointBucket(p.getX(), p.getY(), p.getZ(), i, j, k); }
    void pointBucket(real px, real py, real pz, int* i, int* j, int* k);
    void limitIndex(int *ind);
    real gridSpacingX() { return (xhi - xlo)/N; }
    real gridSpacingY() { return (yhi - ylo)/N; }
    real gridSpacingZ() { return (zhi - zlo)/N; }
    static void calculateExtent(const AtomPointerVector& _atoms, real& _xlo, real& _ylo, real& _zlo, real& _xhi, real& _yhi, real& _zhi);
    static void calculateExtent(const Structure& S, real& _xlo, real& _ylo, real& _zlo, real& _xhi, real& _yhi, real& _zhi);

    bool pointsWithin(const CartesianPoint& c, real dmin, real dmax, vector<int>* list = NULL, bool byTag = false);
    vector<int> getPointsWithin(const CartesianPoint& c, real dmin, real dmax, bool byTag = false);
    int numPointsWithin(const CartesianPoint& c, real dmin, real dmax) {
      vector<int> closeOnes; pointsWithin(c, dmin, dmax, &closeOnes); return closeOnes.size();
    }

    // Returns true if the grid of the current ProximitySearch object overlaps
    // that of the ProximitySearch specified by more than the padding given
    bool overlaps(ProximitySearch& other, real pad = 0);

  protected:
    void setBinWidths();
    void calculateExtent(const AtomPointerVector& _atoms) { ProximitySearch::calculateExtent(_atoms, xlo, ylo, zlo, xhi, yhi, zhi); }

  private:
    int N; // dimension of bucket list is N x N x N

    real xlo, ylo, zlo, xhi, yhi, zhi, xbw, ybw, zbw; // extents of coordinates

    /* Each bucket is a vector of point indices (zero-initiated). These
     * indices are into two vectors: a vector of 3D points (CartesianPoint
     * pointers) and a vector of point tags. In this way, it is easy to go
     * from a bucket into points as well as from a point index into its point
     * or tag. One can also go from a point to its bucket location, by doing
     * a simple computations on the coordinates of the point via poitBucket(). */
    vector<vector<vector<vector<int> > > > buckets;
    vector<CartesianPoint*> pointList;
    vector<int> pointTags;
};

template<class T>
class DecoratedProximitySearch : public ProximitySearch {
  public:
    DecoratedProximitySearch(real _xlo, real _ylo, real _zlo, real _xhi, real _yhi, real _zhi, int _N = 20) :
      ProximitySearch(_xlo, _ylo, _zlo, _xhi, _yhi, _zhi, _N) {}
    DecoratedProximitySearch(AtomPointerVector& _atoms, int _N, vector<T>& _tags, real pad = 0) :
      ProximitySearch(_atoms, _N, true, NULL, pad) {
      tags = _tags;
    }
    DecoratedProximitySearch(AtomPointerVector& _atoms, int _N, real pad = 0) : ProximitySearch(_atoms, _N, false, NULL, pad) {}
    DecoratedProximitySearch(AtomPointerVector& _atoms, real _characteristicDistance, vector<T>& _tags, real pad = 0) :
      ProximitySearch(_atoms, _characteristicDistance, true, NULL, pad) {
      tags = _tags;
    }
    DecoratedProximitySearch(AtomPointerVector& _atoms, real _characteristicDistance, real pad = 0) :
      ProximitySearch(_atoms, _characteristicDistance, false, NULL, pad) { }

    T getPointTag(int i) { return tags[this->ProximitySearch::getPointTag(i)]; }
    void addPoint(CartesianPoint _p, T tag) {
      this->ProximitySearch::addPoint(_p, tags.size());
      tags.push_back(tag);
    }

    vector<T> getPointsWithin(const CartesianPoint& c, real dmin, real dmax) {
      vector<int> inds = this->ProximitySearch::getPointsWithin(c, dmin, dmax, true);
      vector<T> ret(inds.size());
      for (int i = 0; i < inds.size(); i++) ret[i] = tags[inds[i]];
      return ret;
    }
    vector<int> getPointsWithinIndices(const CartesianPoint& c, real dmin, real dmax) {
      return this->ProximitySearch::getPointsWithin(c, dmin, dmax, true);
    }

  private:
    vector<T> tags;
};


class Clusterer {
  public:
    /* Will greedy cluster the given set of units (must all have the same number of atoms),
     * using the given RMSD cutoff, while making sure that no more than Nmax x Nmax RMSD
     * computations are done per iteration. So, if the number of units is below Nmax, a
     * normal greedy clustering will be performed. However, if above Nmax, then will try
     * be further greedy and find best centroids without ever doing all-by-all comparisons. */
    vector<vector<int> > greedyCluster(const vector<vector<Atom*> >& units, real rmsdCut, int Nmax = 10000);


  protected:
    // these functions are protected because they assume that the cache of pre-
    // computed RMSDs is in a good state (so don't want external calls)
    vector<vector<int> > greedyClusterBruteForce(const vector<vector<Atom*> >& units, set<int> remIndices, real rmsdCut, int nClusts = -1);
    vector<int> elementsWithin(const vector<vector<Atom*> >& units, set<int>& remIndices, int from, real rmsdCut);
    set<int> randomSubsample(set<int>& indices, int N);

  private:
    // won't cache for now (need careful memory management)
    map<int, map<int, real > > coputedRMSDs;
};

}

/* Utilities class, with a bunch of useful static functions, is defined outside of the MST namespace because:
 * 1) it really represents a different beast, not an MST type
 * 2) some of its functions (like assert) are likely to clash with function names in other project
 */
class MstUtils {
  public:
    static void openFile (fstream& fs, string filename, ios_base::openmode mode = ios_base::in, string from = "");
    static void fileToArray(string _filename, vector<string>& lines); // reads lines from the file and appends them to the given vector
    static vector<string> fileToArray(string _filename) { vector<string> lines; fileToArray(_filename, lines); return lines; }
    static FILE* openFileC (const char* filename, const char* mode, string from = "");
    static string trim(string str, string delimiters = " \t\n\v\f\r");
    static void warn(string message, string from = "");
    static void error(string message, string from = "", int code = -1);
    static void assert(bool condition, string message = "error: assertion failed", string from = "", int exitCode = -1);
    static string uc(const string& str);                        // returns an upper-case copy of the input string
    static string lc(const string& str);                        // returns an lower-case copy of the input string
    static bool stringsEqual(const string& A, const string& B, bool caseInsensitive = true);
    static string wrapText(string message, int width, int leftSkip = 0, int startingOffset = 0);
    static char* copyStringC(const char* str);
    static int toInt(string num, bool strict = true);
    static bool isInt(string num);
    static MST::real toReal(string num, bool strict = true);
    static bool isReal(string num);
    static MST::real mod(MST::real num, MST::real den);
    static MST::real sign(MST::real val) { return (val > 0) ? 1.0 : ((val < 0) ? -1.0 : 0.0); }
    static string nextToken(string& str, string delimiters = " ", bool skipTrailingDelims = true);
    static vector<string> split(const string& str, string delimiters = " ", bool skipTrailingDelims = true);
    static string readNullTerminatedString(fstream& ifs);
    static string getDate();

    // returns a random number in the range [lower, upper]
    static int randInt(int lower, int upper) { return rand() % (upper - lower + 1) + lower; }
    // returns a random number in the range [0, upper) (convenient for generating random array subscripts)
    static int randInt(int upper) { return randInt(0, upper - 1); }
    // random number in the unit range 0 and 1
    static MST::real randUnit() { return ((MST::real) rand() / RAND_MAX); }

    template <class T>
    static string toString(T obj) { return toString(&obj); }
    template <class T>
    static string toString(T* obj);
    template <class T>
    static vector<int> sortIndices(vector<T>& vec, bool descending = false);
    template <class T1, class T2>
    static vector<T1> keys(map<T1, T2>& _map);
    template <class T>
    static string vecToString(const vector<T>& vec, string del = " ");
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

    // randomly shuffle the array in place using the Fisher-Yates shuffle
    template <class T>
    static void shuffle(vector<T>& vec);
};

template <class T>
string MstUtils::toString(T* obj) {
  stringstream ss;
  ss << *obj;
  return ss.str();
}

template <class T>
vector<int> MstUtils::sortIndices(vector<T>& vec, bool descending) {
  vector<int> sortedIndices(vec.size(), 0);
  for (int i = 0; i < vec.size(); i++) sortedIndices[i] = i;
  if (descending) {
    sort(sortedIndices.begin(), sortedIndices.end(), [&vec](size_t i1, size_t i2) {return vec[i2] < vec[i1];});
  } else {
    sort(sortedIndices.begin(), sortedIndices.end(), [&vec](size_t i1, size_t i2) {return vec[i1] < vec[i2];});
  }
  return sortedIndices;
}

template <class T1, class T2>
vector<T1> MstUtils::keys(map<T1, T2>& _map) {
  vector<T1> K(_map.size());
  int k = 0;
  for (typename map<T1, T2>::iterator it = _map.begin(); it != _map.end(); ++it) {
    K[k] = it->first;
  }
  return K;
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
  MstUtils::assert(vec.size() > 0, "empty vector passed!", "MstUtils::min(vector<T>&)");
  if (beg < 0) beg = 0;
  if (end < 0) end = vec.size()-1;
  MstUtils::assert(beg <= end, "beg = " + MstUtils::toString(beg) + " and end = " + MstUtils::toString(end), "MstUtils::min(vector<T>&)");
  MstUtils::assert(end < vec.size(), "end index out of bounds, " + MstUtils::toString(end), "MstUtils::min(vector<T>&)");

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
  MstUtils::assert(vec.size() > 0, "empty vector passed!", "MstUtils::max(vector<T>&)");
  if (beg < 0) beg = 0;
  if (end < 0) end = vec.size()-1;
  MstUtils::assert(beg <= end, "beg = " + MstUtils::toString(beg) + " and end = " + MstUtils::toString(end), "MstUtils::max(vector<T>&)");
  MstUtils::assert(end < vec.size(), "end index out of bounds, " + MstUtils::toString(end), "MstUtils::max(vector<T>&)");

  T mv = vec[beg];
  if (maxIndex != NULL) *maxIndex = beg;
  for (int i = beg; i <= end; i++) {
    if (vec[i] > mv) mv = vec[i];
    if (maxIndex != NULL) *maxIndex = i;
  }
  return mv;
}

template <class T>
bool MstUtils::closeEnough(const T& a, const T& b, const T& epsilon) {
  return (a - b > -epsilon) && (a - b < epsilon);
}

template <class T>
void MstUtils::shuffle(vector<T>& vec) {
  int i, j, tmp;
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

#endif
