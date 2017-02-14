#ifndef _MSTTYPES_H
#define _MSTTYPES_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <locale>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <map>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

namespace MST {

// forward declarations
class Chain;
class Residue;
class Atom;
class Structure;

typedef double real;
typedef Structure System;                // for interchangability with MSL
typedef int residueID[2];
typedef int atomID[3];

class Structure {
  friend class Chain;

  public:
    Structure();
    Structure(string pdbFile, string options = "");
    Structure(Structure& S);
    Structure(Chain& C);
    Structure(Residue& R);
    ~Structure();

    void readPDB(string pdbFile, string options = "");
    void writePDB(string pdbFile, string options = "");
    void writePDB(fstream& ofs, string options = "");
    void reset();

    int chainSize() const { return chains.size(); }
    int residueSize() { return numResidues; }
    int positionSize() { return residueSize(); }  // for interchangability with MSL
    int atomSize() { return numAtoms; }
    Chain* getChainByID(string id) { return (chainsByID.find(id) != chainsByID.end()) ? chainsByID[id] : NULL; }
    Chain* getChainBySegID(string id) { return (chainsBySegID.find(id) != chainsBySegID.end()) ? chainsBySegID[id] : NULL; }
    Chain& getChain(int i) { return (*this)[i]; }
    Residue& getResidue(int i);
    Chain& operator[](int i) const { return *(chains[i]); }
    vector<Atom*> getAtoms();
    vector<Residue*> getResidues();
    void setName(string _name) { name = _name; }
    void renumber(); // make residue numbering consequitive in each chain and atom index consequitive throughout

    /* ----- functions that grow/shrink structure ----- */
    /* returns false if the chain name collides with existing chain names and no suitable single-letter
     * chain name was found as replacement OR if renaiming was not allowed. This could still mean that
     * a multi-character name is picked that is unique, but that's not technically correct for output,
     * so false will still be returned. Note that in the latter case, the segment ID will be renamed as
     * well to be the same multi-character name. So, although in the output chain names will repeat,
     * segment names will still be unique. If it fails to find an even multi-character name, errors. */
    bool appendChain(Chain* C, bool allowRename = true);
    Chain* appendChain(string cid, bool allowRename = true);

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
    Chain(Chain& C);
    Chain(string chainID, string segID);
    ~Chain();

    int residueSize() { return residues.size(); }
    int positionSize() { return residueSize(); }  // for interchangability with MSL
    int atomSize() { return numAtoms; }
    Residue& operator[](int i) { return *(residues[i]); }
    Residue& getResidue(int i) { return (*this)[i]; }
    vector<Residue*> getResidues() { return residues; }
    vector<Atom*> getAtoms();
    string getID() { return cid; }
    string getSegID() { return sid; }
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
    Residue(Residue& R, bool copyAlt = true);
    Residue(string _resname, int _resnum, char _icode = ' ');
    ~Residue();

    int atomSize() { return atoms.size(); }
    vector<Atom*> getAtoms() { return atoms; }
    Atom& operator[](int i) { return *(atoms[i]); }
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

  public:
    Atom();
    Atom(Atom& A, bool copyAlt = true);
    Atom(int _index, string _name, real _x, real _y, real _z, real _B, real _occ, bool _het, char _alt = ' ', Residue* _parent = NULL);
    ~Atom();

    real getX() const { return x; }
    real getY() const{ return y; }
    real getZ() const{ return z; }
    real operator[](int i) const;
    vector<real> getCoor() { vector<real> coor; coor.push_back(x); coor.push_back(y); coor.push_back(z); return coor; }
    vector<real> getAltCoor(int altInd);
    real getB() { return B; }
    real getOcc() { return occ; }
    string getName() { return string(name); }
    char* getNameC() { return name; }
    bool isHetero() { return het; }
    int getIndex() { return index; }
    char getAlt() { return alt; }
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
    void setCoor(vector<real> xyz) { x = xyz[0]; y = xyz[1]; z = xyz[2]; }
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

    real distance(Atom& another);
    real distance(Atom* another) { return distance(*another); }
    real distance2(Atom& another);
    real distance2(Atom* another) { return distance2(*another); }

    friend ostream & operator<<(ostream &_os, Atom& _atom) {
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
    CartesianPoint cross(CartesianPoint other) const;
    real dot(CartesianPoint other) const;
    CartesianPoint getUnit() const { double L = norm(); return (*this/L); };

    // a few special access operations
    real getX() const { return (*this)[0]; }
    real getY() const { return (*this)[1]; }
    real getZ() const { return (*this)[2]; }

    real distance(CartesianPoint& another);
    real distance(CartesianPoint* another) { return distance(*another); }
    real distance2(CartesianPoint& another);
    real distance2(CartesianPoint* another) { return distance2(*another); }

    friend ostream & operator<<(ostream &_os, CartesianPoint& _p) {
      for (int i = 0; i < _p.size(); i++) {
        _os << _p[i];
        if (i != _p.size() - 1) _os << " ";
      }
      return _os;
    }
};

class CartesianGeometry {
  public:
    static real dihedralRadians(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4);
    static real dihedralRadians(const CartesianPoint * _p1, const CartesianPoint * _p2, const CartesianPoint * _p3, const CartesianPoint * _p4);
    static real dihedral(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4);
    static real dihedral(const CartesianPoint * _p1, const CartesianPoint * _p2, const CartesianPoint * _p3, const CartesianPoint * _p4);
};

class AtomPointerVector : public vector<Atom*> {
  public:
    // inherit a bunch of useful constructors from vector
    AtomPointerVector() : vector<Atom*>() { }
    AtomPointerVector(size_t sz, Atom* val) : vector<Atom*>(sz, val) { }
    AtomPointerVector(const AtomPointerVector& other) : vector<Atom*>(other) { }
    AtomPointerVector(const vector<Atom*>& other) : vector<Atom*>(other) { }

    CartesianPoint getGeometricCenter();
    real radiusOfGyration();
    void deletePointers();

    friend ostream & operator<<(ostream &_os, AtomPointerVector& _atoms) {
      for (int i = 0; i < _atoms.size(); i++) {
        _os << _atoms[i]->pdbLine() << endl;
      }
      return _os;
    }
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
    bool align(vector<Atom*> &_align, vector<Atom*> &_ref, vector<Atom*>& _moveable);

    // quickly calculate RMSD upon optimal superposition without generating the rotation matrix
    real bestRMSD(vector<Atom*> &_align, vector<Atom*> &_ref, bool* _suc = NULL, bool setTransRot = false);

    // in-place RMSD (no transformations)
    static real rmsd(vector<Atom*>& A, vector<Atom*>& B);
    static real rmsd(Structure& A, Structure& B);

 protected:
    // implemetation of Kabsch algoritm for optimal superposition
    bool Kabsch(vector<Atom*> &_align, vector<Atom*> &_ref, int mode);

 private:
    real _rmsd;
    real t[3];    // translation vector
    real u[3][3]; // rotation matrix

};

class ProximitySearch {
  public:
    ProximitySearch(real _xlo, real _ylo, real _zlo, real _xhi, real _yhi, real _zhi, int _N = 20);
    ProximitySearch(AtomPointerVector& _atoms, int _N, bool _addAtoms = true, vector<int>* tags = NULL, real pad = 0);
    ProximitySearch(AtomPointerVector& _atoms, real _characteristicDistance, bool _addAtoms = true, vector<int>* tags = NULL, real pad = 0);
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
    static void calculateExtent(AtomPointerVector& _atoms, real& _xlo, real& _ylo, real& _zlo, real& _xhi, real& _yhi, real& _zhi);
    static void calculateExtent(Structure& S, real& _xlo, real& _ylo, real& _zlo, real& _xhi, real& _yhi, real& _zhi);

    bool pointsWithin(CartesianPoint c, real dmin, real dmax, vector<int>* list = NULL, bool byTag = false);
    vector<int> getPointsWithin(CartesianPoint c, real dmin, real dmax, bool byTag = false);

    // Returns true if the grid of the current ProximitySearch object overlaps
    // that of the ProximitySearch specified by more than the padding given
    bool overlaps(ProximitySearch& other, real pad = 0);

  protected:
    void setBinWidths();
    void calculateExtent(AtomPointerVector& _atoms) { ProximitySearch::calculateExtent(_atoms, xlo, ylo, zlo, xhi, yhi, zhi); }

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

    vector<T> getPointsWithin(CartesianPoint c, real dmin, real dmax) {
      vector<int> inds = this->ProximitySearch::getPointsWithin(c, dmin, dmax, true);
      vector<T> ret(inds.size());
      for (int i = 0; i < inds.size(); i++) ret[i] = tags[inds[i]];
      return ret;
    }
    vector<int> getPointsWithinIndices(CartesianPoint c, real dmin, real dmax) {
      return this->ProximitySearch::getPointsWithin(c, dmin, dmax, true);
    }

  private:
    vector<T> tags;
};

}

/* Utilities class, with a bunch of useful static functions, is defined outside of the MST namespace because:
 * 1) it really represents a different beast, not an MST type
 * 2) some of its functions (like assert) are likely to clash with function names in other project
 */
class MstUtils {
  public:
    static void openFile (fstream& fs, string filename, ios_base::openmode mode, string from = "");
    static void fileToArray(string _filename, vector<string>& lines); // reads lines from the file and appends them to the given vector
    static vector<string> fileToArray(string _filename) { vector<string> lines; fileToArray(_filename, lines); return lines; }
    static FILE* openFileC (const char* filename, const char* mode, string from = "");
    static string trim(string str);
    static void warn(string message, string from = "");
    static void error(string message, string from = "", int code = -1);
    static void assert(bool condition, string message = "error: assertion failed", string from = "", int exitCode = -1);
    static string uc(string& str);                        // returns an upper-case copy of the input string
    static string wrapText(string message, int width, int leftSkip = 0, int startingOffset = 0);
    static char* copyStringC(const char* str);
    static int toInt(string num, bool strict = true);
    static MST::real toReal(string num, bool strict = true);
    static MST::real mod(MST::real num, MST::real den);
    static MST::real sign(MST::real val) { return (val > 0) ? 1.0 : ((val < 0) ? -1.0 : 0.0); }
    static string pathBase(string fn); // gets the base name of the path (removes the extension)
    static bool fileExists(const char *filename);
    static bool isDir(const char *filename);

    static string readNullTerminatedString(fstream& ifs);

    // returns a random number in the range [lower, upper]
    static int randInt(int lower, int upper) { return rand() % (upper - lower + 1) + lower; }
    // returns a random number in the range [0, upper) (convenient for generating random array subscripts)
    static int randInt(int upper) { return randInt(0, upper - 1); }

    template <class T>
    static string toString(T obj) { return toString(&obj); }
    template <class T>
    static string toString(T* obj);
    template <class T>
    static vector<int> sortIndices(vector<T>& vec, bool descending = false);
    template <class T1, class T2>
    static vector<T1> keys(map<T1, T2>& _map);
    template <class T>
    static string vecToString(vector<T>& vec, string del = " ");
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
string MstUtils::vecToString(vector<T>& vec, string del) {
  string str;
  for (int i = 0; i < vec.size(); i++) {
    str += MstUtils::toString(vec[i]);
    if (i != vec.size() - 1) str += del;
  }
  return str;
}

#endif
