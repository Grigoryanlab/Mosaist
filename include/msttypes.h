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
using namespace std;

namespace MST {

// forward declarations
class Chain;
class Residue;
class Atom;
class Structure;

typedef double real;
typedef Structure System;                // for interchangability with MSL

class Structure {
  friend class Chain;

  public:
    Structure();
    Structure(string pdbFile, string options = "");
    Structure(Structure& S);
    ~Structure();
    
    void readPDB(string pdbFile, string options = "");
    void writePDB(string pdbFile, string options = "");
    void reset();

    int chainSize() const { return chains.size(); }
    int residueSize() { return numResidues; }
    int positionSize() { return residueSize(); }  // for interchangability with MSL
    int atomSize() { return numAtoms; }
    Chain* getChainByID(string id) { return (chainsByID.find(id) != chainsByID.end()) ? chainsByID[id] : NULL; }
    Chain* getChainBySegID(string id) { return (chainsBySegID.find(id) != chainsBySegID.end()) ? chainsBySegID[id] : NULL; }
    Chain& getChain(int i) { return (*this)[i]; }
    Chain& operator[](int i) const { return *(chains[i]); }
    vector<Atom*> getAtoms();
    vector<Residue*> getResidues();

    /* returns false if the chain name collides with existing chain names and no suitable single-letter
     * chain name was found as replacement OR if renaiming was not allowed. This could still mean that
     * a multi-character name is picked that is unique, but that's not technically correct for output,
     * so false will still be returned. Note that in the latter case, the segment ID will be renamed as
     * well to be the same multi-character name. So, although in the output chain names will repeat,
     * segment names will still be unique. If it fails to find an even multi-character name, errors. */
    bool appendChain(Chain* C, bool allowRename = true);

    /* makes a copy of the atom, then decides where it is supposed to go and inserts it
     * into the Structure, creating a new Chain and/or Residue as needed. */
    void addAtom(Atom* A);
    void addAtom(Atom& A) { addAtom(&A); }
    void addAtoms(vector<Atom*> atoms) { addAtoms(&atoms); }
    void addAtoms(vector<Atom*>* atoms);

  protected:
    void incrementNumAtoms(int delta = 1) { numAtoms += delta; }
    void incrementNumResidues(int delta = 1) { numResidues += delta; }
    void deletePointers();

  private:
    vector<Chain*> chains;
    string sourceFile;
    int numResidues, numAtoms;
    // NOTE: thse two maps are maintained for convenience and will not guarantee no collisions. That is,
    // If more than one chain use the same ID or segment ID, these maps will only store the last one added.
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
    Structure& getStructure() { return *parent; }
    string getID() { return cid; }
    string getSegID() { return sid; }

    /* convenience functoins, not efficient (linear search). If you need to do this a lot,
     * call getResidues() and construct your own data structure (e.g., a map<>) for fast lookups. */
    Residue* findResidue(string resname, int resnum);
    Residue* findResidue(string resname, int resnum, char icode);

    void setID(string _cid) { cid = _cid; }
    void setSegID(string _sid) { sid = _sid; }
    void appendResidue(Residue* R);

  protected:
    void setParent(Structure* p) { parent = p; } // will not itself update residue/atom counts in parent
    Structure* getParent() { return parent; }
    void incrementNumAtoms(int delta = 1);

  private:
    vector<Residue*> residues;
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
    Residue(Residue& R);
    Residue(string _resname, int _resnum, char _icode = ' ');
    ~Residue();

    int atomSize() { return atoms.size(); }
    vector<Atom*> getAtoms() { return atoms; }
    Atom& operator[](int i) { return *(atoms[i]); }
    Chain& getChain() { return *parent; }
    string getChainID(bool strict = true);
    string getName() { return resname; }
    int getNum() { return resnum; }
    char getIcode() { return icode; }
    bool isNamed(string& _name) { return (resname.compare(_name) == 0); }
    bool isNamed(const char* _name) { return (strcmp(resname.c_str(), _name) == 0); }
    Atom* findAtom(string _name, bool strict = true); // returns NULL if not found and if strict is false
    bool atomExists(string _name) { return (findAtom(_name, false) != NULL); } // mostly for interchangeability with MSL, better to use findAtom and check for NULL

    void setName(const char* _name) { resname = (string) _name; }
    void setName(string _name) { resname = _name; }
    char setIcode(char _icode) { icode = _icode; }
    void appendAtom(Atom* A);

    friend ostream & operator<<(ostream &_os, Residue& _res) {
      if (_res.getParent() != NULL) {
        _os << _res.getParent()->getID();
      }
      _os << _res.getNum() << " " << _res.getName();
    }

  protected:
    void setParent(Chain* _parent) { parent = _parent; } // will not itself update residue/atom counts in parents
    Chain* getParent() { return parent; }

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
    Atom(Atom& A);
    Atom(int _index, string _name, real _x, real _y, real _z, real _B, real _occ, bool _het, char _alt = ' ', Residue* _parent = NULL);
    ~Atom();

    Residue& getResidue() { return *parent; }
    real getX() const { return x; }
    real getY() const{ return y; }
    real getZ() const{ return z; }
    vector<real> getCoor() { vector<real> coor; coor.push_back(x); coor.push_back(y); coor.push_back(z); }
    real getB() { return B; }
    real getOcc() { return occ; }
    string getName() { return string(name); }
    char* getNameC() { return name; }
    bool isHetero() { return het; }
    int getIndex() { return index; }
    char getAlt() { return alt; }
    bool isNamed(const char* _name) { return (strcmp(name, _name) == 0); }
    bool isNamed(string& _name) { return isNamed(_name.c_str()); }
    int numAlternatives() { return (alternatives == NULL) ? 0 : alternatives->size(); }

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

    void addAlternative(real _x, real _y, real _z, real _B, real _occ, char _alt);

    string pdbLine() { return pdbLine((this->parent == NULL) ? 1 : this->parent->getNum(), index); }
    string pdbLine(int resIndex, int atomIndex);

    real distance(Atom& another);
    real distance(Atom* another) { return distance(*another); }
    real distance2(Atom& another);
    real distance2(Atom* another) { return distance2(*another); }

  protected:
    void setParent(Residue* _parent) { parent = _parent; } // will not itself update residue/atom counts in parents
    Residue* getParent() { return parent; }

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
 * use is absolutely optional, and none of the basic MST datastructures use this. On
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

    friend ostream & operator<<(ostream &_os, AtomPointerVector& _atoms) {
      for (int i = 0; i < _atoms.size(); i++) {
        _os << _atoms[i]->pdbLine() << endl;
      }
    }
};

}

/* The utilities class, with a bunch of useful static functions, is defined outside of the MST namespace because:
 * 1) it really represents a different beast, not an MST type
 * 2) some of its functions (like assert) are likely to clash with function names in other project
 */
class MstUtils {
  public:
    static void openFile (fstream& fs, string filename, ios_base::openmode mode, string from = "");
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

    template <class T>
    static string toString(T& obj);
};

template <class T>
string MstUtils::toString(T& obj) {
  stringstream ss;
  ss << obj;
  return ss.str();
}

#endif
