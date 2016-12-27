#ifndef _STRUCTLIB_H
#define _STRUCTLIB_H

// forward declarations
class Chain;
class Residue;
class Atom;

typedef double real;

class Structure {
  public:
    Structure();
    Structure(string pdbFile);
    Structure(Structure& S);
    ~Structure();
    
    void readPDB(string pdbFile); // TODO (sets numResidues and numAtoms)
    void writePDB(string pdbFile); // TODO

    int chainSize() { return chains.size(); }
    int residueSize() { return numResidues; }
    int atomSize() { return numAtoms; }

    Chain& operator[](int i) { return *(chains[i]); }
    vector<Atom*> getAtoms();
    
  private:
    vector<Chain*> chains;
    string sourceFile;
    int numResidues, numAtoms;
};

class Chain {
  public:
    Chain();
    Chain(Chain& C);
  
    int residueSize() { return residues.size(); }
    int atomSize() { return numAtoms; }
    Residue& operator[](int i) { return *(residues[i]); }
    Structure& getStructure() { return *parent; }

  private:
    vector<Residue*> residues;
    Structure* parent;
    int numAtoms;
};

class Residue {
  public:
    Residue();
    Residue(Residue& R);

    int atomSize() { return atoms.size(); }
    vector<Atoms*> getAtoms() { return atoms; }
    Atom& operator[](int i) { return *(atoms[i]); }
    Chain& getChain() { return *parent; }

  private:
    vector<Atom*> atoms;
    Chain* parent;
};

class Atom {
  public:
    Atom();
    Atom(Atom& A);

    Residue& getResidue() { return *parent; }
    real getX() { return x; }
    real getY() { return y; }
    real getZ() { return z; }
    vector<real> getCoor() { vector<real> coor; coor.push_back(x); coor.push_back(y); coor.push_back(z); }
    real getB() { return B; }
    real getOcc() { return occ; }
    string getName() { return string(name); }
    string getNamePtr() { return name; }
    bool isHetero() { return het; }
    int getIndex() { return index; }

    void setName(char* _name);
    void setName(string _name);
    void setX(real _x) { x = _x; }
    void setY(real _y) { y = _y; }
    void setY(real _z) { z = _z; }
    void setCoor(real _x, real _y, real _z) { x = _x; y = _y; z = _z; }
    void setCoor(vector<real> xyz) { x = xyz[0]; y = xyz[1]; z = xyz[2]; }
    void setOcc(real _occ) { occ = _occ; }
    void setB(real _B) { B = _B; }
    void seetHet(bool _het) { het = _het; }
    void setIndex(int _index) { index = _index; }

  private:
    real x, y, z, occ, B;
    char* name;
    Residue* parent;
    bool het;
    int index, 
};
#endif
