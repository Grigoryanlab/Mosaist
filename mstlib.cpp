#include "structlib.h"

/* --------- Structure --------- */
Structure::Structure() {
  numResidues = numAtoms = 0;
}

Structure::Structure(string pdbFile) {
  sourceFile = pdbFile;
  readPDB(pdbFile);
}

Structure::Structure(Structure& S) {
  sourceFile = S.sourceFile;
  numResidues = S.numResidues;
  numAtoms = S.numAtoms;
  for (int i = 0; i < S.chainSize(); i++) {
    chains.push_back(new Chain(S[i]));
  }
}

/* The assumption is that if a Structure is deleted, all
 * of its children objects are no longer needed and should
 * go away. If a user needs to hold on to them, they
 * should generate copies as needed via copy constructors.
 */
Structure::~Structure() {
  for (int i = 0; i < chains.size(); i++) delete(chains[i]);
}

vector<Atom*> Structure::getAtoms() {
  vector<Atoms*> vec;
  
  for (int i = 0; i < this->chainSize(); i++) {
    Chain& c = (*this)[i];
    for (int j = 0; j < c.residueSize(); j++) {
      Residue& r = c[j];
      vec.push_back(r.getAtoms());
    }
  }

  return vec;
}

/* --------- Chain --------- */
Chain::Chain() {
  numAtoms = 0;
  parent = NULL;
}

Chain::Chain(Chain& C) {
  numAtoms = C.numAtoms;
  parent = C.parent;
  for (int i = 0; i < C.residueSize(); i++) {
    residues.push_back(new Residue(C[i]));
  }
}

/* --------- Residue --------- */
Residue::Residue() {
  parent = NULL;
}

Residue::Residue(Residue& R) {
  parent = R.parent;
  for (int i = 0; i < R.atomSize(); i++) {
    atoms.push_back(R[i]);
  }
}

/* --------- Atom --------- */
Atom::Atom() {
  parent = NULL;
  het = false;
  name = NULL;
  setName("");
}

Atom::Atom(Atom& A) {
  x = A.x;
  y = A.y;
  z = A.z;
  occ = A.occ;
  B = A.B;
  setName(A.name);
  parent = A.parent;
  het = A.het;
  idnex = A.index;
}

void Atom::setName(char* _name) {
  if (name != NULL) delete[] name;
  name = new char[strlen(_name)+1];
  strcpy(name, _name);
}

void Atom::setName(string _name) {
  if (name != NULL) delete[] name;
  name = new char[_name.size()+1];
  strcpy(name, _name.c_str());
}
