#include "msttypes.h"

using namespace MST;

/* --------- Structure --------- */
Structure::Structure() {
  numResidues = numAtoms = 0;
  initExtraMaps();
}

Structure::Structure(string pdbFile, string options) {
  sourceFile = pdbFile;
  numResidues = numAtoms = 0;
  initExtraMaps();
  readPDB(pdbFile, options);
}

Structure::Structure(Structure& S) {
  sourceFile = S.sourceFile;
  numResidues = S.numResidues;
  numAtoms = S.numAtoms;
  for (int i = 0; i < S.chainSize(); i++) {
    chains.push_back(new Chain(S[i]));
    chains.back()->setParent(this);
    chainsByID[S[i].getID()] = chains.back();
    chainsBySegID[S[i].getSegID()] = chains.back();
  }
  copyExtraMaps(S);
}

/* The assumption is that if a Structure is deleted, all
 * of its children objects are no longer needed and should
 * go away. If a user needs to hold on to these, they
 * should generate copies as needed via copy constructors. */
Structure::~Structure() {
  deletePointers();
  destroyExtraMaps();
}

void Structure::deletePointers() {
  for (int i = 0; i < chains.size(); i++) delete(chains[i]);
}

void Structure::reset() {
  deletePointers();
  chains.resize(0);
  chainsByID.clear();
  chainsBySegID.clear();
  sourceFile = "";
  numResidues = numAtoms = 0;
}

void Structure::readPDB(string pdbFile, string options) {
  sourceFile = pdbFile;
  int lastresnum = -999999;
  string lastresname = "";
  string lasticode = "";
  string lastchainID = "";
  Chain* chain = NULL;
  Residue* residue = NULL;

  // various parsing options (the wonders of dealing with the good-old PDB format)
  bool ter = true;                   // flag to indicate that chain terminus was reached. Initialize to true so as to create a new chain upon reading the first atom.
  bool usesegid = false;             // use segment IDs to name chains instead of chain IDs? (useful when the latter are absent OR when too many chains, so need multi-letter names)
  bool skipHetero = false;           // skip hetero-atoms?
  bool charmmFormat = false;         // the PDB file was written by CHARMM (slightly different format)
  bool charmm19Format = false;       // upon reading, convert from all-hydrogen topology (param22 and higher) to the CHARMM19 united-atom topology (matters for HIS protonation states)
  bool fixIleCD1 = true;             // rename CD1 in ILE to CD (as is standard in MM packages)
  bool iCodesAsSepResidues = false;  // consequtive residues that differ only in their insertion code will be treated as separate residues
  bool uniqChainIDs = true;          // make sure chain IDs are unique, even if they are not unique in the read file

  // user-specified custom parsing options
  options = MstUtils::uc(options);
  if (options.find("USESEGID") != string::npos) usesegid = true;
  if (options.find("SKIPHETERO") != string::npos) skipHetero = true;
  if (options.find("CHARMM") != string::npos) charmmFormat = true;
  if (options.find("CHARMM19") != string::npos) charmm19Format = true;
  if (options.find("ALLOW DUPLICATE CIDS") != string::npos) uniqChainIDs = false;
  if (options.find("ALLOW ILE CD1") != string::npos) fixIleCD1 = false;
  if (options.find("ICODES AS RESIDUES") != string::npos) iCodesAsSepResidues = true;

  // read line by line
  string line;
  fstream ifh; MstUtils::openFile(ifh, pdbFile, fstream::in, "Structure::readPDB");

  while (getline(ifh, line)) {
    if (line.find("END") == 0) break;
    if (line.find("TER") == 0) { ter = true; continue; }
    if ((skipHetero && (line.find("ATOM") != 0)) || (!skipHetero && (line.find("ATOM") != 0) && (line.find("HETATM") != 0))) continue;

    /* Now read atom record. Sometimes PDB lines are too short (if they do not contain some
     * of the last optional columns). We don't want to read past the end of the string! */
    line += string(" ", 100);
    int atominx = MstUtils::toInt(line.substr(6, 5));
    string atomname = MstUtils::trim(line.substr(12, 4));
    string alt = line.substr(16, 1);
    string resname = MstUtils::trim(line.substr(17, 4));
    string chainID = MstUtils::trim(line.substr(21, 1));
    int resnum = charmmFormat ? MstUtils::toInt(line.substr(23, 4)) : MstUtils::toInt(line.substr(22, 4));
    string icode = charmmFormat ? " " : line.substr(26, 1);
    real x = MstUtils::toReal(line.substr(30, 8));
    real y = MstUtils::toReal(line.substr(38, 8));
    real z = MstUtils::toReal(line.substr(46, 8));
    string segID = MstUtils::trim(line.substr(72, 4));
    real B = MstUtils::toReal(line.substr(60, 6));
    real occ = MstUtils::toReal(line.substr(54, 6));
    bool het = (line.find("HETATM") == 0);

    // use segment ID's instead of chain ID's?
    if (usesegid) {
      chainID = segID;
    } else if ((chainID.compare("") == 0) && (segID.size() > 0) && (isalnum(segID[0]))) {
      // use first character of segment name if no chain name is specified, a segment ID is specified, and the latter starts with an alphanumeric character
      chainID = segID.substr(0, 1);
    }

    // create a new chain object, if necessary
    if ((chainID.compare(lastchainID) != 0) || ter) {
      chain = new Chain(chainID, segID);
      this->appendChain(chain, uniqChainIDs);
      // non-unique chains will be automatically renamed (unless the user specified not to rename chains), BUT we need to
      // remember the name that was actually read, since this name is what will be used to determine when the next chain comes
      if (chainID.compare(chain->getID())) {
        MstUtils::warn("chain name '" + chainID + "' was repeated in '" + pdbFile + "', renaming the chain to '" + chain->getID() + "'", "Structure::readPDB");
      }

      // start to count residue numbers in this chain
      lastresnum = -999999;
      lastresname = "";
      ter = false;
    }

    if (charmm19Format) {
      if (resname.compare("HSE") == 0) resname = "HSD";   // neutral HIS, proton on ND1
      if (resname.compare("HSD") == 0) resname = "HIS";   // neutral HIS, proton on NE2
      if (resname.compare("HSC") == 0) resname = "HSP";   // doubley-protonated +1 HIS
    }
    // many PDB files in the Protein Data Bank call the delta carbon of isoleucine CD1, but
    // the convention in basically all MM packages is to call it CD, since there is only one
    if (fixIleCD1 && atomname.compare("CD1") == 0) atomname = "CD"; 

    // if necessary, make a new residue
    bool reallyNewAtom = true; // is this a truely new atom, as opposed to an alternative position?
    if ((resnum != lastresnum) || (resname.compare(lastresname)) || (iCodesAsSepResidues && (icode.compare(lasticode))))  {
      residue = new Residue(resname, resnum, icode[0]);
      chain->appendResidue(residue);
    } else if (alt.compare(" ")) {
      // if this is not a new residue AND the alternative location flag is specified,
      // figure out if another location for this atom has already been given. If not,
      // then treat this as the "primary" location, and whatever other locations
      // are specified will be treated as alternatives.
      Atom* a = residue->findAtom(atomname, false);
      if (a) {
        reallyNewAtom = false;
        a->addAlternative(x, y, z, B, occ, alt[0]);
      }
    }
    // if necessary, make a new atom
    if (reallyNewAtom) {
      residue->appendAtom(new Atom(atominx, atomname, x, y, z, B, occ, het, alt[0]));
    }

    // remember previous values for determining whether something interesting happens next
    lastresnum = resnum;
    lasticode = icode;
    lastresname = resname;
    lastchainID = chainID;
  }
  ifh.close();
}

void Structure::writePDB(string pdbFile, string options) {
  fstream ofs; MstUtils::openFile(ofs, pdbFile, fstream::out, "Structure::writePDB");
  options = MstUtils::uc(options);

///  my $chainstr = shift; // probably want to implement this eventually. Or maybe some more generic selection mechanism based on regular expressions applied onto full atom strings.

  // various formating options (the wonders of dealing with the good-old PDB format)
  bool charmmFormat = false;         // the PDB file is intended for use in CHARMM or some other MM package
  bool charmm19Format = false;       // upon writing, convert from all-hydrogen topology (param 22 and higher) to CHARMM19 united-atom topology (matters for HIS protonation states)
  bool charmm22Format = false;       // upon writing, convert from CHARMM19 united-atom topology to all-hydrogen param 22 topology (matters for HIS protonation states). Also works for converting generic PDB files downloaded from the PDB.
  bool genericFormat = false;        // upon writing, convert to a generic PDB naming convention (no protonation state specified for HIS)
  bool fixIleCD1 = true;             // rename CD1 in ILE to CD in the output file (as is standard in MM packages)
  bool renumber = false;             // upon writing, renumber residue and atom names to start from 1 and go in order
  bool noend = false;                // do not write END at the end of the PDB file (e.g., useful for concatenating chains from several structures)
  bool noter = false;                // do not demark the end of each chain with TER (this is not _really_ necessary, assuming chain names are unique, and it is sometimes nice not to have extra lines other than atoms)

  // user-specified custom parsing options
  options = MstUtils::uc(options);
  if (options.find("CHARMM") != string::npos) charmmFormat = true;
  if (options.find("CHARMM19") != string::npos) charmm19Format = true;
  if (options.find("CHARMM22") != string::npos) charmm22Format = true;
  if (options.find("ALLOW ILE CD1") != string::npos) fixIleCD1 = false;
  if (options.find("ALLOW ILE CD1") != string::npos) fixIleCD1 = false;
  if (options.find("RENUMBER") != string::npos) renumber = false;
  if (options.find("NOEND") != string::npos) noend = true;
  if (options.find("NOTER") != string::npos) noter = true;
  if (charmm19Format && charmm22Format) MstUtils::error("CHARMM 19 and 22 formatting options cannot be specified together", "Structure::writePDB");

  int atomIndex = 0;
  for (int ci = 0; ci < this->chainSize(); ci++) {
    Chain& chain = (*this)[ci];
    for (int ri = 0; ri < chain.residueSize(); ri++) {
      Residue residue = chain[ri]; // NOTE: using a copy constructor here, in case residue details will be changed for formatting reasons upon writing
      for (int ai = 0; ai < residue.atomSize(); ai++) {
        Atom& atom = residue[ai]; // NOTE: no need to copy atom, since residue copying does a deep copy
        atomIndex++;
        // dirty details of formating for MM purposes converting
        if (charmmFormat) {
          if (residue.isNamed("ILE") && atom.isNamed("CD1")) atom.setName("CD");
          if (atom.isNamed("O") && (ri == chain.residueSize() - 1)) atom.setName("OT1");
          if (atom.isNamed("OXT") && (ri == chain.residueSize() - 1)) atom.setName("OT2");
          if (residue.isNamed("HOH")) residue.setName("TIP3");
        }
        if (charmm19Format) {
          if (residue.isNamed("HSD")) residue.setName("HIS"); // neutral HIS, proton on ND1
          if (residue.isNamed("HSE")) residue.setName("HSD"); // neutral HIS, proton on NE2
          if (residue.isNamed("HSC")) residue.setName("HSP"); // doubley-protonated +1 HIS
        } else if (charmm22Format) {
          /* this will convert from CHARMM19 to CHARMM22 as well as from a generic downlodaded
           * PDB file to one ready for use in CHARMM22. The latter is because in the all-hydrogen
           * topology, HIS protonation state must be explicitly specified, so there is no HIS per se.
           * Whereas in typical downloaded PDB files HIS is used for all histidines (usually, one
           * does not even really know the protonation state). Whether sometimes people do specify it
           * nevertheless, and what naming format they use to do so, I am not sure (welcome to the
           * PDB file format). But certainly almost always it is just HIS. Below HIS is renamed to
           * HSD, the neutral form with proton on ND1. This is an assumption; not a perfect one, but
           * something needs to be assumed. Doing this renaming will make the PDB file work in MM
           * packages with the all-hydrogen model. */
          if (residue.isNamed("HSD")) residue.setName("HSE"); // neutral HIS, proton on NE2
          if (residue.isNamed("HIS")) residue.setName("HSD"); // neutral HIS, proton on ND1
          if (residue.isNamed("HSP")) residue.setName("HSC"); // doubley-protonated +1 HIS
        } else if (genericFormat) {
          if (residue.isNamed("HSD")) residue.setName("HIS");
          if (residue.isNamed("HSP")) residue.setName("HIS");
          if (residue.isNamed("HSE")) residue.setName("HIS");
          if (residue.isNamed("HSC")) residue.setName("HIS");
          if (residue.isNamed("ILE") && atom.isNamed("CD")) atom.setName("CD1");
        }

        // write the atom line
        if (renumber) {
          ofs << atom.pdbLine(ri+1, atomIndex) << endl;
        } else {
          ofs << atom.pdbLine() << endl;
        }

      }
      if (!noter && (ri == chain.residueSize() - 1)) {
        ofs << "TER" << endl;
      }
    }
    if (!noend && (ci == this->chainSize() - 1)) {
      ofs << "END" << endl;
    }
  }
  ofs.close();
}

bool Structure::appendChain(Chain* C, bool allowRename) {
  chains.push_back(C);
  bool cidUnique = (chainsByID.find(C->getID()) == chainsByID.end());

  // if allowed to rename and there is a name clash, try to pick a unique chain name
  if (allowRename && !cidUnique) {
    string goodNames = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz";
    int f = -1;
    for (int i = 0; i < goodNames.size(); i++) {
      if (chainsByID.find(goodNames.substr(i, 1)) == chainsByID.end()) { f = i; break; }
    }
    if (f < 0) {
      MstUtils::warn("ran out of reasonable single-letter chain names, will use more than one character (PDB sctructure may be repeating chain IDs upon writing, but should still have unique segment IDs)!", "Structure::appendChain");
      string longName; bool found = false;
      for (int i = 0; i < goodNames.size(); i++) {
        for (int k = 0; k < 999; k++) {
          longName = goodNames.substr(i, 1) + MstUtils::toString(k);
          if (chainsByID.find(longName) == chainsByID.end()) { found = true; break; }
        }
        if (found) break;
      }
      MstUtils::assert(found, "ran out of even multi-character chain names -- your PDB structure really has more than 36,000 chains???", "Structure::appendChain");
      C->setID((string) longName);
      C->setSegID(C->getID());
    } else {
      C->setID(goodNames.substr(f, 1));
      cidUnique = true;
    }
  }

  chainsByID[C->getID()] = C;
  chainsBySegID[C->getSegID()] = C;
  C->setParent(this);
  numAtoms += C->atomSize();
  numResidues += C->residueSize();
  updateExtraMapsChainAdd(C);
  return cidUnique;
}

vector<Atom*> Structure::getAtoms() {
  vector<Atom*> vec;
  
  for (int i = 0; i < this->chainSize(); i++) {
    Chain& c = (*this)[i];
    for (int j = 0; j < c.residueSize(); j++) {
      Residue& r = c[j];
      vector<Atom*> resAtoms = r.getAtoms();
      vec.insert(vec.end(), resAtoms.begin(), resAtoms.end());
    }
  }

  return vec;
}

vector<Residue*> Structure::getResidues() {
  vector<Residue*> vec;
  
  for (int i = 0; i < this->chainSize(); i++) {
    Chain& c = (*this)[i];
    vector<Residue*> chainResidues = c.getResidues();
    vec.insert(vec.end(), chainResidues.begin(), chainResidues.end());
  }

  return vec;
}

void Structure::addAtom(Atom* A) {
  if ((A->getParent() == NULL) || (A->getParent()->getParent() == NULL)) MstUtils::error("cannot add a disembodies Atom", "Structure::addAtom");
  Residue* oldResidue = A->getParent();
  Chain* oldChain = oldResidue->getParent();
  Chain* newChain; Residue* newResidue; Atom* newAtom;

  // is there a chain matching the Atom's parent chain? If not, create one.
  newChain = getChainByID(oldChain->getID());
  if (newChain == NULL) {
    newChain = new Chain(oldChain->getID(), oldChain->getSegID());
    this->appendChain(newChain);
  }

  // is there a residue matching the Atoms parent residue? If not, create one.
  newResidue = newChain->findResidue(oldResidue->getName(), oldResidue->getNum());
  if (newResidue == NULL) {
    newResidue = new Residue(oldResidue->getName(), oldResidue->getNum(), oldResidue->getIcode());
    newChain->appendResidue(newResidue);
  }

  // finally, insert the atom into the correct residue
  newAtom = new Atom(*A);
  newResidue->appendAtom(newAtom);
}

void Structure::addAtoms(vector<Atom*>* atoms) {
  for (int i = 0; i < atoms->size(); i++) addAtom((*atoms)[i]);
}

void Structure::initExtraMaps() {
  residueIndexInChain = NULL;
  residueAtomsByName = NULL;
}

void Structure::destroyExtraMaps() {
  unmapResidueOrder();
  unmapAtomNames();
}

void Structure::copyExtraMaps(Structure& S) {
  if (S.residueOrderMap() != NULL) mapResidueOrder();
  if (S.atomNameMap() != NULL) mapAtomNames();
}

void Structure::mapResidueOrder() {
  unmapResidueOrder();
  residueIndexInChain = new map<Residue*, int>();
  for (int ci = 0; ci < chainSize(); i++) {
    Chain& chain = this->getChain(ci);
    for (int ri = 0; ri < chain.residueSize(); ri++) {
      (*residueIndexInChain)[&(chain[ri])] = ri;
    }
  }
}

void Structure::unmapResidueOrder() {
  if (residueIndexInChain != NULL) delete residueIndexInChain;
  residueIndexInChain = NULL;
}

void Structure::mapAtomNames() {
  unmapAtomNames();
  residueAtomsByName = new map<Residue*, map<string, Atom*> >();
  for (int ci = 0; ci < chainSize(); i++) {
    Chain& chain = this->getChain(ci);
    for (int ri = 0; ri < chain.residueSize(); ri++) {
      Residue& res = chain.getResidue(ri);
      for (int ai = 0; ai < res.atomSize(); ai++) {
        Atom& a = res.getAtom(ai);
        residueAtomsByName[&res][a.getName()] = &a;
      }
    }
  }
}

void Structure::unmapAtomNames() {
  if (residueAtomsByName != NULL) delete residueAtomsByName;
  residueAtomsByName = NULL;
}

void Structure::updateExtraMapsAtomAdd(Atom* atom) {
  if (residueAtomsByName != NULL) {
    Residue* res = newAtom->getResidue();
    if (res != NULL) residueAtomsByName[res][atom->getName()] = atom;
  }
}

void Structure::updateExtraMapsAtomRemove(Atom* atom) {
  if (residueAtomsByName != NULL) {
    Residue* res = newAtom->getResidue();
    if (res != NULL) residueAtomsByName[res].erase(atom->getName());
  }
}

void Structure::updateExtraMapsResidueAdd(Residue* res) {
  if (residueIndexInChain != NULL) residueIndexInChain.erase(res);
  if (residueAtomsByName != NULL) residueAtomsByName.erase(res);
}

void Structure::updateExtraMapsResidueRemove(Residue* res) {
  if (residueIndexInChain != NULL) residueIndexInChain.erase(res);
  if (residueAtomsByName != NULL) residueAtomsByName.erase(res);
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
    residues.back()->setParent(this);
  }
}

Chain::Chain(string chainID, string segID) {
  numAtoms = 0;
  parent = NULL;
  cid = chainID;
  sid = segID;
}

Chain::~Chain() {
  for (int i = 0; i < residues.size(); i++) delete(residues[i]);
}

vector<Atom*> Chain::getAtoms() {
  vector<Atom*> vec;

  for (int i = 0; i < residues.size(); i++) {
    Residue& r = getResidue(i);
    vector<Atom*> resAtoms = r.getAtoms();
    vec.insert(vec.end(), resAtoms.begin(), resAtoms.end());
  }

  return vec;
}

void Chain::incrementNumAtoms(int delta) {
  numAtoms += delta;
  if (parent != NULL) {
    parent->incrementNumAtoms(delta);
  }
}

void Chain::appendResidue(Residue* R) {
  incrementNumAtoms(R->atomSize());
  if (parent != NULL) {
    parent->incrementNumResidues();
  }
  residues.push_back(R);
  R->setParent(this);
  Structure* S = getStructure();
  if (S != NULL) updateExtraMapsResidueAdd();
}

Residue* Chain::findResidue(string resname, int resnum) {
  for (int i = 0; i < residueSize(); i++) {
    Residue& res = (*this)[i];
    if ((res.getNum() == resnum) && (res.isNamed(resname))) return &res;
  }
  return NULL;
}

Residue* Chain::findResidue(string resname, int resnum, char icode) {
  for (int i = 0; i < residueSize(); i++) {
    Residue& res = (*this)[i];
    if ((res.getNum() == resnum) && (res.isNamed(resname)) && (res.getIcode() == icode)) return &res;
  }
  return NULL;
}

/* --------- Residue --------- */
Residue::Residue() {
  resname = "UNK";
  resnum = 1;
  parent = NULL;
  icode = ' ';
}

Residue::Residue(Residue& R) {
  parent = R.parent;
  for (int i = 0; i < R.atomSize(); i++) {
    atoms.push_back(new Atom(R[i]));
    atoms.back()->setParent(this);
  }
  resnum = R.resnum;
  resname = R.resname;
  icode = R.icode;
}

Residue::Residue(string _resname, int _resnum, char _icode) {
  resname = _resname;
  resnum = _resnum;
  parent = NULL;
  icode = _icode;
}

Residue::~Residue() {
  for (int i = 0; i < atoms.size(); i++) delete(atoms[i]);
}

void Residue::appendAtom(Atom* A) {
  atoms.push_back(A);
  A->setParent(this);
  if (parent != NULL) {
    parent->incrementNumAtoms();
  }
}

void Residue::deleteAtom(int i) {
  if ((i < 0) || (i >= atoms.size())) {
    MstUtils::error("index out of range of atom vector in residue", "Residue::deleteAtom");
  }
  Structure* S = getStructure(); if (S != NULL) S->updateExtraMapsAtomRemove(atoms[i]);
  atoms.erase(atoms.begin() + i);
  if (parent != NULL) {
    parent->incrementNumAtoms(-1);
  }
}

void Residue::deleteAtoms() {
  if (parent != NULL) {
    parent->incrementNumAtoms(-atoms.size());
  }
  Structure* S = getStructure();
  for (int i = 0; i < atoms.size(); i++) {
    if (S != NULL) S->updateExtraMapsAtomRemove(atoms[i]);
    delete atoms[i];
  }
  atoms.resize(0);
}

void Residue::replaceAtoms(vector<Atom*>& newAtoms, vector<int>* toRemove) {
  bool delAll = (toRemove == NULL);
  int N = delAll ? atoms.size() : toRemove->size();
  if (parent != NULL) {
    parent->incrementNumAtoms(newAtoms.size() - N);
  }
  Structure* S = getStructure();

  // delete those atoms needing deletion
  for (int i = 0; i < N; i++) {
    int ai = delAll ? i : toRemove->[i];
    if ((ai < 0) || (ai >= atoms.size())) {
      MstUtils::error("index out of range of atom vector in residue", "Residue::replaceAtoms");
    }
    if (S != NULL) S->updateExtraMapsAtomRemove(atoms[ai]);
    delete atoms[ai];
    atoms[ai] = NULL
  }

  // create a new atom vector without them
  vector<Atom*> oldAtoms = atoms;
  atoms.resize(atoms.size() - toRemove->size());
  int k = 0;
  for (int i = 0; i < oldAtoms.size(); i++) {
    if (oldAtoms[i] == NULL) continue;
    atoms[k] = oldAtoms[i];
    k++;
  }

  // and now append the new atoms
  for (int i = 0; i < newAtoms.size(); i++) {
    atoms[k] = newAtoms[i];
    atoms[k]->setParent(this);
    if (S != NULL) S->updateExtraMapsAtomAdd(atoms[k]);
    k++;
  }
}

Residue* Residue::iPlusDelta(int off) {
  Structure* S = this->getStructure();
  if ((S != NULL) && (S->residueOrderMap() != NULL)) {
    map<Residue*, int>* residueIndexInChain = S->residueOrderMap();
    if (residueIndexInChain.find(this) == residueIndexInChain.end()) {
      MstUtils::error("Residue " + MstUtils::toString(*this) + " pointed to by " + this + " does not appear in the residue index map of its parent!", "Residue::nextResidue");
    }
    int i = residueIndexInChain[this];
    Chain& chain = getChain();
    if ((i + off >= chain.residueSize()) || (i + off < 0)) return NULL;
    return chain.getResidue(i + off);    
  } else {
    Chain& chain = getChain();
    if (chain == NULL) {
      MstUtils::error("Residue " + MstUtils::toString(*this) + " is disembodied!", "Residue::iPlusDelta");
    }
    // find the index of the current residue in the chain by scanning
    for (int i = 0; i < chain.residueSize(); i++) {
      if (this == chain.getResidue(i)) {
        if ((i + off >= chain.residueSize()) || (i + off < 0)) return NULL;
        return chain.getResidue(i + off);    
      }
    }
    MstUtils::error("Residue " + MstUtils::toString(*this) + " does not appear to be in the chain it points back to, something went wrong!", "Residue::iPlusDelta");
  }
}

Residue* Residue::nextResidue() {
  return iPlusDelta(1);
}

Residue* Residue::previousResidue() {
  return iPlusDelta(-1);
}

// C- N  CA  C 
real Residue::getPhi() {
  Residue* res1 = previousResidue();
  if (res1 == NULL) return badDihedral;
  Atom* A = res1->findAtom("C", false);
  Atom* B = findAtom("N", false);
  Atom* C = findAtom("CA", false);
  Atom* D = findAtom("C", false);
  if ((A == NULL) || (B == NULL) || (C == NULL) || (D == NULL)) {
    MstUtils::error("not all backbone atoms present to compute PHI for residue " + MstUtils::toString(*this));
  }

  return CartesianGeometry::dihedral(*A, *B, *C, *D);
}

// N  CA  C  N+
real Residue::getPsi(Residue* res, bool strict) {
  Residue* res1 = nextResidue();
  if (res1 == NULL) return badDihedral;
  Atom* A = findAtom("N", false);
  Atom* B = findAtom("CA", false);
  Atom* C = findAtom("C", false);
  Atom* D = res1->findAtom("N", false);
  if ((A == NULL) || (B == NULL) || (C == NULL) || (D == NULL)) {
    MstUtils::error("not all backbone atoms present to compute PHI for residue " + MstUtils::toString(*this));
  }
}

Atom* Residue::findAtom(string _name, bool strict) {
  for (int i = 0; i < atoms.size(); i++) {
    if (atoms[i]->isNamed(_name)) return atoms[i];
  }
  if (strict) MstUtils::error("could not find atom named '" + _name + "' in residue " + MstUtils::toString(*this), "Residue::findAtom");
  return NULL;
}

string Residue::getChainID(bool strict) {
  if (parent == NULL) {
    if (strict) MstUtils::error("residue has no parent", "Residue::getChainID");
    return "";
  }
  return parent->getID();
}

/* --------- Atom --------- */
Atom::Atom() {
  parent = NULL;
  het = false;
  name = MstUtils::copyStringC("UNK");
  setName("");
  alternatives = NULL;
}

Atom::Atom(Atom& A) {
  index = A.index;
  name = NULL;
  setName(A.name);
  x = A.x;
  y = A.y;
  z = A.z;
  B = A.B;
  occ = A.occ;
  het = A.het;
  alt = A.alt;
  parent = A.parent;
  if (A.alternatives != NULL) {
    alternatives = new vector<altInfo>(A.alternatives->size());
    for (int i = 0; i < A.alternatives->size(); i++) {
      (*alternatives)[i] = (*(A.alternatives))[i];
    }
  } else {
    alternatives = NULL;
  }
}

Atom::Atom(int _index, string _name, real _x, real _y, real _z, real _B, real _occ, bool _het, char _alt, Residue* _parent) {
  index = _index;
  name = NULL;
  setName(_name);
  x = _x;
  y = _y;
  z = _z;
  B = _B;
  occ = _occ;
  het = _het;
  alt = _alt;
  parent = _parent;
  alternatives = NULL;
}

Atom::~Atom() {
  if (name != NULL) delete[] name;
  if (alternatives != NULL) delete alternatives;
}

void Atom::setName(const char* _name) {
  if (name != NULL) delete[] name;
  name = new char[strlen(_name)+1];
  strcpy(name, _name);
}

void Atom::setName(string _name) {
  if (name != NULL) delete[] name;
  name = new char[_name.size()+1];
  strcpy(name, _name.c_str());
}

void Atom::swapWithAlternative(int altInd) {
  if ((alternatives == NULL) || (altInd >= alternatives->size())) MstUtils::error("specified alternative index out of bounds", "Atom::swapWithAlternative");
  altInfo& targ = (*alternatives)[altInd];
  altInfo temp = targ;
  targ.x = x; targ.y = y; targ.z = z; targ.occ = occ; targ.B = B; targ.alt = alt;
  x = temp.x; y = temp.y; z = temp.z; occ = temp.occ; B = temp.B; alt = temp.alt;
}

void Atom::addAlternative(real _x, real _y, real _z, real _B, real _occ, char _alt) {
  if (alternatives == NULL) {
    alternatives = new vector<altInfo>(0);
  }
  alternatives->push_back(altInfo(_x, _y, _z, _B, _occ, _alt));
}

string Atom::pdbLine(int resIndex, int atomIndex) {
  char line[100]; // a PDB line is at most 80 characters, so this is plenty
  string resname = "UNK"; string chainID = "?"; string segID = "?";
  int resnum = 1; char icode = ' ';

  // chain and residue info
  if (parent != NULL) {
    resname = parent->getName();
    if (resname.length() > 4) resname = resname.substr(0, 4);
    resnum = parent->getNum();
    icode = parent->getIcode();
    Chain* chain = parent->getParent();
    if (chain != NULL) {
      chainID = chain->getID();
      if (chainID.length() > 1) chainID = chainID.substr(0, 1);
      segID = chain->getSegID();
      if (segID.length() > 4) segID = segID.substr(0, 4);
    }
  }

  // atom name placement is different when it is 4 characters long
  char atomname[5];
  if (strlen(name) < 4) { sprintf(atomname, " %-.3s", name); }
  else { sprintf(atomname, "%.4s", name); }

  // moduli are used to make sure numbers do not go over prescribe field widths (this is not enforced by sprintf like with strings)
  sprintf(line, "%6s%5d %-4s%c%-4s%.1s%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %.4s", 
          isHetero() ? "HETATM" : "ATOM  ", atomIndex % 100000, atomname, alt, resname.c_str(), chainID.c_str(),
          resnum % 10000, icode, x, y, z, occ, B, segID.c_str());

  return (string) line;
}

real Atom::distance(Atom& another) {
  return sqrt((x - another.x)*(x - another.x) + (y - another.y)*(y - another.y) + (z - another.z)*(z - another.z));
}

real Atom::distance2(Atom& another) {
  return (x - another.x)*(x - another.x) + (y - another.y)*(y - another.y) + (z - another.z)*(z - another.z);
}

/* --------- AtomPointerVector --------- */

CartesianPoint AtomPointerVector::getGeometricCenter() {
  CartesianPoint C(3, 0);
  for (int i = 0; i < this->size(); i++) {
    C += CartesianPoint(*((*this)[i]));
  }
  C /= this->size();
  return C;
}

/* --------- CartesianPoint --------- */

CartesianPoint::CartesianPoint(const Atom& A) {
  this->resize(3, 0);
  (*this)[0] = A.getX();
  (*this)[1] = A.getY();
  (*this)[2] = A.getZ();
}

CartesianPoint& CartesianPoint::operator+=(const CartesianPoint &rhs) {
  if (this->size() != rhs.size()) MstUtils::error("points of different dimensionality!", "CartesianPoint::operator+=");
  for (int i = 0; i < this->size(); i++) {
    (*this)[i] += rhs[i];
  }
  return *this;
}

CartesianPoint& CartesianPoint::operator-=(const CartesianPoint &rhs) {
  if (this->size() != rhs.size()) MstUtils::error("points of different dimensionality!", "CartesianPoint::operator-=");
  for (int i = 0; i < this->size(); i++) {
    (*this)[i] -= rhs[i];
  }
  return *this;
}

CartesianPoint& CartesianPoint::operator*=(const real& s) {
  for (int i = 0; i < this->size(); i++) {
    (*this)[i] *= s;
  }
  return *this;
}

CartesianPoint& CartesianPoint::operator/=(const real& s) {
  for (int i = 0; i < this->size(); i++) {
    (*this)[i] /= s;
  }
  return *this;
}

const CartesianPoint CartesianPoint::operator+(const CartesianPoint &other) const {
  CartesianPoint result = *this;
  result += other;
  return result;
}

const CartesianPoint CartesianPoint::operator-(const CartesianPoint &other) const {
  CartesianPoint result = *this;
  result -= other;
  return result;
}

const CartesianPoint CartesianPoint::operator*(const real& s) const {
  CartesianPoint result = *this;
  result *= s;
  return result;
}

const CartesianPoint CartesianPoint::operator/(const real& s) const {
  CartesianPoint result = *this;
  result /= s;
  return result;
}

CartesianPoint& CartesianPoint::operator=(const Atom& A) {
  *this = CartesianPoint(A);
}

real CartesianPoint::norm() const {
  return sqrt((*this)[0]*(*this)[0] + (*this)[1]*(*this)[1] + (*this)[2]*(*this)[2]);
}

CartesianPoint CartesianPoint::cross(CartesianPoint other) const {
  if (size() != 3) MstUtils::error("don't know how to compute cross produces for dimensions other than 3", "CartesianPoint::cross");
  if (size() != other.size()) MstUtils::error("vector size mismatch", "CartesianPoint::cross");

  CartesianPoint C(3, 0);
  C[0] = getY()*other.getZ() - getZ()*other.getY();
  C[1] = getZ()*other.getX() - getX()*other.getZ();
  C[2] = getX()*other.getY() - getY()*other.getX();

  return C;
}

real CartesianPoint::dot(CartesianPoint other) const {
  if (size() != other.size()) MstUtils::error("vector size mismatch", "CartesianPoint::dot");

  real d = 0;
  for (int i = 0; i < size(); i++) {
    d += (*this)[i] * other[i];
  }

  return d;
}

/* --------- CartesianGeometry --------- */
real CartesianGeometry::dihedralRadians(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4) {
  CartesianPoint AB = _p1 - _p2;
  CartesianPoint CB = _p3 - _p2;
  CartesianPoint DC = _p4 - _p3;

  if (AB.norm() == 0.0 || CB.norm() == 0.0 || DC.norm() == 0.0) MstUtils::error("some points coincide in dihedral calculation", "CartesianGeometry::dihedralRadians");

  CartesianPoint ABxCB = AB.cross(CB).getUnit();
  CartesianPoint DCxCB = DC.cross(CB).getUnit();

  // the following is necessary for values very close to 1 but just above
  double dotp = ABxCB * DCxCB;
  if (dotp > 1.0) {
    dotp = 1.0;
  } else if (dotp < -1.0) {
    dotp = -1.0;
  }

  double angle = acos(dotp);
  if (ABxCB * DC > 0) {
    angle *= -1;
  }
  return angle;
}

real CartesianGeometry::dihedralRadians(const CartesianPoint * _p1, const CartesianPoint * _p2, const CartesianPoint * _p3, const CartesianPoint * _p4) {
  return dihedralRadians(*_p1, *_p2, *_p3, *_p4);
}

real CartesianGeometry::dihedral(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4) {
  return dihedralRadians(_p1, _p2, _p3, _p4)*180/M_PI;
}

real CartesianGeometry::dihedral(const CartesianPoint * _p1, const CartesianPoint * _p2, const CartesianPoint * _p3, const CartesianPoint * _p4) {
  return dihedralRadians(*_p1, *_p2, *_p3, *_p4)*180/M_PI;
}

/* --------- MstUtils --------- */
void MstUtils::openFile (fstream& fs, string filename, ios_base::openmode mode, string from) {
  fs.open(filename.c_str(), mode);
  if (!fs.is_open()) {
    if (!from.empty()) from =+ " -> ";
    MstUtils::error("could not open file " + filename, from + "MstUtils::openFile");
  }
}

FILE* MstUtils::openFileC (const char* filename, const char* mode, string from) {
  FILE * fp;
  fp = fopen (filename, mode);
  if (fp == NULL) {
    if (!from.empty()) from =+ " -> ";
    MstUtils::error("could not open file " + (string) filename, from + "MstUtils::openFileC");
  }
  return fp;
}

string MstUtils::uc(string& str){
  string ret = str;
  for (int i = 0; i < ret.length(); i++) {
    ret[i] = toupper(ret[i]);
  }
  return ret;
}

string MstUtils::trim(string str) {
  if (str.empty()) return str;

  int b = 0; int e = str.length();
  for (int i = 0; i < str.length(); i++) {
    if (!isspace(str[i])) break;
    b = i + 1;
  }
  if (b == str.length()) return "";
  for (int i = str.length()-1; i >= 0; i--) {
    if (!isspace(str[i])) break;
    e = i;
  }
  return str.substr(b, e - b);
}

void MstUtils::warn(string message, string from) {
  string head = from.empty() ? "Warning: " : "Warning in " + from + ": ";
  cerr << head << wrapText(message, 100, 0, head.length()) << endl;
}

void MstUtils::error(string message, string from, int code) {
  string head = from.empty() ? "Error: " : "Error in " + from + ": ";
  cerr << head << wrapText(message, 100, 0, head.length()) << endl;
  exit(code);
}

void MstUtils::assert(bool condition, string message, string from, int exitCode) {
  if(!condition) {
    MstUtils::error(message, from, exitCode);
  }
}

string MstUtils::wrapText(string message, int width, int leftSkip, int startingOffset) {
  string text;
  int b = 0, e, off = startingOffset, word = 0;
  while (b < message.size()) {
    // find the end of the next word
    e = message.find_first_of(" ", b);
    if (e == string::npos) e = message.size();
    int n = e - b;
    // if including the next word on this line will go over the width, start a new line
    if ((off + n >= width) && (word > 0)) {
      text += "\n" + string(leftSkip, ' ');
      off = leftSkip;
      word = 0;
    }
    text += message.substr(b, n) + " ";
    off += n + 1;
    b = e+1;
    word++;
  }
  return text;
}

char* MstUtils::copyStringC(const char* str) {
  char* copy = new char[strlen(str) + 1];
  strcpy(copy, str);
  return copy;
}

int MstUtils::toInt(string num, bool strict) {
  int ret;
  if ((sscanf(num.c_str(), "%d", &ret) == 0) && strict) MstUtils::error("failed to convert '" + num + "' to integer", "MstUtils::toInt");
  return ret;
}

MST::real MstUtils::toReal(string num, bool strict) {
  double ret;
  if ((sscanf(num.c_str(), "%lf", &ret) == 0) && strict) MstUtils::error("failed to convert '" + num + "' to real", "MstUtils::toReal");
  return (real) ret;
}

MST::real MstUtils::mod(MST::real num, MST::real den) {
  return num - ((MST::real) floor((double) num / (double) den)) * den;
}
