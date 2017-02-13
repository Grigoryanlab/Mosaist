#include "msttypes.h"

using namespace MST;

/* --------- Structure --------- */
Structure::Structure() {
  numResidues = numAtoms = 0;
}

Structure::Structure(string pdbFile, string options) {
  name = pdbFile;
  numResidues = numAtoms = 0;
  readPDB(pdbFile, options);
}

Structure::Structure(Structure& S) {
  name = S.name;
  numResidues = S.numResidues;
  numAtoms = S.numAtoms;
  for (int i = 0; i < S.chainSize(); i++) {
    chains.push_back(new Chain(S[i]));
    chains.back()->setParent(this);
    chainsByID[S[i].getID()] = chains.back();
    chainsBySegID[S[i].getSegID()] = chains.back();
  }
}

Structure::Structure(Chain& C) {
  numResidues = numAtoms = 0;
  appendChain(new Chain(C));
}

Structure::Structure(Residue& R) {
  numResidues = numAtoms = 0;
  Chain* newChain = appendChain("A", true);
  newChain->appendResidue(new Residue(R));
}

/* The assumption is that if a Structure is deleted, all
 * of its children objects are no longer needed and should
 * go away. If a user needs to hold on to these, they
 * should generate copies as needed via copy constructors. */
Structure::~Structure() {
  deletePointers();
}

void Structure::deletePointers() {
  for (int i = 0; i < chains.size(); i++) delete(chains[i]);
}

void Structure::reset() {
  deletePointers();
  chains.resize(0);
  chainsByID.clear();
  chainsBySegID.clear();
  name = "";
  numResidues = numAtoms = 0;
}

void Structure::readPDB(string pdbFile, string options) {
  name = pdbFile;
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
  bool ignoreTER = false;            // if true, will not pay attention to TER lines in deciding when chains end/begin

  // user-specified custom parsing options
  options = MstUtils::uc(options);
  if (options.find("USESEGID") != string::npos) usesegid = true;
  if (options.find("SKIPHETERO") != string::npos) skipHetero = true;
  if (options.find("CHARMM") != string::npos) charmmFormat = true;
  if (options.find("CHARMM19") != string::npos) charmm19Format = true;
  if (options.find("ALLOW DUPLICATE CIDS") != string::npos) uniqChainIDs = false;
  if (options.find("ALLOW ILE CD1") != string::npos) fixIleCD1 = false;
  if (options.find("ICODES AS RESIDUES") != string::npos) iCodesAsSepResidues = true;
  if (options.find("IGNORE-TER") != string::npos) ignoreTER = true;

  // read line by line
  string line;
  fstream ifh; MstUtils::openFile(ifh, pdbFile, fstream::in, "Structure::readPDB");

  while (getline(ifh, line)) {
    if (line.find("END") == 0) break;
    if ((line.find("TER") == 0) && !ignoreTER) { ter = true; continue; }
    if ((skipHetero && (line.find("ATOM") != 0)) || (!skipHetero && (line.find("ATOM") != 0) && (line.find("HETATM") != 0))) continue;

    /* Now read atom record. Sometimes PDB lines are too short (if they do not contain some
     * of the last optional columns). We don't want to read past the end of the string! */
    line += string(100, ' ');
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
  fstream ofs; MstUtils::openFile(ofs, pdbFile, fstream::out, "Structure::writePDB(string, string)");
  writePDB(ofs, options);
  ofs.close();
}

void Structure::writePDB(fstream& ofs, string options) {
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
  if (options.find("RENUMBER") != string::npos) renumber = true;
  if (options.find("NOEND") != string::npos) noend = true;
  if (options.find("NOTER") != string::npos) noter = true;
  if (charmm19Format && charmm22Format) MstUtils::error("CHARMM 19 and 22 formatting options cannot be specified together", "Structure::writePDB");

  int atomIndex = 0;
  for (int ci = 0; ci < this->chainSize(); ci++) {
    Chain& chain = (*this)[ci];
    for (int ri = 0; ri < chain.residueSize(); ri++) {
      Residue residue = chain[ri]; // NOTE: using a copy constructor here, in case residue details will be changed for formatting reasons upon writing
      residue.setParent(&chain);   // by default, objects are copied as being disembodied (so as not to create inconsistent states)
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
  return cidUnique;
}

Chain* Structure::appendChain(string cid, bool allowRename) {
  Chain* newChain = new Chain(cid, cid);
  this->appendChain(newChain, allowRename);
  return newChain;
}

Residue& Structure::getResidue(int i) {
  if ((i < 0) && (i >= residueSize()))
    MstUtils::error("residue index " + MstUtils::toString(i) + " out of range for Structure", "Structure::getResidue(int)");
  for (int ci = 0; ci < chainSize(); ci++) {
    Chain& chain = getChain(ci);
    if (i >= chain.residueSize()) {
      i -= chain.residueSize();
    } else {
      return chain[i];
    }
  }
  MstUtils::error("something strange happened; most likely, various counters are inconsistent in Structure object", "Structure::getResidue(int)");
  return *(new Residue()); // just to make the compiler happy and not throw a warning; this is never reached
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

void Structure::renumber() {
  int index = 1;
  for (int i = 0; i < chains.size(); i++) {
    Chain& chain = (*this)[i];
    for (int j = 0; j < chain.residueSize(); j++) {
      Residue& res = chain[j];
      res.setNum(j+1);
      for (int k = 0; k < res.atomSize(); k++) {
        Atom& a = res[k];
        a.setIndex(index);
        index++;
      }
    }
  }
}

void Structure::addAtom(Atom* A) {
  if ((A->getParent() == NULL) || (A->getParent()->getParent() == NULL)) MstUtils::error("cannot add a disembodied Atom", "Structure::addAtom");
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


/* --------- Chain --------- */
Chain::Chain() {
  numAtoms = 0;
  parent = NULL;
}

Chain::Chain(Chain& C) {
  numAtoms = C.numAtoms;
  parent = NULL;
  for (int i = 0; i < C.residueSize(); i++) {
    residues.push_back(new Residue(C[i]));
    residues.back()->setParent(this);
    residueIndexInChain[residues.back()] = i;
  }
  cid = C.cid;
  sid = C.sid;
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

int Chain::getResidueIndex(Residue* res) {
  if (residueIndexInChain.find(res) == residueIndexInChain.end())
    MstUtils::error("passed residue does not appear in chain's index map", "Chain::residueIndexInChain");
  return residueIndexInChain[res];
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
  residueIndexInChain[R] = residues.size() - 1;
}

void Chain::insertResidue(Residue* R, int index) {
  if ((index < 0) || (index >= residues.size()))
    MstUtils::error("residue index '" + MstUtils::toString(index) + "' out of range", "Chain::insertResidue(Residue*, int)");
  incrementNumAtoms(R->atomSize());
  if (parent != NULL) {
    parent->incrementNumResidues();
  }
  residues.insert(residues.begin() + index, R);
  R->setParent(this);
  residueIndexInChain[R] = index;
  for (int i = index + 1; i < residues.size(); i++) {
    residueIndexInChain[residues[i]]++;
  }
}
Residue* Chain::insertResidueCopy(Residue* R, int index) {
  Residue* newRes = new Residue(*R);
  if (index == -1) this->appendResidue(newRes);
  else this->insertResidue(newRes, index);
  return newRes;
}

Residue* Chain::insertResidueCopy(Residue& R, int index) {
  Residue* newRes = new Residue(R);
  if (index == -1) this->appendResidue(newRes);
  else this->insertResidue(newRes, index);
  return newRes;
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

const real Residue::badDihedral = 999.0;

Residue::Residue() {
  resname = "UNK";
  resnum = 1;
  parent = NULL;
  icode = ' ';
}

Residue::Residue(Residue& R, bool copyAlt) {
  parent = NULL;
  for (int i = 0; i < R.atomSize(); i++) {
    atoms.push_back(new Atom(R[i], copyAlt));
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

void Residue::appendAtoms(vector<Atom*>& A) {
  for (int i = 0; i < A.size(); i++) appendAtom(A[i]);
}

void Residue::deleteAtom(int i) {
  if ((i < 0) || (i >= atoms.size())) {
    MstUtils::error("index out of range of atom vector in residue", "Residue::deleteAtom");
  }
  atoms.erase(atoms.begin() + i);
  if (parent != NULL) {
    parent->incrementNumAtoms(-1);
  }
}

void Residue::copyAtoms(Residue& R, bool copyAlt) {
  deleteAtoms();
  for (int i = 0; i < R.atomSize(); i++) {
    atoms.push_back(new Atom(R[i], copyAlt));
    atoms.back()->setParent(this);
  }
  if (parent != NULL) parent->incrementNumAtoms(R.atomSize());
}

void Residue::makeAlternativeMain(int altInd) {
  for (int i = 0; i < atomSize(); i++) {
    (*this)[i].makeAlternativeMain(altInd);
  }
}

void Residue::deleteAtoms() {
  if (parent != NULL) {
    parent->incrementNumAtoms(-atoms.size());
  }
  Structure* S = getStructure();
  for (int i = 0; i < atoms.size(); i++) {
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
    int ai = delAll ? i : (*toRemove)[i];
    if ((ai < 0) || (ai >= atoms.size())) {
      MstUtils::error("index out of range of atom vector in residue", "Residue::replaceAtoms");
    }
    delete atoms[ai];
    atoms[ai] = NULL;
  }

  // create a new atom vector without them
  vector<Atom*> oldAtoms = atoms;
  atoms.resize(atoms.size() - N + newAtoms.size());
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
    k++;
  }
}

Residue* Residue::iPlusDelta(int off) {
  Chain* chain = this->getChain();
  if (chain == NULL) {
    MstUtils::error("Residue " + MstUtils::toString(*this) + " is disembodied!", "Residue::iPlusDelta");
  }
  int i = chain->getResidueIndex(this);
  if ((i + off >= chain->residueSize()) || (i + off < 0)) return NULL;
  return &(chain->getResidue(i + off));
}

Residue* Residue::nextResidue() {
  return iPlusDelta(1);
}

Residue* Residue::previousResidue() {
  return iPlusDelta(-1);
}

// C- N  CA  C
real Residue::getPhi(bool strict) {
  Residue* res1 = previousResidue();
  if (res1 == NULL) return badDihedral;
  Atom* A = res1->findAtom("C", false);
  Atom* B = findAtom("N", false);
  Atom* C = findAtom("CA", false);
  Atom* D = findAtom("C", false);
  if ((A == NULL) || (B == NULL) || (C == NULL) || (D == NULL)) {
    if (strict) MstUtils::error("not all backbone atoms present to compute PHI for residue " + MstUtils::toString(*this), "Residue::getPhi");
    return badDihedral;
  }

  return CartesianGeometry::dihedral(*A, *B, *C, *D);
}

// N  CA  C  N+
real Residue::getPsi(bool strict) {
  Residue* res1 = nextResidue();
  if (res1 == NULL) return badDihedral;
  Atom* A = findAtom("N", false);
  Atom* B = findAtom("CA", false);
  Atom* C = findAtom("C", false);
  Atom* D = res1->findAtom("N", false);
  if ((A == NULL) || (B == NULL) || (C == NULL) || (D == NULL)) {
    if (strict) MstUtils::error("not all backbone atoms present to compute PSI for residue " + MstUtils::toString(*this), "Residue::getPsi");
    return badDihedral;
  }

  return CartesianGeometry::dihedral(*A, *B, *C, *D);
}

// CA-  C-  N  CA
real Residue::getOmega(bool strict) {
  Residue* res1 = previousResidue();
  if (res1 == NULL) return badDihedral;
  Atom* A = res1->findAtom("CA", false);
  Atom* B = res1->findAtom("C", false);
  Atom* C = findAtom("N", false);
  Atom* D = findAtom("CA", false);
  if ((A == NULL) || (B == NULL) || (C == NULL) || (D == NULL)) {
    if (strict) MstUtils::error("not all backbone atoms present to compute OMEGA for residue " + MstUtils::toString(*this), "Residue::getOmega");
    return badDihedral;
  }

  return CartesianGeometry::dihedral(*A, *B, *C, *D);
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

Atom::Atom(Atom& A, bool copyAlt) {
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
  if (copyAlt && (A.alternatives != NULL)) {
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

real Atom::operator[](int i) const {
  switch(i) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      MstUtils::error("invalid coordinate index " + MstUtils::toString(i), "Atom::operator[](int)");
  }
  return 0.0; // just to silence the warning from some compilres; in reality, this is never reached
}

vector<real> Atom::getAltCoor(int altInd) {
  if ((alternatives == NULL) || (altInd >= alternatives->size()) || (altInd < 0)) MstUtils::error("alternative index " + MstUtils::toString(altInd) + " out of bounds (" + MstUtils::toString(alternatives->size()) + " alternatives available)", "Atom::swapWithAlternative");
  altInfo& targ = (*alternatives)[altInd];
  vector<real> coor;
  coor.push_back(targ.x); coor.push_back(targ.y); coor.push_back(targ.z);
  return coor;
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
  if ((alternatives == NULL) || (altInd >= alternatives->size())) MstUtils::error("alternative index " + MstUtils::toString(altInd) + " out of bounds (" + MstUtils::toString(alternatives->size()) + " alternatives available)", "Atom::swapWithAlternative");
  altInfo& targ = (*alternatives)[altInd];
  altInfo temp = targ;
  targ.x = x; targ.y = y; targ.z = z; targ.occ = occ; targ.B = B; targ.alt = alt;
  x = temp.x; y = temp.y; z = temp.z; occ = temp.occ; B = temp.B; alt = temp.alt;
}

void Atom::makeAlternativeMain(int altInd) {
  if ((alternatives == NULL) || (altInd >= alternatives->size())) MstUtils::error("alternative index " + MstUtils::toString(altInd) + " out of bounds (" + MstUtils::toString(alternatives->size()) + " alternatives available)", "Atom::swapWithAlternative");
  altInfo& targ = (*alternatives)[altInd];
  x = targ.x; y = targ.y; z = targ.z; occ = targ.occ; B = targ.B; alt = targ.alt;
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

real AtomPointerVector::radiusOfGyration() {
  CartesianPoint center = getGeometricCenter();
  real s = 0;
  for (int i = 0; i < size(); i++) {
      s += CartesianPoint((*this)[i]).distance2(center);
  }
  return sqrt(s / size());
}

void AtomPointerVector::deletePointers() {
  for (int i = 0; i < size(); i++) delete (*this)[i];
  resize(0);
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

const CartesianPoint CartesianPoint::operator-() const {
  return CartesianPoint(this->size(), 0) - *this;
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
  return *this;
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

real CartesianPoint::distance(CartesianPoint& another) {
  if (this->size() != another.size()) MstUtils::error("point dimensions disagree", "CartesianPoint::distance(CartesianPoint&)");
  real d = 0;
  for (int i = 0; i < this->size(); i++) {
    d += ((*this)[i] - another[i])*((*this)[i] - another[i]);
  }
  return sqrt(d);
}

real CartesianPoint::distance2(CartesianPoint& another) {
  if (this->size() != another.size()) MstUtils::error("point dimensions disagree", "CartesianPoint::distance2(CartesianPoint&)");
  real d = 0;
  for (int i = 0; i < this->size(); i++) {
    d += ((*this)[i] - another[i])*((*this)[i] - another[i]);
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


/* --------- RMSDCalculator --------- */

vector<real> RMSDCalculator::lastTranslation() {
    vector<real> trans(3, 0.0);
    for (int i = 0; i < 3; i++) trans[i] = t[i];
    return trans;
}

vector<vector<real> > RMSDCalculator::lastRotation() {
    vector<vector<real> > rot(3);
    for (int i = 0; i < 3; i++) {
        rot[i].resize(3, 0);
        for (int j = 0; j < 3; j++ ) {
            rot[i][j] = u[i][j];
        }
    }
    return rot;
}

real RMSDCalculator::bestRMSD(vector<Atom*> &_align, vector<Atom*> &_ref, bool* _suc, bool setTransRot) {
    _rmsd = 999999.0;
    if (Kabsch(_align, _ref, setTransRot)) { if (_suc != NULL) *_suc = true; }
    else { if (_suc != NULL) *_suc = false; }
    return _rmsd;
}

bool RMSDCalculator::align(vector<Atom*> &_align, vector<Atom*> &_ref, vector<Atom*>& _moveable) {
    _rmsd = 999999.0;
    bool suc = Kabsch(_align, _ref, 1);

    if (suc) {
        real x[3],x1[3];
        for(int k=0; k<_moveable.size(); k++) {
            x[0]=_moveable[k]->getX();
            x[1]=_moveable[k]->getY();
            x[2]=_moveable[k]->getZ();
            x1[0] = t[0]+u[0][0]*x[0]+u[0][1]*x[1]+u[0][2]*x[2];
            x1[1] = t[1]+u[1][0]*x[0]+u[1][1]*x[1]+u[1][2]*x[2];
            x1[2] = t[2]+u[2][0]*x[0]+u[2][1]*x[1]+u[2][2]*x[2];
            _moveable[k]->setCoor(x1[0],x1[1],x1[2]);
        }
    }
    return suc;
}


/**************************************************************************
  Implemetation of Kabsch algoritm for finding the best rotation matrix
---------------------------------------------------------------------------
  x    - x(i,m) are coordinates of atom m in set x            (input)
  y    - y(i,m) are coordinates of atom m in set y            (input)
  n    - n is number of atom pairs                            (input)
  mode  - 0:calculate rmsd only                               (input)
          1:calculate rmsd,u,t                                (takes longer)
  rms   - sum of w*(ux+t-y)**2 over all atom pairs            (output)
  u    - u(i,j) is   rotation  matrix for best superposition  (output)
  t    - t(i)   is translation vector for best superposition  (output)
**************************************************************************/
bool RMSDCalculator::Kabsch(vector<Atom*> &_align, vector<Atom*> &_ref, int mode) {
    int i, j, m, m1, l, k;
    real e0, rms1, d, h, g;
    real cth, sth, sqrth, p, det, sigma;
    real xc[3], yc[3];
    real a[3][3], b[3][3], r[3][3], e[3], rr[6], ss[6];
    real sqrt3=1.73205080756888, tol=0.01;
    int ip[]={0, 1, 3, 1, 2, 4, 3, 4, 5};
    int ip2312[]={1, 2, 0, 1};

    int a_failed=0, b_failed=0;
    real epsilon=0.00000001;

    int n=_ref.size();
    if(n != _align.size()) {
        cout << "Two proteins have different length!" << endl;
        return false;
    }

    //initializtation
    _rmsd=0;
    rms1=0;
    e0=0;
    for (i=0; i<3; i++) {
        xc[i]=0.0;
        yc[i]=0.0;
        t[i]=0.0;
        for (j=0; j<3; j++) {
            u[i][j]=0.0;
            r[i][j]=0.0;
            a[i][j]=0.0;
            if (i==j) {
                u[i][j]=1.0;
                a[i][j]=1.0;
            }
        }
    }

    if (n<1) {
        cout << "Protein length is zero!" << endl;
        return false;
    }

    //compute centers for vector sets x, y
    for(i=0; i<n; i++){
        xc[0] += _align[i]->getX();
        xc[1] += _align[i]->getY();
        xc[2] += _align[i]->getZ();

        yc[0] += _ref[i]->getX();
        yc[1] += _ref[i]->getY();
        yc[2] += _ref[i]->getZ();
    }
    for(i=0; i<3; i++){
        xc[i] = xc[i]/(real)n;
        yc[i] = yc[i]/(real)n;
    }

    //compute e0 and matrix r
    for (m=0; m<n; m++) {
        e0 += (_align[m]->getX()-xc[0])*(_align[m]->getX()-xc[0]) \
          +(_ref[m]->getX()-yc[0])*(_ref[m]->getX()-yc[0]);
        e0 += (_align[m]->getY()-xc[1])*(_align[m]->getY()-xc[1]) \
          +(_ref[m]->getY()-yc[1])*(_ref[m]->getY()-yc[1]);
        e0 += (_align[m]->getZ()-xc[2])*(_align[m]->getZ()-xc[2]) \
          +(_ref[m]->getZ()-yc[2])*(_ref[m]->getZ()-yc[2]);
        r[0][0] += (_ref[m]->getX() - yc[0])*(_align[m]->getX() - xc[0]);
        r[0][1] += (_ref[m]->getX() - yc[0])*(_align[m]->getY() - xc[1]);
        r[0][2] += (_ref[m]->getX() - yc[0])*(_align[m]->getZ() - xc[2]);
        r[1][0] += (_ref[m]->getY() - yc[1])*(_align[m]->getX() - xc[0]);
        r[1][1] += (_ref[m]->getY() - yc[1])*(_align[m]->getY() - xc[1]);
        r[1][2] += (_ref[m]->getY() - yc[1])*(_align[m]->getZ() - xc[2]);
        r[2][0] += (_ref[m]->getZ() - yc[2])*(_align[m]->getX() - xc[0]);
        r[2][1] += (_ref[m]->getZ() - yc[2])*(_align[m]->getY() - xc[1]);
        r[2][2] += (_ref[m]->getZ() - yc[2])*(_align[m]->getZ() - xc[2]);
    }
    //compute determinat of matrix r
    det = r[0][0] * ( r[1][1]*r[2][2] - r[1][2]*r[2][1] )       \
    - r[0][1] * ( r[1][0]*r[2][2] - r[1][2]*r[2][0] )       \
    + r[0][2] * ( r[1][0]*r[2][1] - r[1][1]*r[2][0] );
    sigma=det;

    //compute tras(r)*r
    m = 0;
    for (j=0; j<3; j++) {
        for (i=0; i<=j; i++) {
            rr[m]=r[0][i]*r[0][j]+r[1][i]*r[1][j]+r[2][i]*r[2][j];
            m++;
        }
    }

    real spur=(rr[0]+rr[2]+rr[5]) / 3.0;
    real cof = (((((rr[2]*rr[5] - rr[4]*rr[4]) + rr[0]*rr[5]) \
          - rr[3]*rr[3]) + rr[0]*rr[2]) - rr[1]*rr[1]) / 3.0;
    det = det*det;

    for (i=0; i<3; i++){
        e[i]=spur;
    }

    if (spur>0) {
        d = spur*spur;
        h = d - cof;
        g = (spur*cof - det)/2.0 - spur*h;

        if (h>0) {
            sqrth = sqrt(h);
            d = h*h*h - g*g;
            if(d<0.0) d=0.0;
            d = atan2( sqrt(d), -g ) / 3.0;
            cth = sqrth * cos(d);
            sth = sqrth*sqrt3*sin(d);
            e[0]= (spur + cth) + cth;
            e[1]= (spur - cth) + sth;
            e[2]= (spur - cth) - sth;

            if (mode!=0) {//compute a
                for (l=0; l<3; l=l+2) {
                    d = e[l];
                    ss[0] = (d-rr[2]) * (d-rr[5])  - rr[4]*rr[4];
                    ss[1] = (d-rr[5]) * rr[1]      + rr[3]*rr[4];
                    ss[2] = (d-rr[0]) * (d-rr[5])  - rr[3]*rr[3];
                    ss[3] = (d-rr[2]) * rr[3]      + rr[1]*rr[4];
                    ss[4] = (d-rr[0]) * rr[4]      + rr[1]*rr[3];
                    ss[5] = (d-rr[0]) * (d-rr[2])  - rr[1]*rr[1];

                    if (fabs(ss[0])<=epsilon) ss[0]=0.0;
                    if (fabs(ss[1])<=epsilon) ss[1]=0.0;
                    if (fabs(ss[2])<=epsilon) ss[2]=0.0;
                    if (fabs(ss[3])<=epsilon) ss[3]=0.0;
                    if (fabs(ss[4])<=epsilon) ss[4]=0.0;
                    if (fabs(ss[5])<=epsilon) ss[5]=0.0;

                    if (fabs(ss[0]) >= fabs(ss[2])) {
                        j=0;
                        if( fabs(ss[0]) < fabs(ss[5])){
                            j = 2;
                        }
                    } else if ( fabs(ss[2]) >= fabs(ss[5]) ){
                        j = 1;
                    } else {
                        j = 2;
                    }

                    d = 0.0;
                    j = 3 * j;
                    for (i=0; i<3; i++) {
                        k=ip[i+j];
                        a[i][l] = ss[k];
                        d = d + ss[k]*ss[k];
                    }


                    //if( d > 0.0 ) d = 1.0 / sqrt(d);
                    if (d > epsilon) d = 1.0 / sqrt(d);
                    else d=0.0;
                    for (i=0; i<3; i++) {
                        a[i][l] = a[i][l] * d;
                    }
                }//for l

                d = a[0][0]*a[0][2] + a[1][0]*a[1][2] + a[2][0]*a[2][2];
                if ((e[0] - e[1]) > (e[1] - e[2])) {
                    m1=2;
                    m=0;
                } else {
                    m1=0;
                    m=2;
                }
                p=0;
                for(i=0; i<3; i++){
                    a[i][m1] = a[i][m1] - d*a[i][m];
                    p = p + a[i][m1]*a[i][m1];
                }
                if (p <= tol) {
                    p = 1.0;
                    for (i=0; i<3; i++) {
                        if (p < fabs(a[i][m])){
                            continue;
                        }
                        p = fabs( a[i][m] );
                        j = i;
                    }
                    k = ip2312[j];
                    l = ip2312[j+1];
                    p = sqrt( a[k][m]*a[k][m] + a[l][m]*a[l][m] );
                    if (p > tol) {
                        a[j][m1] = 0.0;
                        a[k][m1] = -a[l][m]/p;
                        a[l][m1] =  a[k][m]/p;
                    } else {//goto 40
                        a_failed=1;
                    }
                } else {//if p<=tol
                    p = 1.0 / sqrt(p);
                    for(i=0; i<3; i++){
                        a[i][m1] = a[i][m1]*p;
                    }
                }//else p<=tol
                if (a_failed!=1) {
                    a[0][1] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
                    a[1][1] = a[2][2]*a[0][0] - a[2][0]*a[0][2];
                    a[2][1] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
                }
            }//if(mode!=0)
        }//h>0

        //compute b anyway
        if (mode!=0 && a_failed!=1) {//a is computed correctly
            //compute b
            for (l=0; l<2; l++) {
                d=0.0;
                for(i=0; i<3; i++){
                    b[i][l] = r[i][0]*a[0][l] + r[i][1]*a[1][l] + r[i][2]*a[2][l];
                    d = d + b[i][l]*b[i][l];
                }
                //if( d > 0 ) d = 1.0 / sqrt(d);
                if (d > epsilon) d = 1.0 / sqrt(d);
                else d=0.0;
                for (i=0; i<3; i++) {
                    b[i][l] = b[i][l]*d;
                }
            }
            d = b[0][0]*b[0][1] + b[1][0]*b[1][1] + b[2][0]*b[2][1];
            p=0.0;

            for (i=0; i<3; i++) {
                b[i][1] = b[i][1] - d*b[i][0];
                p += b[i][1]*b[i][1];
            }

            if (p <= tol) {
                p = 1.0;
                for (i=0; i<3; i++) {
                    if (p<fabs(b[i][0])) {
                        continue;
                    }
                    p = fabs( b[i][0] );
                    j=i;
                }
                k = ip2312[j];
                l = ip2312[j+1];
                p = sqrt( b[k][0]*b[k][0] + b[l][0]*b[l][0] );
                if (p > tol) {
                    b[j][1] = 0.0;
                    b[k][1] = -b[l][0]/p;
                    b[l][1] =  b[k][0]/p;
                } else {
                    //goto 40
                    b_failed=1;
                }
            } else {//if( p <= tol )
                p = 1.0 / sqrt(p);
                for(i=0; i<3; i++){
                    b[i][1]=b[i][1]*p;
                }
            }
            if (b_failed!=1){
                b[0][2] = b[1][0]*b[2][1] - b[1][1]*b[2][0];
                b[1][2] = b[2][0]*b[0][1] - b[2][1]*b[0][0];
                b[2][2] = b[0][0]*b[1][1] - b[0][1]*b[1][0];
                //compute u
                for (i=0; i<3; i++){
                    for(j=0; j<3; j++){
                        u[i][j] = b[i][0]*a[j][0] + b[i][1]*a[j][1] \
                                + b[i][2]*a[j][2];
                    }
                }
            }

            //compute t
            for(i=0; i<3; i++){
                t[i] = ((yc[i] - u[i][0]*xc[0]) - u[i][1]*xc[1])    \
                        - u[i][2]*xc[2];
            }
        }//if(mode!=0 && a_failed!=1)
    } else {//spur>0, just compute t and errors
        //compute t
        for (i=0; i<3; i++) {
            t[i] = ((yc[i] - u[i][0]*xc[0]) - u[i][1]*xc[1]) - u[i][2]*xc[2];
        }
    } //else spur>0

    //compute rmsd
    for(i=0; i<3; i++){
        if( e[i] < 0 ) e[i] = 0;
        e[i] = sqrt( e[i] );
    }
    d = e[2];
    if( sigma < 0.0 ){
        d = - d;
    }
    d = (d + e[1]) + e[0];
    rms1 = (e0 - d) - d;
    if( rms1 < 0.0 ) rms1 = 0.0;

    _rmsd=sqrt(rms1/(real)n);

    return true;
}

real RMSDCalculator::rmsd(vector<Atom*>& A, vector<Atom*>& B) {
  if (A.size() != B.size())
    MstUtils::error("atom vectors of different length (" + MstUtils::toString(A.size()) + " and " + MstUtils::toString(B.size()) + ")", "RMSDCalculator::rmsd(vector<Atom*>&, vector<Atom*>&)");

  real ret = 0;
  for (int i = 0; i < A.size(); i++) {
    ret += A[i]->distance2(B[i]);
  }
  return sqrt(ret/A.size());
}

real RMSDCalculator::rmsd(Structure& A, Structure& B) {
  vector<Atom*> atomsA = A.getAtoms();
  vector<Atom*> atomsB = B.getAtoms();
  return rmsd(atomsA, atomsB);
}


/* --------- ProximitySearch --------- */

ProximitySearch::ProximitySearch(real _xlo, real _ylo, real _zlo, real _xhi, real _yhi, real _zhi, int _N) {
  xlo = _xlo; ylo = _ylo; zlo = _zlo;
  xhi = _xhi; yhi = _yhi; zhi = _zhi;
  reinitBuckets(_N);
  setBinWidths();
}

ProximitySearch::ProximitySearch(AtomPointerVector& _atoms, int _N, bool _addAtoms, vector<int>* tags, real pad) {
  calculateExtent(_atoms);
  xlo -= pad; ylo -= pad; zlo -= pad;
  xhi += pad; yhi += pad; zhi += pad;
  reinitBuckets(_N);
  setBinWidths();
  if (_addAtoms) {
    for (int i = 0; i < _atoms.size(); i++) {
      addPoint(_atoms[i]->getCoor(), (tags == NULL) ? i : (*tags)[i]);
    }
  }
}

ProximitySearch::ProximitySearch(AtomPointerVector& _atoms, real _characteristicDistance, bool _addAtoms, vector<int>* tags, real pad) {
  calculateExtent(_atoms);
  if (xlo == xhi) { xlo -= _characteristicDistance/2; xhi += _characteristicDistance/2; }
  if (ylo == yhi) { ylo -= _characteristicDistance/2; yhi += _characteristicDistance/2; }
  if (zlo == zhi) { zlo -= _characteristicDistance/2; zhi += _characteristicDistance/2; }
  xlo -= pad; ylo -= pad; zlo -= pad;
  xhi += pad; yhi += pad; zhi += pad;
  int _N = int(ceil(max(max((xhi - xlo), (yhi - ylo)), (zhi - zlo))/_characteristicDistance));
  reinitBuckets(_N);
  setBinWidths();
  if (_addAtoms) {
    for (int i = 0; i < _atoms.size(); i++) {
      addPoint(_atoms[i]->getCoor(), (tags == NULL) ? i : (*tags)[i]);
    }
  }
}

void ProximitySearch::setBinWidths() {
  xbw = (xhi - xlo)/(N - 1);
  ybw = (yhi - ylo)/(N - 1);
  zbw = (zhi - zlo)/(N - 1);
}

ProximitySearch::~ProximitySearch() {
  for (int i = 0; i < pointList.size(); i++) delete(pointList[i]);
}

void ProximitySearch::calculateExtent(Structure& S, real& _xlo, real& _ylo, real& _zlo, real& _xhi, real& _yhi, real& _zhi) {
  AtomPointerVector atoms = S.getAtoms();
  calculateExtent(atoms, _xlo, _ylo, _zlo, _xhi, _yhi, _zhi);
}

void ProximitySearch::calculateExtent(AtomPointerVector& _atoms, real& _xlo, real& _ylo, real& _zlo, real& _xhi, real& _yhi, real& _zhi) {
  if (_atoms.size() == 0) { cout << "Error in nnclass::calculateExtent() -- empty atom vector passed!\n"; exit(-1); }
  _xlo = _xhi = _atoms[0]->getX();
  _ylo = _yhi = _atoms[0]->getY();
  _zlo = _zhi = _atoms[0]->getZ();
  for (int i = 0; i < _atoms.size(); i++) {
    if (_xlo > _atoms[i]->getX()) _xlo = _atoms[i]->getX();
    if (_ylo > _atoms[i]->getY()) _ylo = _atoms[i]->getY();
    if (_zlo > _atoms[i]->getZ()) _zlo = _atoms[i]->getZ();
    if (_xhi < _atoms[i]->getX()) _xhi = _atoms[i]->getX();
    if (_yhi < _atoms[i]->getY()) _yhi = _atoms[i]->getY();
    if (_zhi < _atoms[i]->getZ()) _zhi = _atoms[i]->getZ();
  }
}

void ProximitySearch::reinitBuckets(int _N) {
  N = _N;
  buckets.resize(N);
  for (int i = 0; i < N; i++) {
    buckets[i].resize(N);
    for (int j = 0; j < N; j++) {
      buckets[i][j].resize(N);
      for (int k = 0; k < N; k++) { buckets[i][j][k].resize(0); }
    }
  }
  pointList.resize(0);
  pointTags.resize(0);
}

void ProximitySearch::addPoint(CartesianPoint _p, int tag) {
  int i, j, k;
  pointBucket(&_p, &i, &j, &k);
  if ((i < 0) || (j < 0) || (k < 0) || (i > N-1) || (j > N-1) || (k > N-1)) { cout << "Error: point " << _p << " out of range for ProximitySearch object!\n"; exit(-1); }
  CartesianPoint* p = new CartesianPoint(_p);
  buckets[i][j][k].push_back(pointList.size());
  pointList.push_back(p);
  pointTags.push_back(tag);
}

void ProximitySearch::addAtoms(AtomPointerVector& apv, vector<int>* tags) {
  if ((tags != NULL) && (apv.size() != tags->size())) MstUtils::error("different number of atoms and tags specified!", "ProximitySearch::addAtoms");
  for (int i = 0; i < apv.size(); i++) {
    addPoint(apv[i], (tags == NULL) ? i : (*tags)[i]);
  }
}

bool ProximitySearch::isPointWithinGrid(CartesianPoint _p) {
  int i, j, k;
  pointBucket(&_p, &i, &j, &k);
  if ((i < 0) || (j < 0) || (k < 0) || (i > N-1) || (j > N-1) || (k > N-1)) return false;
  return true;
}

void ProximitySearch::pointBucket(real px, real py, real pz, int* i, int* j, int* k) {
  *i = (int) floor((px - xlo)/xbw + 0.5);
  *j = (int) floor((py - ylo)/ybw + 0.5);
  *k = (int) floor((pz - zlo)/zbw + 0.5);
}

void ProximitySearch::limitIndex(int *ind) {
  if (*ind < 0) *ind = 0;
  if (*ind > N-1) *ind = N-1;
}

bool ProximitySearch::pointsWithin(CartesianPoint c, real dmin, real dmax, vector<int>* list, bool byTag) {
  real cx = c.getX(); real cy = c.getY(); real cz = c.getZ();
  // first check if the point is outside of the bounding box of the point cloud by a sufficient amount
  if ((cx < xlo - dmax) || (cy < ylo - dmax) || (cz < zlo - dmax) || (cx > xhi + dmax) || (cy > yhi + dmax) || (cz > zhi + dmax)) return false;

  real d2, dmin2, dmax2;
  int ci, cj, ck;
  int imax1, jmax1, kmax1, imax2, jmax2, kmax2; // external box (no point in looking beyond it, points there are too far)
  int imin1, jmin1, kmin1, imin2, jmin2, kmin2; // internal box (no point in looking within it, points there are too close)
  pointBucket(cx, cy, cz, &ci, &cj, &ck);
  pointBucket(cx - dmax, cy - dmax, cz - dmax, &imax1, &jmax1, &kmax1);
  pointBucket(cx + dmax, cy + dmax, cz + dmax, &imax2, &jmax2, &kmax2);
  if (dmin > 0) {
    real sr3 = sqrt(3);
    pointBucket(cx - dmin/sr3, cy - dmin/sr3, cz - dmin/sr3, &imin1, &jmin1, &kmin1);
    pointBucket(cx + dmin/sr3, cy + dmin/sr3, cz + dmin/sr3, &imin2, &jmin2, &kmin2);
    // need to trim the internal box to make sure it is fully contained within the sphere of radius dmin from the central point
    if (imin1 != ci) imin1++;
    if (jmin1 != cj) jmin1++;
    if (kmin1 != ck) kmin1++;
    if (imin2 != ci) imin2--;
    if (jmin2 != cj) jmin2--;
    if (kmin2 != ck) kmin2--;
  } else {
    imin1 = imin2 = ci;
    jmin1 = jmin2 = cj;
    kmin1 = kmin2 = ck;
  }
  limitIndex(&imin1); limitIndex(&imin2); limitIndex(&jmin1); limitIndex(&jmin2); limitIndex(&kmin1); limitIndex(&kmin2);
  limitIndex(&imax1); limitIndex(&imax2); limitIndex(&jmax1); limitIndex(&jmax2); limitIndex(&kmax1); limitIndex(&kmax2);

  // search only within the boxes where points of interest can be, in principle
  if (list != NULL) list->clear();
  bool found = false;
  bool yesno = (list == NULL);
  dmin2 = dmin*dmin; dmax2 = dmax*dmax;
  for (int i = imax1; i <= imax2; i++) {
    bool insi = (i >= imin1) && (i <= imin2);
    vector<vector<vector<int> > >& Bi = buckets[i];
    for (int j = jmax1; j <= jmax2; j++) {
      bool ins = insi && (j >= jmin1) && (j <= jmin2);
      vector<vector<int> >& Bij = Bi[j];
      for (int k = kmax1; k <= kmax2; k++) {
        vector<int>& Bijk = Bij[k];
        // check all points in bucket i, j, k
        for (int ii = 0; ii < Bijk.size(); ii++) {
          int pi = Bijk[ii];
          d2 = c.distance2(*(pointList[pi]));
          if ((d2 >= dmin2) && (d2 <= dmax2)) {
            if (yesno) return true;
            list->push_back(byTag ? pointTags[Bijk[ii]] : Bijk[ii]);
            found = true;
          }
        }
        // skip the range from kmin1 to kmin2 (too close)
        if (ins && (k == kmin1) && (kmin1 != kmin2)) k = kmin2 - 1;
      }
    }
  }
  return found;
}

vector<int> ProximitySearch::getPointsWithin(CartesianPoint c, real dmin, real dmax, bool byTag) {
  vector<int> closeOnes;
  pointsWithin(c, dmin, dmax, &closeOnes, byTag);
  return closeOnes;
}

bool ProximitySearch::overlaps(ProximitySearch& other, real pad) {
  if ((other.xlo > xhi + pad) || (xlo > other.xhi + pad) || (other.ylo > yhi + pad) || (ylo > other.yhi + pad) || (other.zlo > zhi + pad) || (zlo > other.zhi + pad)) {
    return false;
  }
  return true;
}

/* --------- MstUtils --------- */
void MstUtils::openFile (fstream& fs, string filename, ios_base::openmode mode, string from) {
  fs.open(filename.c_str(), mode);
  if (!fs.is_open()) {
    if (!from.empty()) from += " -> ";
    MstUtils::error("could not open file '" + filename + "'", from + "MstUtils::openFile");
  }
}

void MstUtils::fileToArray(string _filename, vector<string>& lines) {
  fstream inp;
  MstUtils::openFile(inp, _filename, ios_base::in, "MstUtils::fileToArray");
  string line;
  while (true) {
    getline(inp, line);
    if (inp.eof()) break; // if eof set upon trying to read the line, then it was not really a valid line, so we are done
    lines.push_back(line);
  }
}

string MstUtils::pathBase(string fn) {
  if (fn.find_last_of(".") == string::npos) return fn;
  else return fn.substr(0, fn.find_last_of("."));
}

bool MstUtils::fileExists(const char *filename) {
  struct stat buffer ;
  if (stat( filename, &buffer) == 0) return true;
  return false;
}

bool MstUtils::isDir(const char *filename) {
  struct stat buffer ;
  if (stat( filename, &buffer) < 0) return false;
  return (buffer.st_mode & S_IFDIR);
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
  if ((sscanf(num.c_str(), "%lf", &ret) != 1) && strict) MstUtils::error("failed to convert '" + num + "' to real", "MstUtils::toReal");
  return (real) ret;
}

MST::real MstUtils::mod(MST::real num, MST::real den) {
  return num - ((MST::real) floor((double) num / (double) den)) * den;
}

string MstUtils::readNullTerminatedString(fstream& ifs) {
  string str;
  char c;
  while (ifs.get(c)) {
    if (c == '\0') break;
    str += c;
  }
  return str;
}
