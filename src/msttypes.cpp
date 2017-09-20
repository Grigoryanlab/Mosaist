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

Structure::Structure(const Structure& S) {
  copy(S);
}

void Structure::copy(const Structure& S) {
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

Structure::Structure(const vector<Atom*>& atoms) {
  numResidues = numAtoms = 0;
  addAtoms(atoms);
}

Structure::Structure(const vector<Residue*>& residues) {
  numResidues = numAtoms = 0;
  for (int i = 0; i < residues.size(); i++) addResidue(residues[i]);
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

Structure& Structure::operator=(const Structure& A) {
  reset();
  copy(A);
  return *this;
}

void Structure::readPDB(string pdbFile, string options) {
  name = pdbFile;
  int lastresnum = -999999;
  string lastresname = "XXXXXX";
  string lasticode = "";
  string lastchainID = "";
  string lastalt = " ";
  Chain* chain = NULL;
  Residue* residue = NULL;

  // various parsing options (the wonders of dealing with the good-old PDB format)
  bool ter = true;                   // flag to indicate that chain terminus was reached. Initialize to true so as to create a new chain upon reading the first atom.
  bool usesegid = false;             // use segment IDs to name chains instead of chain IDs? (useful when the latter are absent OR when too many chains, so need multi-letter names)
  bool skipHetero = false;           // skip hetero-atoms?
  bool charmmFormat = false;         // the PDB file was written by CHARMM (slightly different format)
  bool charmm19Format = false;       // upon reading, convert from all-hydrogen topology (param22 and higher) to the CHARMM19 united-atom topology (matters for HIS protonation states)
  bool fixIleCD1 = true;             // rename CD1 in ILE to CD (as is standard in MM packages)
  bool iCodesAsSepResidues = true;   // consequtive residues that differ only in their insertion code will be treated as separate residues
  bool uniqChainIDs = true;          // make sure chain IDs are unique, even if they are not unique in the read file
  bool ignoreTER = false;            // if true, will not pay attention to TER lines in deciding when chains end/begin
  bool verbose = true;               // report various warnings when weird things are found and fixed?

  // user-specified custom parsing options
  options = MstUtils::uc(options);
  if (options.find("USESEGID") != string::npos) usesegid = true;
  if (options.find("SKIPHETERO") != string::npos) skipHetero = true;
  if (options.find("CHARMM") != string::npos) charmmFormat = true;
  if (options.find("CHARMM19") != string::npos) charmm19Format = true;
  if (options.find("ALLOW DUPLICATE CIDS") != string::npos) uniqChainIDs = false;
  if (options.find("ALLOW ILE CD1") != string::npos) fixIleCD1 = false;
  if (options.find("IGNORE-ICODES") != string::npos) iCodesAsSepResidues = true;
  if (options.find("IGNORE-TER") != string::npos) ignoreTER = true;
  if (options.find("QUIET") != string::npos) verbose = false;

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
    mstreal x = MstUtils::toReal(line.substr(30, 8));
    mstreal y = MstUtils::toReal(line.substr(38, 8));
    mstreal z = MstUtils::toReal(line.substr(46, 8));
    string segID = MstUtils::trim(line.substr(72, 4));
    mstreal B = MstUtils::toReal(line.substr(60, 6));
    mstreal occ = MstUtils::toReal(line.substr(54, 6));
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
      if (verbose && chainID.compare(chain->getID())) {
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
    if ((resnum != lastresnum) || resname.compare(lastresname) || (iCodesAsSepResidues && (icode.compare(lasticode))))  {
      // this corresponds to a case, where the alternative location flag is being used to
      // designate two (or more) different possible amino acids at a particular position
      // (e.g., where the density is not clear to assign one). In this case, we shall keep
      // only the first option, because we don't know any better.
      if ((resnum == lastresnum) && resname.compare(lastresname) && (alt != lastalt)) {
        continue;
      }
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
    lastalt = alt;
  }
  ifh.close();
}

void Structure::writePDB(string pdbFile, string options) const {
  fstream ofs; MstUtils::openFile(ofs, pdbFile, fstream::out, "Structure::writePDB(string, string)");
  writePDB(ofs, options);
  ofs.close();
}

void Structure::writePDB(fstream& ofs, string options) const {
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
    string goodNames = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
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

void Structure::deleteChain(Chain* chain) {
  // make sure the chain is from this Structure
  int i = 0;
  for (i = 0; i < chains.size(); i++) {
    if (chains[i] == chain) break;
  }
  if (i == chains.size()) MstUtils::error("chain not from this structure", "Structure::deleteChain");

  chains.erase(chains.begin() + i);
  numResidues -= chain->residueSize();
  numAtoms -= chain->atomSize();
  if (chainsByID.find(chain->getID()) != chainsByID.end()) chainsByID.erase(chain->getID());
  if (chainsBySegID.find(chain->getSegID()) != chainsBySegID.end()) chainsBySegID.erase(chain->getSegID());
  delete chain;
}

Residue& Structure::getResidue(int i) const {
  if ((i < 0) && (i >= residueSize()))
    MstUtils::error("residue index " + MstUtils::toString(i) + " out of range for Structure", "Structure::getResidue(int)");
  int io = i;
  for (int ci = 0; ci < chainSize(); ci++) {
    Chain& chain = getChain(ci);
    if (i >= chain.residueSize()) {
      i -= chain.residueSize();
    } else {
      return chain[i];
    }
  }
  MstUtils::error("something strange happened when fetching residue " + MstUtils::toString(io) + " from Structure object that reports " + MstUtils::toString(this->residueSize()) + " residues; most likely, various counters are inconsistent in Structure object", "Structure::getResidue(int)");
  return *(new Residue()); // just to make the compiler happy and not throw a warning; this is never reached
}

vector<Atom*> Structure::getAtoms() const {
  vector<Atom*> vec(this->atomSize());
  int ii = 0;
  for (int i = 0; i < this->chainSize(); i++) {
    Chain& c = (*this)[i];
    for (int j = 0; j < c.residueSize(); j++) {
      Residue& r = c[j];
      for (int k = 0; k < r.atomSize(); k++) {
        vec[ii] = &(r[k]);
        ii++;
      }
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

Structure Structure::reassignChainsByConnectivity(mstreal maxPeptideBond) {
  Structure S;
  reassignChainsByConnectivity(S, maxPeptideBond);
  return S;
}

void Structure::reassignChainsByConnectivity(Structure& dest, mstreal maxPeptideBond) {
  vector<Residue*> residues = this->getResidues();
	Chain* chain = dest.appendChain("A");
	for (int i = 0; i < residues.size() - 1; i++) {
    chain->appendResidue(new Residue(*residues[i]));
		Atom* atomC = residues[i]->findAtom("C", true);
		Atom* atomN = residues[i + 1]->findAtom("N", true);
    if ((atomC == NULL) || (atomN == NULL)) MstUtils::error("cannot break into disjoint segments as some C or N backbone atoms are missing", "Structure::reassignChainsByConnectivity");
		if (atomC->distance(atomN) > maxPeptideBond) {
      chain = dest.appendChain("A");
    }
	}
	chain->appendResidue(new Residue(*residues[residues.size() - 1]));
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
  newResidue = newChain->findResidue(oldResidue->getName(), oldResidue->getNum(), oldResidue->getIcode());
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

Residue* Structure::addResidue(Residue* res) {
  if (res->getParent() == NULL) MstUtils::error("cannot add a disembodied Residue", "Structure::addResidue");
  Chain* oldChain = res->getParent();

  // is there a chain matching the Residue's parent chain? If not, create one.
  Chain* newChain = getChainByID(oldChain->getID());
  if (newChain == NULL) {
    newChain = new Chain(oldChain->getID(), oldChain->getSegID());
    this->appendChain(newChain);
  }

  // append a copy of the given residue into the correct chain
  Residue* newResidue = new Residue(*res);
  newChain->appendResidue(newResidue);
  return newResidue;
}

int Structure::getResidueIndex(Residue* res) {
  Chain* parentChain = res->getParent();
  if (parentChain == NULL)
    MstUtils::error("cannot find index of a disembodied residue '" + MstUtils::toString(res) + "'", "Structure::getResidueIndex(Residue*)");

  int n = 0;
  bool found = false;
  for (int i = 0; i < chainSize(); i++) {
    Chain& chain = (*this)[i];
    if (&chain  == parentChain) {
      n += chain.getResidueIndex(res);
      found = true;
      break;
    }
    n += chain.residueSize();
  }
  if (!found) MstUtils::error("residue not from Structure '" + MstUtils::toString(res) + "'", "Structure::getResidueIndex(Residue*)");

  return n;
}

/* --------- Chain --------- */
Chain::Chain() {
  numAtoms = 0;
  parent = NULL;
}

Chain::Chain(const Chain& C) {
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

Chain::Chain(const string& chainID, const string& segID) {
  numAtoms = 0;
  parent = NULL;
  cid = chainID;
  sid = segID;
}

Chain::Chain(const string& chainID, const string& segID, const vector<Residue*>& residues) {
	numAtoms = 0;
	parent = NULL;
	cid = chainID;
	sid = segID;
	appendResidueCopies(residues);
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

int Chain::getResidueIndex(const Residue* res) {
  if (residueIndexInChain.find((Residue*) res) == residueIndexInChain.end())
    MstUtils::error("passed residue does not appear in chain's index map", "Chain::residueIndexInChain");
  return residueIndexInChain[(Residue*) res];
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

void Chain::appendResidueCopies(const vector<Residue*>& residues) {
	for (int i = 0; i < residues.size(); i++) {
		insertResidueCopy(residues[i], -1);
	}
}

void Chain::insertResidue(Residue* R, int index) {
  if ((index < 0) || (index > residues.size()))
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

const mstreal Residue::badDihedral = 999.0;

Residue::Residue() {
  resname = "UNK";
  resnum = 1;
  parent = NULL;
  icode = ' ';
}

Residue::Residue(const Residue& R, bool copyAlt) {
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
  delete atoms[i];
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

void Residue::replaceAtoms(const vector<Atom*>& newAtoms, vector<int>* toRemove) {
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
mstreal Residue::getPhi(bool strict) {
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
mstreal Residue::getPsi(bool strict) {
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
mstreal Residue::getOmega(bool strict) {
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

int Residue::getResidueIndex() const {
  Chain* parentChain = NULL; Structure* parentStructure = NULL;
  parentChain = getParent();
  if (parentChain != NULL) parentStructure = parentChain->getParent();
  if ((parentChain == NULL) || (parentStructure == NULL))
    MstUtils::error("cannot find index of a disembodied residue '" + MstUtils::toString(this) + "'", "Residue::getResidueIndex()");

  int n = 0;
  bool found = false;
  for (int i = 0; i < parentStructure->chainSize(); i++) {
    Chain& chain = (*parentStructure)[i];
    if (&chain == parentChain) {
      n += chain.getResidueIndex(this);
      found = true;
      break;
    }
    n += chain.residueSize();
  }
  if (!found) MstUtils::error("residue not in its parent Structure '" + MstUtils::toString(this) + "'", "Residue::getResidueIndex()");

  return n;
}

/* --------- Atom --------- */
Atom::Atom() {
  parent = NULL;
  het = false;
  name = MstUtils::copyStringC("UNK");
  setName("");
  alternatives = NULL;
}

Atom::Atom(const Atom& A, bool copyAlt) {
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

Atom::Atom(int _index, string _name, mstreal _x, mstreal _y, mstreal _z, mstreal _B, mstreal _occ, bool _het, char _alt, Residue* _parent) {
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

mstreal& Atom::operator[](int i) {
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
  return x; // just to silence the warning from some compilres; in reality, this is never reached
}

CartesianPoint Atom::getCoor() const {
  CartesianPoint coor(x, y, z); return coor;
}

CartesianPoint Atom::getAltCoor(int altInd) const {
  if ((alternatives == NULL) || (altInd >= alternatives->size()) || (altInd < 0)) MstUtils::error("alternative index " + MstUtils::toString(altInd) + " out of bounds (" + MstUtils::toString(alternatives->size()) + " alternatives available)", "Atom::swapWithAlternative");
  altInfo& targ = (*alternatives)[altInd];
  CartesianPoint coor(targ.x, targ.y, targ.z);
  return coor;
}

void Atom::setCoor(const CartesianPoint& xyz) {
  x = xyz[0]; y = xyz[1]; z = xyz[2];
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

void Atom::addAlternative(mstreal _x, mstreal _y, mstreal _z, mstreal _B, mstreal _occ, char _alt) {
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

mstreal Atom::distance(const Atom& another) const {
  return sqrt((x - another.x)*(x - another.x) + (y - another.y)*(y - another.y) + (z - another.z)*(z - another.z));
}

mstreal Atom::distance2(const Atom& another) const {
  return (x - another.x)*(x - another.x) + (y - another.y)*(y - another.y) + (z - another.z)*(z - another.z);
}

mstreal Atom::angle(const Atom& A, const Atom& B, bool radians) const {
  return CartesianGeometry::angle(*this, A, B, radians);
}

mstreal Atom::angle(const Atom* A, const Atom* B, bool radians) const {
  return CartesianGeometry::angle(this, A, B, radians);
}

mstreal Atom::dihedral(const Atom& A, const Atom& B, const Atom& C, bool radians) const {
  return CartesianGeometry::dihedral(*this, A, B, C, radians);
}

mstreal Atom::dihedral(const Atom* A, const Atom* B, const Atom* C, bool radians) const {
  return CartesianGeometry::dihedral(this, A, B, C, radians);
}

bool Atom::build(const Atom& diA, const Atom& anA, const Atom& thA, mstreal di, mstreal an, mstreal th, bool radians) {
  if (!radians) {
    an *= M_PI/180;
    th *= M_PI/180;
  }

  // unit vector from diA to anA (B - C)
  CartesianPoint uCB = ((CartesianPoint) diA - (CartesianPoint) anA).getUnit();

  // vector from anA to thA (C - D)
  CartesianPoint dDC = (CartesianPoint) anA - (CartesianPoint) thA;

  mstreal an2 = M_PI - an;
  mstreal th2 = M_PI + th;
  mstreal rsin = di * sin(an2);
  mstreal rcos = di * cos(an2);
  mstreal rsinsin = rsin * sin(th2);
  mstreal rsincos = rsin * cos(th2);

  CartesianPoint c1 = uCB.cross(dDC);

  // when the first three atoms of the dihedral are co-linear, can't interpret
  // the dihedral angle, so just place the atom di distance away from diA, but
  // along some arbitrary direction
  if (c1.norm() < 0.00000001) { setCoor(diA.getCoor() + CartesianPoint(di, 0, 0)); return false; }
  c1 *= rsinsin / c1.norm();
  CartesianPoint c2 = (-uCB * dDC.dot(uCB) + dDC).getUnit() * rsincos;
  CartesianPoint dd = uCB * rcos + c1 + c2;

  // set coordinate of placed atom
  setCoor(diA.getCoor() + dd);
  return true;
}

/* --------- AtomPointerVector --------- */

void AtomPointerVector::push_back(const Residue& R) {
  int sz = this->size();
  this->resize(sz + R.atomSize());
  for (int i = 0; i < R.atomSize(); i++) {
    (*this)[i + sz] = &(R[i]);
  }
}

CartesianPoint AtomPointerVector::getGeometricCenter() {
  CartesianPoint C(3, 0);
  for (int i = 0; i < this->size(); i++) {
    C += CartesianPoint(*((*this)[i]));
  }
  C /= this->size();
  return C;
}

void AtomPointerVector::getGeometricCenter(mstreal& xc, mstreal& yc, mstreal& zc) {
  xc = 0; yc = 0; zc = 0;
  Atom* a;
  for (int i = 0; i < this->size(); i++) {
    a = (*this)[i];
    xc += (*a)[0]; yc += (*a)[1]; zc += (*a)[2];
  }
  xc /= size(); yc /= size(); zc /= size();
}

void AtomPointerVector::center() {
  CartesianPoint C = getGeometricCenter();
  for (int i = 0; i < this->size(); i++) {
    Atom& a = *((*this)[i]);
    for (int k = 0; k < 3; k++) a[k] -= C[k];
  }
}

mstreal AtomPointerVector::radiusOfGyration() {
  CartesianPoint center = getGeometricCenter();
  mstreal s = 0;
  for (int i = 0; i < size(); i++) {
    s += CartesianPoint((*this)[i]).distance2(center);
  }
  return sqrt(s / size());
}

mstreal AtomPointerVector::boundingSphereRadiusCent() {
  CartesianPoint center = getGeometricCenter();
  mstreal r = 0;
  for (int i = 0; i < size(); i++) {
    mstreal dist = CartesianPoint((*this)[i]).distance(center);
    if (dist > r) r = dist;
  }
  return r;
}

void AtomPointerVector::deletePointers() {
  for (int i = 0; i < size(); i++) delete (*this)[i];
  resize(0);
}

AtomPointerVector AtomPointerVector::clone() {
  AtomPointerVector into;
  clone(into);
  return into;
}

AtomPointerVector AtomPointerVector::subvector(int beg, int end) {
  return AtomPointerVector(vector<Atom*>(begin() + beg, begin() + end));
}

void AtomPointerVector::clone(AtomPointerVector& into) {
  int L = into.size();
  into.resize(L + size());
  for (int i = 0; i < size(); i++) {
    Atom* newAtom = new Atom(*((*this)[i]));
    newAtom->setParent(NULL);
    into[L + i] = newAtom;
  }
}

/* --------- CartesianPoint --------- */

CartesianPoint::CartesianPoint(const Atom& A) {
  this->resize(3, 0);
  (*this)[0] = A.getX();
  (*this)[1] = A.getY();
  (*this)[2] = A.getZ();
}

CartesianPoint& CartesianPoint::operator+=(const CartesianPoint &rhs) {
  if (size() != rhs.size()) MstUtils::error("points of different dimensionality!", "CartesianPoint::operator+=");
  for (int i = 0; i < size(); i++) {
    (*this)[i] += rhs[i];
  }
  return *this;
}

CartesianPoint& CartesianPoint::operator-=(const CartesianPoint &rhs) {
  if (size() != rhs.size()) MstUtils::error("points of different dimensionality!", "CartesianPoint::operator-=");
  for (int i = 0; i < size(); i++) {
    (*this)[i] -= rhs[i];
  }
  return *this;
}

CartesianPoint& CartesianPoint::operator*=(const mstreal& s) {
  for (int i = 0; i < size(); i++) {
    (*this)[i] *= s;
  }
  return *this;
}

CartesianPoint& CartesianPoint::operator/=(const mstreal& s) {
  for (int i = 0; i < size(); i++) {
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
  return CartesianPoint(size(), 0) - *this;
}

const CartesianPoint CartesianPoint::operator*(const mstreal& s) const {
  CartesianPoint result = *this;
  result *= s;
  return result;
}

const CartesianPoint CartesianPoint::operator/(const mstreal& s) const {
  CartesianPoint result = *this;
  result /= s;
  return result;
}

CartesianPoint& CartesianPoint::operator=(const Atom& A) {
  *this = CartesianPoint(A);
  return *this;
}

mstreal CartesianPoint::norm() const {
  mstreal n = 0;
  for (int i = 0; i < size(); i++) n += (*this)[i]*(*this)[i];
  return sqrt(n);
}

mstreal CartesianPoint::norm2() const {
  mstreal n = 0;
  for (int i = 0; i < size(); i++) n += (*this)[i]*(*this)[i];
  return n;
}

mstreal CartesianPoint::mean() const {
  return sum()/size();
}

mstreal CartesianPoint::sum() const {
  mstreal s = 0;
  for (int i = 0; i < size(); i++) s += (*this)[i];
  return s;
}

mstreal CartesianPoint::median() const {
  if (size() == 0) return 0; // could also throw an error in this case
  vector<mstreal> vec = *this;
  sort(vec.begin(), vec.end());
  if (vec.size() % 2 == 0) return (vec[vec.size() / 2] + vec[vec.size() / 2 - 1])/2;
  return vec[vec.size() / 2];
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

mstreal CartesianPoint::dot(CartesianPoint other) const {
  if (size() != other.size()) MstUtils::error("vector size mismatch", "CartesianPoint::dot");

  mstreal d = 0;
  for (int i = 0; i < size(); i++) {
    d += (*this)[i] * other[i];
  }

  return d;
}

mstreal CartesianPoint::distance(const CartesianPoint& another) const {
  if (this->size() != another.size()) MstUtils::error("point dimensions disagree", "CartesianPoint::distance(CartesianPoint&)");
  mstreal d = 0;
  for (int i = 0; i < this->size(); i++) {
    d += ((*this)[i] - another[i])*((*this)[i] - another[i]);
  }
  return sqrt(d);
}

mstreal CartesianPoint::distancenc(const CartesianPoint& another) const {
  mstreal d = 0, dd;
  const CartesianPoint& p = *this;
  for (int i = 0; i < p.size(); i++) {
    dd = p[i] - another[i];
    d += dd*dd;
  }
  return sqrt(d);
}

mstreal CartesianPoint::distance2(const CartesianPoint& another) const {
  if (this->size() != another.size()) MstUtils::error("point dimensions disagree", "CartesianPoint::distance2(CartesianPoint&)");
  mstreal d = 0, dd;
  for (int i = 0; i < this->size(); i++) {
    dd = (*this)[i] - another[i];
    d += dd*dd;
  }
  return d;
}

mstreal CartesianPoint::distance2nc(const CartesianPoint& another) const {
  mstreal d = 0, dd;
  const CartesianPoint& p = *this;
  for (int i = 0; i < p.size(); i++) {
    dd = p[i] - another[i];
    d += dd*dd;
  }
  return d;
}

/* --------- CartesianGeometry --------- */
mstreal CartesianGeometry::dihedral(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4, bool radians) {
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
  if (!radians) angle *= 180/M_PI;
  return angle;
}

mstreal CartesianGeometry::dihedral(const CartesianPoint * _p1, const CartesianPoint * _p2, const CartesianPoint * _p3, const CartesianPoint * _p4, bool radians) {
  return dihedral(*_p1, *_p2, *_p3, *_p4, radians);
}

mstreal CartesianGeometry::angle(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, bool radians) {
  CartesianPoint v21 = (_p1 - _p2).getUnit();
  CartesianPoint v23 = (_p3 - _p2).getUnit();
  mstreal c = v21.dot(v23);
  return atan2(sqrt(1 - c*c), c) * (radians ? 1 : 180/M_PI);
}

mstreal CartesianGeometry::angle(const CartesianPoint * _p1, const CartesianPoint * _p2, const CartesianPoint * _p3, bool radians) {
  return angle(*_p1, *_p2, *_p3, radians);
}

mstreal CartesianGeometry::angleDiff(mstreal A, mstreal B, bool radians) {
  mstreal PI = radians ? M_PI : 180.0;
  mstreal TWOPI = 2*PI;
  mstreal da = MstUtils::mod((MstUtils::mod(A, TWOPI) - MstUtils::mod(B, TWOPI)), TWOPI);
  if (da > PI) da -= TWOPI;
  return da;
}

mstreal CartesianGeometry::angleDiffCCW(mstreal A, mstreal B, bool radians) {
  mstreal TWOPI = radians ? 2*M_PI : 360.0;
  mstreal da = MstUtils::mod((MstUtils::mod(A, TWOPI) - MstUtils::mod(B, TWOPI)), TWOPI);
  return da;
}

pair<mstreal, mstreal> CartesianGeometry::angleRange(const vector<mstreal>& angles, bool radians) {
  int bestMin = 0; int bestMax = 0; double bestArc = (radians ? 2*M_PI : 360.0) + 1;
  for (int minInd = 0; minInd < angles.size(); minInd++) {
    int maxInd = 0; double arc = 0;
    for (int j = 0; j < angles.size(); j++) {
      if (CartesianGeometry::angleDiffCCW(angles[j], angles[minInd]) > arc) {
        arc = CartesianGeometry::angleDiffCCW(angles[j], angles[minInd]);
        maxInd = j;
      }
    }
    if (arc < bestArc) {
      bestMin = minInd; bestMax = maxInd; bestArc = arc;
    }
  }
  return pair<mstreal, mstreal>(angles[bestMin], angles[bestMax]);
}

mstreal CartesianGeometry::angleMean(const vector<mstreal>& angles, bool radians) {
  pair<mstreal, mstreal> minmax = angleRange(angles, radians);
  mstreal m = 0;
  for (int i = 0; i < angles.size(); i++) {
    m += angleDiffCCW(angles[i], minmax.first);
  }
  m /= angles.size();
  return angleDiff(minmax.first + m, 0, radians);
}

/* --------- selector --------------- */

selector::selector(Structure& S) {
  atoms = S.getAtoms();
  residues.resize(atoms.size());
  chains.resize(atoms.size());
  for (int i = 0; i < atoms.size(); i++) {
    residues[i] = atoms[i]->getParent();
    if (residues[i] != NULL) chains[i] = residues[i]->getParent();
    if ((residues[i] == NULL) || (chains[i] == NULL)) MstUtils::error("internally inconsistent Structure given", "selector::selector");
  }
}

AtomPointerVector selector::select(string selStr) {
  expressionTree* tree = buildExpressionTree(selStr);
  AtomPointerVector sel;
  select(tree, sel);
  delete tree;
  return sel;
}

vector<Residue*> selector::selectRes(string selStr) {
  AtomPointerVector sel = select(selStr);
  map<Residue*, bool> selResMap;
  for (int i = 0; i < sel.size(); i++) selResMap[sel[i]->getParent()] = true;
  vector<Residue*> selRes(selResMap.size());
  int i = 0;
  for (map<Residue*, bool>::iterator it = selResMap.begin(); it != selResMap.end(); ++it, i++) {
    selRes[i] = it->first;
  }
  return selRes;
}

void selector::select(expressionTree* tree, AtomPointerVector& sel) {
  if (tree->numChildren() == 0) {
    // this is a terminal node, so just do the selection
    for (int i = 0; i < atoms.size(); i++) {
      switch(tree->getProperty()) {
        case (expressionTree::selProperty::RESID):
          if (residues[i]->getNum() == tree->getNum()) sel.push_back(atoms[i]);
          break;
        case (expressionTree::selProperty::RESNAME):
          if (residues[i]->getName() == tree->getString()) sel.push_back(atoms[i]);
          break;
        case (expressionTree::selProperty::CHAIN):
          if (chains[i]->getID() == tree->getString()) sel.push_back(atoms[i]);
          break;
        case (expressionTree::selProperty::SEGID):
          if (chains[i]->getSegID() == tree->getString()) sel.push_back(atoms[i]);
          break;
        case (expressionTree::selProperty::NAME):
          if (atoms[i]->getName() == tree->getString()) sel.push_back(atoms[i]);
          break;
        default:
          MstUtils::error("uknown selectable property " + MstUtils::toString(tree->getProperty()), "selector::select");
      }
    }
  } else {
    AtomPointerVector selA, selB;
    switch(tree->getLogicalOperator()) {
      case (expressionTree::logicalOp::AND):
        if (tree->numChildren() != 2)
          MstUtils::error("poorly parsed expressoin: expected two operands for AND", "selector::select(expressionTree* )");
        select(tree->getChild(0), selA);
        select(tree->getChild(1), selB);
        sel = intersect(selA, selB);
        break;
      case (expressionTree::logicalOp::OR):
        if (tree->numChildren() != 2)
          MstUtils::error("poorly parsed expressoin: expected two operands for OR", "selector::select(expressionTree* )");
        select(tree->getChild(0), selA);
        select(tree->getChild(1), selB);
        sel = combine(selA, selB);
        break;
      case (expressionTree::logicalOp::NOT):
        if (tree->numChildren() != 1)
          MstUtils::error("poorly parsed expressoin: expected one operand for NOT", "selector::select(expressionTree* )");
        select(tree->getChild(0), selA);
        sel = invert(selA);
        break;
      case (expressionTree::logicalOp::BYRES):
        if (tree->numChildren() != 1)
          MstUtils::error("poorly parsed expressoin: expected one operand for BYRES", "selector::select(expressionTree* )");
        select(tree->getChild(0), selA);
        sel = byRes(selA);
        break;
      case (expressionTree::logicalOp::BYCHAIN):
        if (tree->numChildren() != 1)
          MstUtils::error("poorly parsed expressoin: expected one operand for BYCHAIN", "selector::select(expressionTree* )");
        select(tree->getChild(0), selA);
        sel = byChain(selA);
        break;
      case (expressionTree::logicalOp::IS):
        if (tree->numChildren() != 1)
          MstUtils::error("poorly parsed expressoin: expected one operand for IS", "selector::select(expressionTree* )");
        select(tree->getChild(0), sel);
        break;
      case (expressionTree::logicalOp::AROUND):
        if (tree->numChildren() != 1)
          MstUtils::error("poorly parsed expressoin: expected one selection operand for AROUND", "selector::select(expressionTree* )");
        select(tree->getChild(0), selA);
        sel = around(selA, tree->getVal());
        break;
      default:
        MstUtils::error("uknown selectable property " + MstUtils::toString(tree->getProperty()), "selector::select");
    }
  }
}

expressionTree* selector::buildExpressionTree(string selStr) {
  expressionTree* tree = new expressionTree();
  string token = getNextSelectionToken(selStr); // either something in () or the next space-separated word

  if (MstUtils::stringsEqual(token, "not")) {
    tree->setLogicalOperator(expressionTree::logicalOp::NOT);
    tree->addChild(buildExpressionTree(selStr));
    return tree;
  } else if (MstUtils::stringsEqual(token, "byres")) {
    tree->setLogicalOperator(expressionTree::logicalOp::BYRES);
    tree->addChild(buildExpressionTree(selStr));
    return tree;
  } else if (MstUtils::stringsEqual(token, "bychain")) {
    tree->setLogicalOperator(expressionTree::logicalOp::BYCHAIN);
    tree->addChild(buildExpressionTree(selStr));
    return tree;
  } else if (MstUtils::stringsEqual(token, "resid")) {
    string str = getNextSelectionToken(selStr);
    if (!MstUtils::isInt(str)) MstUtils::error("bad selection, expected number when saw " + str, "selector::buildExpressionTree(string)");
    tree->setLogicalOperator(expressionTree::logicalOp::IS);
    tree->setProperty(expressionTree::selProperty::RESID);
    tree->setNum(MstUtils::toInt(str));
  } else if (MstUtils::stringsEqual(token, "resname")) {
    string str = getNextSelectionToken(selStr);
    tree->setLogicalOperator(expressionTree::logicalOp::IS);
    tree->setProperty(expressionTree::selProperty::RESNAME);
    tree->setString(str);
  } else if (MstUtils::stringsEqual(token, "chain")) {
    string str = getNextSelectionToken(selStr);
    tree->setLogicalOperator(expressionTree::logicalOp::IS);
    tree->setProperty(expressionTree::selProperty::CHAIN);
    tree->setString(str);
  } else if (MstUtils::stringsEqual(token, "segid")) {
    string str = getNextSelectionToken(selStr);
    tree->setLogicalOperator(expressionTree::logicalOp::IS);
    tree->setProperty(expressionTree::selProperty::SEGID);
    tree->setString(str);
  } else if (MstUtils::stringsEqual(token, "name")) {
    string str = getNextSelectionToken(selStr);
    tree->setLogicalOperator(expressionTree::logicalOp::IS);
    tree->setProperty(expressionTree::selProperty::NAME);
    tree->setString(str);
  } else if (token[0] == '(') {
    tree->setLogicalOperator(expressionTree::logicalOp::IS);
    tree->addChild(buildExpressionTree(token.substr(1, token.size()-2)));
  } else {
    MstUtils::error("do not know how to parse token '" + token + "'", "selector::buildExpressionTree");
  }

  string connector = getNextSelectionToken(selStr);
  if (connector.empty()) {
    return tree;
  }
  expressionTree* root = new expressionTree();
  root->addChild(tree);
  if (MstUtils::stringsEqual(connector, "and")) {
    root->setLogicalOperator(expressionTree::logicalOp::AND);
  } else if (MstUtils::stringsEqual(connector, "or")) {
    root->setLogicalOperator(expressionTree::logicalOp::OR);
  } else if (MstUtils::stringsEqual(connector, "around")) {
    root->setLogicalOperator(expressionTree::logicalOp::AROUND);
    root->setVal(MstUtils::toReal(getNextSelectionToken(selStr)));
    string extra = getNextSelectionToken(selStr);
    MstUtils::assert(extra.empty(), "poorly formed selection expression: extra stuff after an 'around' operator; use parentheses as necessary.");
    return root;
  } else {
    MstUtils::error("bad selection, unrecognized connector keyword '" + connector + "' after token '" + token + "'", "selector::buildExpressionTree(string)");
  }

  root->addChild(buildExpressionTree(selStr));
  return root;
}

string selector::getNextSelectionToken(string& selStr) {
  selStr = MstUtils::trim(selStr);
  if (selStr.empty()) { return ""; }

  if (selStr[0] == '(') { // if parenthetical, find matching paren
    // find matching closing paren
    int n = 1; int i;
    for (i = 1; i < selStr.size(); i++) {
      if (selStr[i] == '(') n++;
      if (selStr[i] == ')') n--;
      if (n == 0) break;
    }
    if (n != 0) MstUtils::error("ill-formed selection expression '" + selStr + "'", "Structure::getNextSelectionToken");
    string token = selStr.substr(0, i+1); // keep trailing ()
    selStr = selStr.substr(i+1);
    return token;
  }

  // otherwise, just get the next space-delimited token
  string token = MstUtils::nextToken(selStr);
  return token;
}

AtomPointerVector selector::invert(AtomPointerVector& selAtoms) {
  map<Atom*, bool> in;
  AtomPointerVector notIn;
  for (int i = 0; i < selAtoms.size(); i++) in[selAtoms[i]] = true;
  for (int i = 0; i < atoms.size(); i++) {
    if (in.find(atoms[i]) == in.end()) notIn.push_back(atoms[i]);
  }
  return notIn;
}

AtomPointerVector selector::intersect(AtomPointerVector& selA, AtomPointerVector& selB) {
  map<Atom*, bool> inA;
  for (int i = 0; i < selA.size(); i++) inA[selA[i]] = true;

  AtomPointerVector common;
  for (int i = 0; i < selB.size(); i++) {
    if (inA.find(selB[i]) != inA.end()) common.push_back(selB[i]);
  }
  return common;
}

AtomPointerVector selector::combine(AtomPointerVector& selA, AtomPointerVector& selB) {
  map<Atom*, bool> inBothMap;
  for (int i = 0; i < selA.size(); i++) inBothMap[selA[i]] = true;
  for (int i = 0; i < selB.size(); i++) inBothMap[selB[i]] = true;
  AtomPointerVector inBoth(inBothMap.size(), NULL);
  int i = 0;
  for (map<Atom*, bool>::iterator it = inBothMap.begin(); it != inBothMap.end(); ++it, i++)
    inBoth[i] = it->first;
  return inBoth;
}

AtomPointerVector selector::around(AtomPointerVector& selAtoms, mstreal dcut) {
  AtomPointerVector within;
  mstreal dcut2 = dcut*dcut;
  for (int i = 0; i < atoms.size(); i++) {
    Atom* a = atoms[i];
    for (int j = 0; j < selAtoms.size(); j++) {
      if (a->distance2(selAtoms[j]) <= dcut2) {
        within.push_back(a);
        break;
      }
    }
  }
  return within;
}

AtomPointerVector selector::byRes(AtomPointerVector& selAtoms) {
  set<Residue*> added;
  AtomPointerVector expanded;
  for (int i = 0; i < selAtoms.size(); i++) {
    Residue* res = selAtoms[i]->getParent();
    if (res == NULL) MstUtils::error("some atoms in selection are parentless", "selector::byRes");
    if (added.find(res) != added.end()) continue;
    added.insert(res); // this will accumulate residues in the "natural" order (according to Residue::operator<)
  }
  for (auto it = added.begin(); it != added.end(); ++it) {
    Residue* res = *it;
    for (int j = 0; j < res->atomSize(); j++) expanded.push_back(&(res->getAtom(j)));
  }

  return expanded;
}

AtomPointerVector selector::byChain(AtomPointerVector& selAtoms) {
  map<Chain*, bool> added;
  AtomPointerVector expanded;
  for (int i = 0; i < selAtoms.size(); i++) {
    Chain* chain = selAtoms[i]->getChain();
    if (chain == NULL) MstUtils::error("some atoms in selection are parentless", "selector::byChain");
    if (added.find(chain) != added.end()) continue;
    for (int j = 0; j < chain->residueSize(); j++) {
      Residue& res = chain->getResidue(j);
      for (int k = 0; k < res.atomSize(); k++) expanded.push_back(&(res[k]));
    }
    added[chain] = true;
  }
  return expanded;
}

/* --------- RMSDCalculator --------- */

vector<mstreal> RMSDCalculator::lastTranslation() {
    vector<mstreal> trans(3, 0.0);
    for (int i = 0; i < 3; i++) trans[i] = t[i];
    return trans;
}

vector<vector<mstreal> > RMSDCalculator::lastRotation() {
    vector<vector<mstreal> > rot(3);
    for (int i = 0; i < 3; i++) {
        rot[i].resize(3, 0);
        for (int j = 0; j < 3; j++ ) {
            rot[i][j] = u[i][j];
        }
    }
    return rot;
}

mstreal RMSDCalculator::bestRMSD(const vector<Atom*> &_align, const vector<Atom*> &_ref, bool setTransRot, bool* _suc) {
    _res = 0.0;
    if (Kabsch(_align, _ref, setTransRot)) { if (_suc != NULL) *_suc = true; }
    else { if (_suc != NULL) *_suc = false; }
    return sqrt(_res/_n);
}

mstreal RMSDCalculator::bestResidual(const vector<Atom*> &_align, const vector<Atom*> &_ref, bool setTransRot, bool* _suc) {
    _res = 0.0;
    if (Kabsch(_align, _ref, setTransRot)) { if (_suc != NULL) *_suc = true; }
    else { if (_suc != NULL) *_suc = false; }
    return _res;
}

bool RMSDCalculator::align(const vector<Atom*> &_align, const vector<Atom*> &_ref, vector<Atom*>& _moveable) {
    bool suc = Kabsch(_align, _ref, 1);

    if (suc) {
        mstreal x[3],x1[3];
        for(int k=0; k<_moveable.size(); k++) {
            x[0]=_moveable[k]->getX();
            x[1]=_moveable[k]->getY();
            x[2]=_moveable[k]->getZ();
            x1[0] = t[0]+u[0][0]*x[0]+u[0][1]*x[1]+u[0][2]*x[2];
            x1[1] = t[1]+u[1][0]*x[0]+u[1][1]*x[1]+u[1][2]*x[2];
            x1[2] = t[2]+u[2][0]*x[0]+u[2][1]*x[1]+u[2][2]*x[2];
            _moveable[k]->setCoor(x1[0], x1[1], x1[2]);
        }
    }
    return suc;
}

bool RMSDCalculator::align(const vector<Atom*> &_align, const vector<Atom*> &_ref, Structure& _moveable) {
  vector<Atom*> atoms = _moveable.getAtoms();
  return align(_align, _ref, atoms);
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
bool RMSDCalculator::Kabsch(const vector<Atom*> &_align, const vector<Atom*> &_ref, int mode) {
    int i, j, m, m1, l, k;
    mstreal e0, rms1, d, h, g;
    mstreal cth, sth, sqrth, p, det, sigma;
    mstreal xc[3], yc[3];
    mstreal a[3][3], b[3][3], r[3][3], e[3], rr[6], ss[6];
    mstreal sqrt3=1.73205080756888, tol=0.01;
    int ip[]={0, 1, 3, 1, 2, 4, 3, 4, 5};
    int ip2312[]={1, 2, 0, 1};

    int a_failed=0, b_failed=0;
    mstreal epsilon=0.00000001;

    int n=_ref.size();
    if(n != _align.size()) {
        cout << "Two proteins have different length!" << endl;
        return false;
    }

    //initializtation
    _res=0;
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
        xc[i] = xc[i]/n;
        yc[i] = yc[i]/n;
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

    mstreal spur=(rr[0]+rr[2]+rr[5]) / 3.0;
    mstreal cof = (((((rr[2]*rr[5] - rr[4]*rr[4]) + rr[0]*rr[5]) \
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

//    _rmsd=sqrt(rms1/(mstreal)n);
    _res = rms1;
    _n = n;

    return true;
}

mstreal RMSDCalculator::rmsd(const vector<Atom*>& A, const vector<Atom*>& B) {
  if (A.size() != B.size())
    MstUtils::error("atom vectors of different length (" + MstUtils::toString(A.size()) + " and " + MstUtils::toString(B.size()) + ")", "RMSDCalculator::rmsd(vector<Atom*>&, vector<Atom*>&)");

  mstreal ret = 0;
  for (int i = 0; i < A.size(); i++) {
    ret += A[i]->distance2(B[i]);
  }
  return sqrt(ret/A.size());
}

mstreal RMSDCalculator::rmsd(const Structure& A, const Structure& B) {
  vector<Atom*> atomsA = A.getAtoms();
  vector<Atom*> atomsB = B.getAtoms();
  return rmsd(atomsA, atomsB);
}

mstreal RMSDCalculator::rmsdCutoff(const vector<int>& L, mstreal rmsdMax, mstreal L0) {
  mstreal a = (mstreal) exp(-1./L0);
  int N = 0, n;
  mstreal c = 0;

  // disjoint segments are counted as independent, so their correlation
  // with respect to each other is zero
  for (int i = 0; i < L.size(); i++) {
    N += L[i];
    n = L[i];
    c = c + (a/(1-a))*(n-1) - pow((a/(1-a)), 2)*(1 - pow(a, n-1));
  }
  double df = N*(1 - (2.0/(N*(N-1)))*c);

  return rmsdMax/sqrt(N/df);
}

mstreal RMSDCalculator::rmsdCutoff(const Structure& S, mstreal rmsdMax, mstreal L0) {
  vector<int> L(S.chainSize());
  for (int i = 0; i < S.chainSize(); i++) L[i] = S[i].residueSize();
  return RMSDCalculator::rmsdCutoff(L, rmsdMax, L0);
}

/* --------- ProximitySearch --------- */

ProximitySearch::ProximitySearch(mstreal _xlo, mstreal _ylo, mstreal _zlo, mstreal _xhi, mstreal _yhi, mstreal _zhi, int _N) {
  xlo = _xlo; ylo = _ylo; zlo = _zlo;
  xhi = _xhi; yhi = _yhi; zhi = _zhi;
  reinitBuckets(_N);
  setBinWidths();
}

ProximitySearch::ProximitySearch(const AtomPointerVector& _atoms, int _N, bool _addAtoms, vector<int>* tags, mstreal pad) {
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

ProximitySearch::ProximitySearch(const AtomPointerVector& _atoms, mstreal _characteristicDistance, bool _addAtoms, vector<int>* tags, mstreal pad) {
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

ProximitySearch::ProximitySearch(const ProximitySearch& ps) {
  *this = ps;
}

void ProximitySearch::setBinWidths() {
  xbw = (xhi - xlo)/(N - 1);
  ybw = (yhi - ylo)/(N - 1);
  zbw = (zhi - zlo)/(N - 1);
}

ProximitySearch::~ProximitySearch() {
  for (int i = 0; i < pointList.size(); i++) delete(pointList[i]);
}

ProximitySearch& ProximitySearch::operator=(const ProximitySearch& ps) {
  N = ps.N;
  xlo = ps.xlo; ylo = ps.ylo; zlo = ps.zlo;
  xhi = ps.xhi; yhi = ps.yhi; zhi = ps.zhi;
  xbw = ps.xbw; ybw = ps.ybw; zbw = ps.zbw;
  buckets = ps.buckets;
  pointTags = pointTags;
  for (int i = 0; i < pointList.size(); i++) delete(pointList[i]);
  pointList.resize(ps.pointList.size());
  for (int i = 0; i < ps.pointList.size(); i++) {
    pointList[i] = new CartesianPoint(*(ps.pointList[i]));
  }
  return *this;
}

void ProximitySearch::calculateExtent(const Structure& S, mstreal& _xlo, mstreal& _ylo, mstreal& _zlo, mstreal& _xhi, mstreal& _yhi, mstreal& _zhi) {
  AtomPointerVector atoms = S.getAtoms();
  calculateExtent(atoms, _xlo, _ylo, _zlo, _xhi, _yhi, _zhi);
}

void ProximitySearch::calculateExtent(const AtomPointerVector& _atoms, mstreal& _xlo, mstreal& _ylo, mstreal& _zlo, mstreal& _xhi, mstreal& _yhi, mstreal& _zhi) {
  if (_atoms.size() == 0) { cout << "Error in ProximitySearch::calculateExtent() -- empty atom vector passed!\n"; exit(-1); }
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
    }
  }
  pointList.resize(0);
  pointTags.resize(0);
}

void ProximitySearch::addPoint(const CartesianPoint& _p, int tag) {
  int i, j, k;
  pointBucket(_p[0], _p[1], _p[2], &i, &j, &k);
  if ((i < 0) || (j < 0) || (k < 0) || (i > N-1) || (j > N-1) || (k > N-1)) { cout << "Error: point " << _p << " out of range for ProximitySearch object!\n"; exit(-1); }
  CartesianPoint* p = new CartesianPoint(_p);
  vector<int>& bijk = buckets[i][j][k];
  if (bijk.empty()) fullBuckets.push_back(&bijk);
  bijk.push_back(pointList.size());
  pointList.push_back(p);
  pointTags.push_back(tag);
}

void ProximitySearch::addPoint(mstreal xc, mstreal yc, mstreal zc, int tag) {
  int i, j, k;
  pointBucket(xc, yc, zc, &i, &j, &k);
  if ((i < 0) || (j < 0) || (k < 0) || (i > N-1) || (j > N-1) || (k > N-1)) { cout << "Error: point " << xc << " " << yc << " " << zc << " out of range for ProximitySearch object!\n"; exit(-1); }
  CartesianPoint* p = new CartesianPoint(xc, yc, zc);
  vector<int>& bijk = buckets[i][j][k];
  if (bijk.empty()) fullBuckets.push_back(&bijk);
  bijk.push_back(pointList.size());
  pointList.push_back(p);
  pointTags.push_back(tag);
}

void ProximitySearch::addAtoms(AtomPointerVector& apv, vector<int>* tags) {
  if ((tags != NULL) && (apv.size() != tags->size())) MstUtils::error("different number of atoms and tags specified!", "ProximitySearch::addAtoms");
  for (int i = 0; i < apv.size(); i++) {
    addPoint(apv[i], (tags == NULL) ? i : (*tags)[i]);
  }
}

void ProximitySearch::dropAllPoints() {
  for (int i = 0; i < pointList.size(); i++) delete(pointList[i]);
  pointList.resize(0);
  pointTags.resize(0);
  for (int i = 0; i < fullBuckets.size(); i++) fullBuckets[i]->resize(0);
  fullBuckets.resize(0);
}

bool ProximitySearch::isPointWithinGrid(CartesianPoint _p) {
  int i, j, k;
  pointBucket(&_p, &i, &j, &k);
  if ((i < 0) || (j < 0) || (k < 0) || (i > N-1) || (j > N-1) || (k > N-1)) return false;
  return true;
}

void ProximitySearch::pointBucket(mstreal px, mstreal py, mstreal pz, int* i, int* j, int* k) {
  *i = (int) floor((px - xlo)/xbw + 0.5);
  *j = (int) floor((py - ylo)/ybw + 0.5);
  *k = (int) floor((pz - zlo)/zbw + 0.5);
}

void ProximitySearch::limitIndex(int *ind) {
  if (*ind < 0) *ind = 0;
  if (*ind > N-1) *ind = N-1;
}

bool ProximitySearch::pointsWithin(const CartesianPoint& c, mstreal dmin, mstreal dmax, vector<int>* list, bool byTag) {
  mstreal cx = c.getX(); mstreal cy = c.getY(); mstreal cz = c.getZ();
  // first check if the point is outside of the bounding box of the point cloud by a sufficient amount
  if ((cx < xlo - dmax) || (cy < ylo - dmax) || (cz < zlo - dmax) || (cx > xhi + dmax) || (cy > yhi + dmax) || (cz > zhi + dmax)) return false;

  mstreal d2, dmin2, dmax2;
  int ci, cj, ck;
  int iOutLo, jOutLo, kOutLo, iOutHi, jOutHi, kOutHi; // external box (no point in looking beyond it, points there are too far)
  int iInLo, jInLo, kInLo, iInHi, jInHi, kInHi;       // internal box (no point in looking within it, points there are too close)
  pointBucket(cx, cy, cz, &ci, &cj, &ck);
  pointBucket(cx - dmax, cy - dmax, cz - dmax, &iOutLo, &jOutLo, &kOutLo);
  pointBucket(cx + dmax, cy + dmax, cz + dmax, &iOutHi, &jOutHi, &kOutHi);
  if (dmin > 0) {
    mstreal sr3 = sqrt(3);
    pointBucket(cx - dmin/sr3, cy - dmin/sr3, cz - dmin/sr3, &iInLo, &jInLo, &kInLo);
    pointBucket(cx + dmin/sr3, cy + dmin/sr3, cz + dmin/sr3, &iInHi, &jInHi, &kInHi);
    // NOTE: I used to think this adjustment was necessary, but upon more thought, I don't think so
    // // need to trim the internal box to make sure it is fully contained within the sphere of radius dmin from the central point
    // if (iInLo != ci) iInLo++;
    // if (jInLo != cj) jInLo++;
    // if (kInLo != ck) kInLo++;
    // if (iInHi != ci) iInHi--;
    // if (jInHi != cj) jInHi--;
    // if (kInHi != ck) kInHi--;
  } else {
    iInLo = iInHi = ci;
    jInLo = jInHi = cj;
    kInLo = kInHi = ck;
  }
  limitIndex(&iInLo); limitIndex(&iInHi); limitIndex(&jInLo); limitIndex(&jInHi); limitIndex(&kInLo); limitIndex(&kInHi);
  limitIndex(&iOutLo); limitIndex(&iOutHi); limitIndex(&jOutLo); limitIndex(&jOutHi); limitIndex(&kOutLo); limitIndex(&kOutHi);

  // search only within the boxes where points of interest can be, in principle
  if (list != NULL) list->clear();
  bool found = false, insi, ins;
  bool yesno = (list == NULL);
  dmin2 = dmin*dmin; dmax2 = dmax*dmax;
  int i, j, k, ii, pi;
  for (i = iOutLo; i <= iOutHi; i++) {
    insi = (i >= iInLo) && (i <= iInHi);
    vector<vector<vector<int> > >& Bi = buckets[i];
    for (j = jOutLo; j <= jOutHi; j++) {
      ins = insi && (j >= jInLo) && (j <= jInHi);
      vector<vector<int> >& Bij = Bi[j];
      for (k = kOutLo; k <= kOutHi; k++) {
        vector<int>& Bijk = Bij[k];
        // check all points in bucket i, j, k
        for (ii = 0; ii < Bijk.size(); ii++) {
          pi = Bijk[ii];
          d2 = c.distance2nc(*(pointList[pi]));
          if ((d2 >= dmin2) && (d2 <= dmax2)) {
            if (yesno) return true;
            list->push_back(byTag ? pointTags[Bijk[ii]] : Bijk[ii]);
            found = true;
          }
        }
        // skip the range from kInLo to kInHi (too close)
        if (ins && (k == kInLo) && (kInLo != kInHi)) k = kInHi - 1;
      }
    }
  }
  return found;
}

vector<int> ProximitySearch::getPointsWithin(const CartesianPoint& c, mstreal dmin, mstreal dmax, bool byTag) {
  vector<int> closeOnes;
  pointsWithin(c, dmin, dmax, &closeOnes, byTag);
  return closeOnes;
}

bool ProximitySearch::overlaps(ProximitySearch& other, mstreal pad) {
  if ((other.xlo > xhi + pad) || (xlo > other.xhi + pad) || (other.ylo > yhi + pad) || (ylo > other.yhi + pad) || (other.zlo > zhi + pad) || (zlo > other.zhi + pad)) {
    return false;
  }
  return true;
}

/* --------- Clusterer --------- */
vector<vector<int> > Clusterer::greedyCluster(const vector<vector<Atom*> >& units, mstreal rmsdCut, int Nmax) {
  coputedRMSDs.clear();
  vector<vector<int> > clusters;
  set<int> remIndices;
  for (int i = 0; i < units.size(); i++) remIndices.insert(i);
  if (remIndices.size() >= Nmax) return Clusterer::greedyClusterBruteForce(units, remIndices, rmsdCut);

  while (remIndices.size() > Nmax) {
    // sub-sample Nmax elements
    set<int> subSample = Clusterer::randomSubsample(remIndices, Nmax);

    // get the top cluster from these and use its centroid
    vector<int> topClustSub = Clusterer::greedyClusterBruteForce(units, subSample, rmsdCut, 1)[0];
    vector<int> topClust = Clusterer::elementsWithin(units, remIndices, topClustSub[0], rmsdCut);

    // now try to grow the cluster by improving the centroid, if there are compute cycles left
    int Ntry = int(round(Nmax * Nmax * 1.0 / remIndices.size()));
    for (int i = 1; i <= MstUtils::min(Ntry, (int) topClust.size()-1); i++) {
      vector<int> newTopClust = Clusterer::elementsWithin(units, remIndices, topClust[i], rmsdCut);
      if (newTopClust.size() > topClust.size()) topClust = newTopClust;
    }

    // keep whatever cluster end up with, exclude its elements
    clusters.push_back(topClust);
    for (int i = 0; i < topClust.size(); i++) remIndices.erase(topClust[i]);
  }

  // brute force through the rest
  vector<vector<int> > remClusters = greedyClusterBruteForce(units, remIndices, rmsdCut);
  clusters.insert(clusters.end(), remClusters.begin(), remClusters.end());
  return clusters;
}

vector<vector<int> > Clusterer::greedyClusterBruteForce(const vector<vector<Atom*> >& units, set<int> remIndices, mstreal rmsdCut, int nClusts) {
  vector<vector<int> > clusters;
  while (remIndices.size() != 0) {
    // pick the best current centroid
    vector<int> bestClust;
    for (auto it = remIndices.begin(); it != remIndices.end(); ++it) {
      vector<int> clust = elementsWithin(units, remIndices, *it, rmsdCut);
      if (clust.size() > bestClust.size()) bestClust = clust;
    }

    // add its corresponding cluster
    clusters.push_back(bestClust);
    for (int i = 0; i < bestClust.size(); i++) remIndices.erase(bestClust[i]);

    // where we asked for at most the top some number of clusters?
    if ((nClusts > 0) && (clusters.size() >= nClusts)) break;
  }
  return clusters;
}

vector<int> Clusterer::elementsWithin(const vector<vector<Atom*> >& units, set<int>& remIndices, int from, mstreal rmsdCut) {
  vector<int> neigh; vector<mstreal> rmsds;
  const vector<Atom*>& fromUnit = units[from];
  RMSDCalculator rCalc;
  for (auto it = remIndices.begin(); it != remIndices.end(); ++it) {
    mstreal r = rCalc.bestRMSD(fromUnit, units[*it]);
    if (r <= rmsdCut) {
      neigh.push_back(*it);
      rmsds.push_back(r);
    }
  }

  // sort by ascending RMSD
  vector<int> si = MstUtils::sortIndices(rmsds);
  vector<int> orderedNeigh = neigh;
  for (int i = 0; i < si.size(); i++) {
    orderedNeigh[i] = neigh[si[i]];
  }

  return orderedNeigh;
}

set<int> Clusterer::randomSubsample(set<int>& indices, int N) {
  if (N > indices.size())
    MstUtils::error("asked for a subsample of " + MstUtils::toString(N) + " elements from an array of " + MstUtils::toString(indices.size()) + " elements", "Clusterer::randomSubsample");

  vector<int> inds(indices.size(), 0);
  int k = 0;
  for (auto it = indices.begin(); it != indices.end(); ++it) {
    inds[k] = *it; k++;
  }
  MstUtils::shuffle(inds);

  set<int> sub;
  for (int i = 0; i < N; i++) sub.insert(inds[i]);
  return sub;
}

vector<vector<int> > Clusterer::kmeans(const vector<CartesianPoint>& points, int k, int Ntrials, int Niter) {
  if (k > points.size()) MstUtils::error("asked for " + MstUtils::toString(k) + " means, but there are only " + MstUtils::toString(points.size()) + " points in the cloud!", "Clusterer::kmeans");
  if (k == 0) return vector<vector<int> >(); // could also error
  mstreal dist, bestDist, err, bestWCSS, wcss;
  vector<vector<int> > bestClusts(k);
  int best;

  for (int t = 0; t < Ntrials; t++) {
    vector<vector<int> > clusts(k);

    // randomly permute and take the first k as initial means
    vector<CartesianPoint> shuffledPoints = points;
    MstUtils::shuffle(shuffledPoints);
    vector<CartesianPoint> means(k);
    for (int i = 0; i < k; i++) means[i] = shuffledPoints[i];

    // now cycle to update according to the Lloyd's algorithm
    bool ok = true;
    for (int c = 0; c < Niter; c++) {
      for (int j = 0; j < k; j++) clusts[j].resize(0);
      wcss = 0; // within-cluster sum of squares
      // assign each point to the cluster with closest mean
      for (int i = 0; i < points.size(); i++) {
        const CartesianPoint& p = points[i];
        best = 0; bestDist = p.distance(means[0]);
        for (int j = 1; j < means.size(); j++) {
          dist = p.distance(means[j]);
          if (dist < bestDist) {
            bestDist = dist; best = j;
          }
        }
        clusts[best].push_back(i);
        wcss += bestDist*bestDist;
      }

      // recompute means
      err = 0;
      for (int i = 0; i < clusts.size(); i++) {
        vector<int>& clust = clusts[i];
        if (clust.empty()) { ok = false; break; }
        CartesianPoint o = means[i];
        CartesianPoint& m = means[i]; m = CartesianPoint(m.size(), 0.0);
        for (int j = 0; j < clust.size(); j++) m += points[clust[j]];
        m /= clust.size();
        err += (m - o).norm2();
      }
      if (!ok) break;
      if (err < 10E-8) break;
    }
    if (!ok) { // if one of the clusters ended up empty, try again
      t--;
    } else if ((t == 0) || (wcss < bestWCSS)) {
      bestWCSS = wcss;
      bestClusts = clusts;
    }
  }
  return bestClusts;
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
    // if eof set upon trying to read the line, and the line is empty, then it
    // was not really a valid line (but just the last eof that was not yet read)
    if (inp.eof() && line.empty()) break;
    lines.push_back(line);
  }
}

string MstUtils::nextToken(string& str, string delimiters, bool skipTrailingDelims) {
  string ret; int i;
  if (!delimiters.empty()) {
    if (skipTrailingDelims) str = trim(str, delimiters);
    i = str.find_first_of(delimiters);
    if (i == string::npos) {
      ret = str;
      str = "";
      return ret;
    }
  } else {
    i = min(1, (int) str.length()); // interpret an empty list of delimiters in the same way as perl's split("", $string)
  }
  ret = str.substr(0, i);
  str = str.substr(i);
  return ret;
}

vector<string> MstUtils::split(const string& str, string delimiters, bool skipTrailingDelims) {
  string strCopy = str;
  vector<string> tokens;
  while (!strCopy.empty()) tokens.push_back(nextToken(strCopy, delimiters, skipTrailingDelims));
  return tokens;
}

string MstUtils::join(const string& delim, const vector<string>& words) {
  string joined;
  if (words.empty()) return joined;
  joined = words[0];
  for (int i = 1; i < words.size(); i++) joined += delim + words[i];
  return joined;
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

string MstUtils::uc(const string& str){
  string ret = str;
  for (int i = 0; i < ret.length(); i++) {
    ret[i] = toupper(ret[i]);
  }
  return ret;
}

string MstUtils::lc(const string& str){
  string ret = str;
  for (int i = 0; i < ret.length(); i++) {
    ret[i] = tolower(ret[i]);
  }
  return ret;
}

bool MstUtils::stringsEqual(const string& A, const string& B, bool caseInsensitive) {
  if (caseInsensitive) return (strcasecmp(A.c_str(), B.c_str()) == 0);
  return (strcmp(A.c_str(), B.c_str()) == 0);
}

string MstUtils::trim(string str, string delimiters) {
  int i = str.find_first_not_of(delimiters);
  if (i == string::npos) return "";
  int j = str.find_last_not_of(delimiters);
  return str.substr(i, j - i + 1);
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

bool MstUtils::isInt(string num) {
  int ret;
  return (sscanf(num.c_str(), "%d", &ret) != 0);
}

bool MstUtils::isReal(string num) {
  double ret;
  return (sscanf(num.c_str(), "%lf", &ret) != 0);
}

MST::mstreal MstUtils::toReal(string num, bool strict) {
  double ret;
  if ((sscanf(num.c_str(), "%lf", &ret) != 1) && strict) MstUtils::error("failed to convert '" + num + "' to mstreal", "MstUtils::toReal");
  return (mstreal) ret;
}

MST::mstreal MstUtils::mod(MST::mstreal num, MST::mstreal den) {
  return num - ((MST::mstreal) floor((double) num / (double) den)) * den;
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

string MstUtils::getDate() {
  time_t rawtime;
  time (&rawtime);
  return string(ctime(&rawtime));
}
