#include "msttypes.h"
mt19937 MstUtils::mt;

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

Structure::Structure(istream& is, string options) {
  name = "";
  numResidues = numAtoms = 0;
  readPDB(is, options);
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

Structure Structure::combine(const Structure& atomsStruct, const Structure& topoStruct, bool renameResidues) {
  Structure combo = topoStruct;
  vector<Residue*> atomResidues = atomsStruct.getResidues();
  vector<Residue*> comboResidues = combo.getResidues();
  if (atomResidues.size() != comboResidues.size()) MstUtils::error("numbers of residues in structures disagrees", "Structure::combine");
  for (int i = 0; i < comboResidues.size(); i++) {
    comboResidues[i]->copyAtoms(*(atomResidues[i]));
    if (renameResidues) comboResidues[i]->setName(atomResidues[i]->getName());
  }

  return combo;
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

void Structure::readPDB(const string& pdbFile, string options) {
  name = pdbFile;
  fstream ifh; MstUtils::openFile(ifh, pdbFile, fstream::in, "Structure::readPDB");
  readPDB(ifh, options);
  ifh.close();
}

void Structure::readPDB(istream& is, string options) {
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
  while (getline(is, line)) {
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
    mstreal B = MstUtils::toReal(line.substr(60, 6), false);
    mstreal occ = MstUtils::toReal(line.substr(54, 6), false);
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
        MstUtils::warn("chain name '" + chainID + "' was repeated in '" + name + "', renaming the chain to '" + chain->getID() + "'", "Structure::readPDB");
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
    if (fixIleCD1 && (resname.compare("ILE") == 0) && (atomname.compare("CD1") == 0)) atomname = "CD";

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
}

void Structure::writePDB(const string& pdbFile, string options) const {
  fstream ofs; MstUtils::openFile(ofs, pdbFile, fstream::out, "Structure::writePDB(string, string)");
  writePDB(ofs, options);
  ofs.close();
}

void Structure::writePDB(ostream& ofs, string options) const {
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

void Structure::writeData(const string& dataFile) const {
  ofstream ofs; MstUtils::openFile(ofs, dataFile, fstream::out | fstream::binary, "Structure::writeData(const string&)");
  writeData(ofs);
  ofs.close();
}

void Structure::writeData(ostream& ofs) const {
  char ter = '\0';
  ofs << getName() << ter;
  MstUtils::writeBin(ofs, chainSize());
  for (int i = 0; i < chainSize(); i++) {
    Chain& chain = getChain(i);
    ofs << chain.getID() << ter << chain.getSegID() << ter;
    MstUtils::writeBin(ofs,chain.residueSize());
    for (int j = 0; j < chain.residueSize(); j++) {
      Residue& res = chain[j];
      ofs << res.getName() << ter << res.getIcode();
      MstUtils::writeBin(ofs, res.getNum());
      MstUtils::writeBin(ofs, res.atomSize());
      for (int k = 0; k < res.atomSize(); k++) {
        Atom& atom = res[k];
        ofs << atom.getName() << ter << atom.getAlt();
        MstUtils::writeBin(ofs, atom.getX());
        MstUtils::writeBin(ofs, atom.getY());
        MstUtils::writeBin(ofs, atom.getZ());
        MstUtils::writeBin(ofs, atom.getOcc());
        MstUtils::writeBin(ofs, atom.getB());
        MstUtils::writeBin(ofs, (char) atom.isHetero());
        MstUtils::writeBin(ofs, atom.getIndex());
        MstUtils::writeBin(ofs, atom.numAlternatives());
        for (int ii = 0; ii < atom.numAlternatives(); ii++) {
          CartesianPoint altCoor = atom.getAltCoor(ii);
          MstUtils::writeBin(ofs, altCoor.getX());
          MstUtils::writeBin(ofs, altCoor.getY());
          MstUtils::writeBin(ofs, altCoor.getZ());
          MstUtils::writeBin(ofs, atom.getAltOcc(ii));
          MstUtils::writeBin(ofs, atom.getAltB(ii));
          ofs << atom.getAltLocID(ii);
        }
      }
    }
  }
}

void Structure::readData(const string& dataFile) {
  ifstream ifs; MstUtils::openFile(ifs, dataFile, fstream::in | fstream::binary, "Structure::readData(const string&)");
  readData(ifs);
  ifs.close();
}

void Structure::readData(istream& ifs) {
  char ter = '\0';
  getline(ifs, name, '\0');
  string resname, atomname;
  char icode, alt;
  int resnum, atominx, nC, nR, nA, numAlt;
  mstreal x, y, z, B, occ;
  char het;

  MstUtils::readBin(ifs, nC);
  for (int i = 0; i < nC; i++) {
    string chainID, segID;
    getline(ifs, chainID, '\0'); getline(ifs, segID, '\0');
    Chain* chain = new Chain(chainID, segID); appendChain(chain);
    MstUtils::readBin(ifs, nR);
    for (int j = 0; j < nR; j++) {
      getline(ifs, resname, '\0');
      MstUtils::readBin(ifs, icode);
      MstUtils::readBin(ifs, resnum);
      Residue* residue = new Residue(resname, resnum, icode);
      chain->appendResidue(residue);
      MstUtils::readBin(ifs, nA);
      for (int k = 0; k < nA; k++) {
        getline(ifs, atomname, '\0');
        MstUtils::readBin(ifs, alt);
        MstUtils::readBin(ifs, x); MstUtils::readBin(ifs, y); MstUtils::readBin(ifs, z);
        MstUtils::readBin(ifs, occ); MstUtils::readBin(ifs, B); MstUtils::readBin(ifs, het);
        MstUtils::readBin(ifs, atominx); MstUtils::readBin(ifs, numAlt);
        Atom* atom = new Atom(atominx, atomname, x, y, z, B, occ, (bool) het, alt);
        residue->appendAtom(atom);
        for (int ii = 0; ii < numAlt; ii++) {
          MstUtils::readBin(ifs, x); MstUtils::readBin(ifs, y); MstUtils::readBin(ifs, z);
          MstUtils::readBin(ifs, occ); MstUtils::readBin(ifs, B); MstUtils::readBin(ifs, alt);
          atom->addAlternative(x, y, z, B, occ, alt);
        }
      }
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
      MstUtils::assertCond(found, "ran out of even multi-character chain names -- your PDB structure really has more than 36,000 chains???", "Structure::appendChain");
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

void Structure::deleteShortChains(int sizeThreshold) {
  vector<Chain*> toDelete;
  for (Chain* C : chains) if (C->residueSize() < sizeThreshold) toDelete.push_back(C);
  for (Chain* C : toDelete) deleteChain(C);
}

Residue& Structure::getResidue(int i) const {
  if ((i < 0) || (i >= residueSize()))
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

vector<Residue*> Structure::getResidues() const {
  vector<Residue*> vec;

  for (int i = 0; i < this->chainSize(); i++) {
    Chain& c = (*this)[i];
    vector<Residue*> chainResidues = c.getResidues();
    vec.insert(vec.end(), chainResidues.begin(), chainResidues.end());
  }

  return vec;
}

void Structure::renumber(int startResNum, int startAtomIndex) {
  int index = startAtomIndex;
  for (int i = 0; i < chains.size(); i++) {
    Chain& chain = (*this)[i];
    for (int j = 0; j < chain.residueSize(); j++) {
      Residue& res = chain[j];
      res.setNum(startResNum+j);
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
  S.setName(getName());
  reassignChainsByConnectivity(S, maxPeptideBond);
  return S;
}

void Structure::reassignChainsByConnectivity(Structure& dest, mstreal maxPeptideBond) {
  if (this->residueSize() == 0) return;
  int k = 0; // index of chain in original structure; used to try to pick the
             // same chain names (in same order) as in the original structure,
             // so that if the connectivity was correct, same chain names remain
  vector<Residue*> residues = this->getResidues();
	Chain* chain = dest.appendChain((k < this->chainSize()) ? (*this)[k].getID() : "A"); k++;
	for (int i = 0; i < residues.size() - 1; i++) {
    chain->appendResidue(new Residue(*residues[i]));
		Atom* atomC = residues[i]->findAtom("C", true);
		Atom* atomN = residues[i + 1]->findAtom("N", true);
    if ((atomC == NULL) || (atomN == NULL)) MstUtils::error("cannot break into disjoint segments as some C or N backbone atoms are missing", "Structure::reassignChainsByConnectivity");
		if (atomC->distance(atomN) > maxPeptideBond) {
      chain = dest.appendChain((k < this->chainSize()) ? (*this)[k].getID() : "A"); k++;
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

int Chain::getIndex() const {
  const Structure* P = getParent();
  if (P == NULL) MstUtils::error("cannot get index of disembodied chain", "Chain::getIndex()");
  int idx = -1;
  for (int i = 0; i < P->chainSize(); i++) {
    if (&(P->getChain(i)) == this) return i;
  }
  MstUtils::error("strange error: Chain does not appear to belong to its parent Structure", "Chain::getIndex()");
  return -1;
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

void Residue::appendAtoms(const vector<Atom*>& A) {
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
  copyAtoms(R.getAtoms(), copyAlt);
}

void Residue::copyAtoms(const vector<Atom*>& _atoms, bool copyAlt) {
  deleteAtoms();
  for (int i = 0; i < _atoms.size(); i++) appendAtom(new Atom(_atoms[i], copyAlt));
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

void Residue::replaceAtoms(const vector<Atom*>& newAtoms, const vector<Atom*>& oldAtoms) {
  vector<int> oldAtomIndices;
  for (int i = 0; i < oldAtoms.size(); i++) {
    for (int k = 0; k < atoms.size(); k++) {
      if (atoms[k] == oldAtoms[i]) {
        oldAtomIndices.push_back(k); break;
      }
    }
  }
  replaceAtoms(newAtoms, &oldAtomIndices);
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

bool Residue::areBonded(const Residue& resN, const Residue& resC, mstreal maxPeptideBond) {
  Atom* atomC = resN.findAtom("C", true);
  Atom* atomN = resC.findAtom("N", true);
  if ((atomC == NULL) || (atomN == NULL)) MstUtils::error("necessary C or N backbone atom(s) missing", "Residue::areBonded");
  if (atomC->distance(atomN) > maxPeptideBond) return false;
  return true;
}

Atom* Residue::findAtom(string _name, bool strict) const {
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

int Residue::getResidueIndexInChain() const {
  if (getParent() == NULL)
    MstUtils::error("cannot find index of a disembodied residue '" + MstUtils::toString(this) + "'", "Residue::getResidueIndexInChain()");
  return getParent()->getResidueIndex(this);
}

/* --------- Atom --------- */
/* --------- Atom --------- */
Atom::atomInfo::atomInfo() {
  parent = NULL;
  het = false;
  name = MstUtils::copyStringC("UNK");
  alternatives = NULL;
  alt = ' ';
  index = 0;
  occ = B = 0;
}

Atom::atomInfo::atomInfo(const atomInfo& other, bool copyAlt) {
  index = other.index;
  name = NULL;
  setName(other.name);
  B = other.B;
  occ = other.occ;
  het = other.het;
  alt = other.alt;
  parent = other.parent;
  if (copyAlt && (other.alternatives != NULL)) {
    alternatives = new vector<altInfo>(other.alternatives->size());
    for (int i = 0; i < other.alternatives->size(); i++) {
      (*alternatives)[i] = (*(other.alternatives))[i];
    }
  } else {
    alternatives = NULL;
  }
}

Atom::atomInfo::atomInfo(int _index, const string& _name, mstreal _B, mstreal _occ, bool _het, char _alt, Residue* _parent) {
  index = _index;
  name = NULL;
  setName(_name);
  B = _B;
  occ = _occ;
  het = _het;
  alt = _alt;
  parent = _parent;
  alternatives = NULL;
}

Atom::atomInfo::~atomInfo() {
  if (name != NULL) delete[] name;
  if (alternatives != NULL) delete alternatives;
}

Atom::Atom() {
  x = y = z = 0;
  info = new atomInfo();
}

Atom::Atom(const Atom& A, bool copyAlt) {
  x = A.x;
  y = A.y;
  z = A.z;
  info = NULL;
  if (A.hasInfo()) info = new atomInfo(*(A.info));
}

Atom::Atom(int _index, const string& _name, mstreal _x, mstreal _y, mstreal _z, mstreal _B, mstreal _occ, bool _het, char _alt, Residue* _parent) {
  x = _x;
  y = _y;
  z = _z;
  info = new atomInfo(_index, _name, _B, _occ, _het, _alt, _parent);
}

Atom::~Atom() {
  if (info != NULL) delete info;
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

const mstreal& Atom::operator[](int i) const {
  switch(i) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      MstUtils::error("invalid coordinate index " + MstUtils::toString(i), "Atom::operator[](int) const");
  }
  return x; // just to silence the warning from some compilres; in reality, this is never reached
}

CartesianPoint Atom::getCoor() const {
  CartesianPoint coor(x, y, z); return coor;
}

void Atom::atomInfo::setName(const char* _name) {
  if (name != NULL) delete[] name;
  name = new char[strlen(_name)+1];
  strcpy(name, _name);
}

CartesianPoint Atom::atomInfo::getAltCoor(int altInd) const {
  if ((alternatives == NULL) || (altInd >= alternatives->size()) || (altInd < 0)) MstUtils::error("alternative index " + MstUtils::toString(altInd) + " out of bounds (" + MstUtils::toString(alternatives->size()) + " alternatives available)", "Atom::getAltCoor");
  altInfo& targ = (*alternatives)[altInd];
  CartesianPoint coor(targ.x, targ.y, targ.z);
  return coor;
}
CartesianPoint Atom::getAltCoor(int altInd) const { return info->getAltCoor(altInd); }

mstreal Atom::atomInfo::getAltB(int altInd) const {
  if ((alternatives == NULL) || (altInd >= alternatives->size()) || (altInd < 0)) MstUtils::error("alternative index " + MstUtils::toString(altInd) + " out of bounds (" + MstUtils::toString(alternatives->size()) + " alternatives available)", "Atom::getAltB");
  return (*alternatives)[altInd].B;
}

mstreal Atom::atomInfo::getAltOcc(int altInd) const {
  if ((alternatives == NULL) || (altInd >= alternatives->size()) || (altInd < 0)) MstUtils::error("alternative index " + MstUtils::toString(altInd) + " out of bounds (" + MstUtils::toString(alternatives->size()) + " alternatives available)", "Atom::getAltOcc");
  return (*alternatives)[altInd].occ;
}

char Atom::atomInfo::getAltLocID(int altInd) const {
  if ((alternatives == NULL) || (altInd >= alternatives->size()) || (altInd < 0)) MstUtils::error("alternative index " + MstUtils::toString(altInd) + " out of bounds (" + MstUtils::toString(alternatives->size()) + " alternatives available)", "Atom::getAltLocID");
  return (*alternatives)[altInd].alt;
}

mstreal Atom::getMass(const char* name) {
  switch (name[0]) {
    case 'N':
      return 14.0067;
    case 'C':
      return 12.0107;
    case 'O':
      return 15.9994;
    case 'P':
      return 30.973762;
    case 'S':
      return 32.065;
    case 'H':
      return 1.00794;
    default:
      MstUtils::error("do not know mass for atom named " + string(name));
  }
  return 0; // to make the compiler happy
}

void Atom::setCoor(const CartesianPoint& xyz) {
  x = xyz[0]; y = xyz[1]; z = xyz[2];
}

void Atom::setCoor(const Atom& a) {
  x = a.getX(); y = a.getY(); z = a.getZ();
}

void Atom::atomInfo::setAltCoor(int ai, mstreal _x, mstreal _y, mstreal _z) {
  altInfo& A = (*alternatives)[ai];
  A.x = _x; A.y = _y; A.z = _z;
}

void Atom::swapWithAlternative(int altInd) {
  if (!hasInfo() || (altInd >= numAlternatives())) MstUtils::error("alternative index " + MstUtils::toString(altInd) + " out of bounds (" + MstUtils::toString(numAlternatives()) + " alternatives available)", "Atom::swapWithAlternative");
  atomInfo::altInfo& targ = (*(info->alternatives))[altInd];
  atomInfo::altInfo temp = targ;
  targ.x = x; targ.y = y; targ.z = z; targ.occ = getOcc(); targ.B = getB(); targ.alt = getAlt();
  x = temp.x; y = temp.y; z = temp.z; setOcc(temp.occ); setB(temp.B); setAlt(temp.alt);
}

void Atom::makeAlternativeMain(int altInd) {
  if (!hasInfo() || (altInd >= numAlternatives())) MstUtils::error("alternative index " + MstUtils::toString(altInd) + " out of bounds (" + MstUtils::toString(numAlternatives()) + " alternatives available)", "Atom::makeAlternativeMain");
  atomInfo::altInfo& targ = (*(info->alternatives))[altInd];
  x = targ.x; y = targ.y; z = targ.z; setOcc(targ.occ); setB(targ.B); setAlt(targ.alt);
}

void Atom::atomInfo::addAlternative(mstreal _x, mstreal _y, mstreal _z, mstreal _B, mstreal _occ, char _alt) {
  if (alternatives == NULL) {
    alternatives = new vector<altInfo>(0);
  }
  alternatives->push_back(altInfo(_x, _y, _z, _B, _occ, _alt));
}

void Atom::atomInfo::removeLastAlternative() {
  alternatives->pop_back();
}

void Atom::atomInfo::removeAlternative(int i) {
  alternatives->erase(alternatives->begin() + i);
}

void Atom::atomInfo::clearAlternatives() {
  if (alternatives != NULL) { delete alternatives; alternatives = NULL; }
}

string Atom::pdbLine(int resIndex, int atomIndex) {
  char line[100]; // a PDB line is at most 80 characters, so this is plenty
  string resname = "UNK"; string chainID = "?"; string segID = "?";
  int resnum = 1; char icode = ' ';

  // chain and residue info
  Residue* parent = getParent();
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
  if (strlen(getNameC()) < 4) { sprintf(atomname, " %-.3s", getNameC()); }
  else { sprintf(atomname, "%.4s", getNameC()); }

  // moduli are used to make sure numbers do not go over prescribe field widths (this is not enforced by sprintf like with strings)
  sprintf(line, "%6s%5d %-4s%c%-4s%.1s%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %.4s",
          isHetero() ? "HETATM" : "ATOM  ", atomIndex % 100000, atomname, getAlt(), resname.c_str(), chainID.c_str(),
          resnum % 10000, icode, x, y, z, getOcc(), getB(), segID.c_str());

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

void Atom::write(ostream& _os) const {
  MstUtils::writeBin(_os, x); MstUtils::writeBin(_os, y); MstUtils::writeBin(_os, z);
  MstUtils::writeBin(_os, getOcc()); MstUtils::writeBin(_os, getB());
  MstUtils::writeBin(_os, getName()); MstUtils::writeBin(_os, getAlt());
  MstUtils::writeBin(_os, isHetero()); MstUtils::writeBin(_os, getIndex());
  if (numAlternatives() == 0) {
    MstUtils::writeBin(_os, false);
  } else {
    MstUtils::writeBin(_os, true);
    MstUtils::writeBin(_os, numAlternatives());
    for (int i = 0; i < numAlternatives(); i++) {
      CartesianPoint coor = getAltCoor(i);
      MstUtils::writeBin(_os, coor[0]); MstUtils::writeBin(_os, coor[1]); MstUtils::writeBin(_os, coor[2]);
      MstUtils::writeBin(_os, getAltOcc(i)); MstUtils::writeBin(_os, getAltB(i)); MstUtils::writeBin(_os, getAltLocID(i));
    }
  }
}

void Atom::read(istream& _os) {
  MstUtils::readBin(_os, x); MstUtils::readBin(_os, y); MstUtils::readBin(_os, z);
  MstUtils::readBin(_os, info->occ); MstUtils::readBin(_os, info->B);
  string _name; MstUtils::readBin(_os, _name); setName(_name); MstUtils::readBin(_os, info->alt);
  MstUtils::readBin(_os, info->het); MstUtils::readBin(_os, info->index);
  bool hasAlt; MstUtils::readBin(_os, hasAlt);
  clearAlternatives();
  if (hasAlt) {
    int nAlts; MstUtils::readBin(_os, nAlts);
    for (int i = 0; i < nAlts; i++) {
      mstreal aX, aY, aZ, aO, aB; char aA;
      MstUtils::readBin(_os, aX); MstUtils::readBin(_os, aY); MstUtils::readBin(_os, aZ);
      MstUtils::readBin(_os, aO); MstUtils::readBin(_os, aB); MstUtils::readBin(_os, aA);
      addAlternative(aX, aY, aZ, aB, aO, aA);
    }
  }
}

ostream& MST::operator<<(ostream &_os, const Atom& _atom) {
  _os << _atom.getName() << _atom.getAlt() << " " << _atom.getIndex() << " " << (_atom.isHetero() ? "HETERO" : "");
  _os << _atom.getX() << " " << _atom.getY() << " " << _atom.getZ() << " : " << _atom.getOcc() << " " << _atom.getB();
  return _os;
}

void Atom::stripInfo() {
  if (info != NULL) {
    delete info;
    info = NULL;
  }
}

/* --------- AtomPointerVector --------- */

void AtomPointerVector::copyCoordinates(const AtomPointerVector& other) {
  if (size() != other.size()) MstUtils::error("vector sizes disagree", "AtomPointerVector::copyCoordinates");
  for (int i = 0; i < size(); i++) {
    (*this)[i]->setCoor(other[i]);
  }
}

void AtomPointerVector::push_back(const Residue& R) {
  int sz = this->size();
  this->resize(sz + R.atomSize());
  for (int i = 0; i < R.atomSize(); i++) {
    (*this)[i + sz] = &(R[i]);
  }
}

void AtomPointerVector::push_back(const vector<Atom*>& atoms) {
  int sz = this->size();
  this->resize(sz + atoms.size());
  for (int i = 0; i < atoms.size(); i++) {
    (*this)[i + sz] = atoms[i];
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

AtomPointerVector AtomPointerVector::clone() const {
  AtomPointerVector into;
  clone(into);
  return into;
}

AtomPointerVector AtomPointerVector::subvector(int beg, int end) {
  return AtomPointerVector(vector<Atom*>(begin() + beg, begin() + end));
}

void AtomPointerVector::clone(AtomPointerVector& into) const {
  int L = into.size();
  into.resize(L + size());
  for (int i = 0; i < size(); i++) {
    Atom* newAtom = new Atom(*((*this)[i]));
    newAtom->setParent(NULL);
    into[L + i] = newAtom;
  }
}

AtomPointerVector& AtomPointerVector::operator+=(const AtomPointerVector& rhs) {
  if (size() != rhs.size()) MstUtils::error("vector sizes disagree", "AtomPointerVector::operator+=");
  for (int i = 0; i < size(); i++) {
    Atom& A = *((*this)[i]);
    const Atom& B = *(rhs[i]);
    A.setCoor(A.getX() + B.getX(), A.getY() + B.getY(), A.getZ() + B.getZ());
  }
  return *this;
}

AtomPointerVector& AtomPointerVector::operator-=(const AtomPointerVector& rhs) {
  if (size() != rhs.size()) MstUtils::error("vector sizes disagree", "AtomPointerVector::operator+=");
  for (int i = 0; i < size(); i++) {
    Atom& A = *((*this)[i]);
    Atom& B = *(rhs[i]);
    A.setCoor(A.getX() - B.getX(), A.getY() - B.getY(), A.getZ() - B.getZ());
  }
  return *this;
}

AtomPointerVector& AtomPointerVector::operator/=(const mstreal& s) {
  for (int i = 0; i < size(); i++) {
    (*this)[i]->setCoor((*this)[i]->getX()/s, (*this)[i]->getY()/s, (*this)[i]->getZ()/s);
  }
  return *this;
}

AtomPointerVector& AtomPointerVector::operator*=(const mstreal& s) {
  for (int i = 0; i < size(); i++) {
    (*this)[i]->setCoor((*this)[i]->getX()*s, (*this)[i]->getY()*s, (*this)[i]->getZ()*s);
  }
  return *this;
}


void AtomPointerVector::write(ostream &_os) const {
  MstUtils::writeBin(_os, (int) size());
  for (int i = 0; i < size(); i++) (*this)[i]->write(_os);
}

void AtomPointerVector::read(istream &_is) {
  int len; MstUtils::readBin(_is, len); resize(len, NULL);
  for (int i = 0; i < size(); i++) {
    (*this)[i] = new Atom();
    (*this)[i]->read(_is);
  }
}

ostream& MST::operator<<(ostream &_os, const AtomPointerVector& _atoms) {
  for (int i = 0; i < _atoms.size(); i++) {
    _os << *(_atoms[i]) << endl;
  }
  return _os;
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

CartesianPoint& CartesianPoint::operator+=(const mstreal& d) {
  for (int i = 0; i < size(); i++) {
    (*this)[i] += d;
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

const CartesianPoint CartesianPoint::operator+(const mstreal& d) const {
  CartesianPoint result = *this;
  result += d;
  return result;
}

const CartesianPoint CartesianPoint::operator/(const mstreal& s) const {
  CartesianPoint result = *this;
  result /= s;
  return result;
}

// CartesianPoint& CartesianPoint::operator=(const Atom& A) {
//   *this = CartesianPoint(A);
//   return *this;
// }

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

mstreal CartesianPoint::stdev() const {
  return sqrt(var());
}

mstreal CartesianPoint::var() const {
  mstreal m = mean();
  return norm2()/size() - m*m;
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

CartesianPoint CartesianPoint::cross(const CartesianPoint& other) const {
  if (size() != 3) MstUtils::error("don't know how to compute cross produces for dimensions other than 3", "CartesianPoint::cross");
  if (size() != other.size()) MstUtils::error("vector size mismatch", "CartesianPoint::cross");

  CartesianPoint C(3, 0);
  C[0] = getY()*other.getZ() - getZ()*other.getY();
  C[1] = getZ()*other.getX() - getX()*other.getZ();
  C[2] = getX()*other.getY() - getY()*other.getX();

  return C;
}

CartesianPoint CartesianPoint::elemProd(const CartesianPoint& other) const {
  if (size() != other.size()) MstUtils::error("vector size mismatch", "CartesianPoint::elemProd");
  CartesianPoint P(size());
  for (int i = 0; i < size(); i++) P[i] = (*this)[i] * other[i];
  return P;
}

mstreal CartesianPoint::dot(const CartesianPoint& other) const {
  if (size() != other.size()) MstUtils::error("vector size mismatch", "CartesianPoint::dot");

  mstreal d = 0;
  for (int i = 0; i < size(); i++) {
    d += (*this)[i] * other[i];
  }

  return d;
}

void CartesianPoint::convertToSphericalCoordinates(mstreal& radius, mstreal& polarAngle, mstreal& azimuthalAngle) {
  if (this->size() != 3) MstUtils::error("Can only convert 3D vector to spherical coordinates","CartesianPoint::convertToSphericalCoordinates");
  radius = this->norm();
  polarAngle = acos((*this)[2]/radius);
  azimuthalAngle = fmod(atan2((*this)[1],(*this)[0]) + 2*M_PI, 2*M_PI);
}

void CartesianPoint::setPositionBySphericalCoordinates(mstreal radius, mstreal polarAngle, mstreal azimuthalAngle) {
  this->resize(3, 0);
  (*this)[0] = radius*cos(azimuthalAngle)*sin(polarAngle);
  (*this)[1] = radius*sin(azimuthalAngle)*sin(polarAngle);
  (*this)[2] = radius*cos(polarAngle);
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

mstreal CartesianPoint::cosineAngle(const CartesianPoint& another) const {
  if (this->size() != another.size()) MstUtils::error("point dimensions disagree", "CartesianPoint::cosineAngle");
  return (*this).dot(another) / ((*this).norm()* another.norm());
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
      if (CartesianGeometry::angleDiffCCW(angles[j], angles[minInd], radians) > arc) {
        arc = CartesianGeometry::angleDiffCCW(angles[j], angles[minInd], radians);
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

// see derivation in http://grigoryanlab.org/docs/dynamics_derivatives.pdf
template <class T>
mstreal CartesianGeometry::distance(const CartesianPoint& atom1, const CartesianPoint& atom2, T& grad) {
  mstreal x12 = atom1.getX() - atom2.getX();
  mstreal y12 = atom1.getY() - atom2.getY();
  mstreal z12 = atom1.getZ() - atom2.getZ();
  mstreal d = sqrt(x12*x12 + y12*y12 + z12*z12);
  if (MstUtils::closeEnough(d, 0.0)) {
    // by convention, set gradient to unity when the two atoms are on top of each other
    for (int i = 0; i < grad.size(); i++) grad[i] = 1.0;
  } else {
    grad[0] = x12/d;
    grad[1] = y12/d;
    grad[2] = z12/d;
    grad[3] = -grad[0];
    grad[4] = -grad[1];
    grad[5] = -grad[2];
  }
  return d;
}

// see derivation in http://grigoryanlab.org/docs/dynamics_derivatives.pdf
template <class T>
mstreal CartesianGeometry::angle(const CartesianPoint& atom1, const CartesianPoint& atom2, const CartesianPoint& atom3, T& grad, bool radians) {
  mstreal x12 = atom1.getX() - atom2.getX();
  mstreal y12 = atom1.getY() - atom2.getY();
  mstreal z12 = atom1.getZ() - atom2.getZ();
  mstreal x32 = atom3.getX() - atom2.getX();
  mstreal y32 = atom3.getY() - atom2.getY();
  mstreal z32 = atom3.getZ() - atom2.getZ();
  mstreal L1 = sqrt(x12*x12 + y12*y12 + z12*z12);
  mstreal L2 = sqrt(x32*x32 + y32*y32 + z32*z32);
  if (MstUtils::closeEnough(L1, 0.0) || MstUtils::closeEnough(L2, 0.0)) {
    // angle is not really defined, when two subsequent atoms coicide, so just
    // set the derivative to something for the sake of code stability
    for (int i = 0; i < grad.size(); i++) grad[i] = 0.0;
  }
  mstreal p = x12*x32 + y12*y32 + z12*z32;
  mstreal d = p/(L1*L2);
  mstreal an = acos(d);
  if (MstUtils::closeEnough(fabs(d), 1.0)) {
    // signs are arbitrary in this case (set to + by convention)
    grad[0] = (sqrt(y12*y12 + z12*z12)/L1)/L1;
    grad[1] = (sqrt(x12*x12 + z12*z12)/L1)/L1;
    grad[2] = (sqrt(x12*x12 + y12*y12)/L1)/L1;
    grad[6] = (sqrt(y32*y32 + z32*z32)/L2)/L2;
    grad[7] = (sqrt(x32*x32 + z32*z32)/L2)/L2;
    grad[8] = (sqrt(x32*x32 + y32*y32)/L2)/L2;
  } else {
    mstreal L12 = L1*L2;
    mstreal C = -1/(sqrt(1 - d*d)*L12*L12);
    mstreal pL12i = p*L1/L2;
    mstreal pL21i = p*L2/L1;
    grad[0] = C*(x32*L12 - pL21i*x12);
    grad[1] = C*(y32*L12 - pL21i*y12);
    grad[2] = C*(z32*L12 - pL21i*z12);
    grad[6] = C*(x12*L12 - pL12i*x32);
    grad[7] = C*(y12*L12 - pL12i*y32);
    grad[8] = C*(z12*L12 - pL12i*z32);
  }

  // gardient with respect to the middle point is always the minus sum of that
  // of the first and the third points
  grad[3] = -(grad[0] + grad[6]);
  grad[4] = -(grad[1] + grad[7]);
  grad[5] = -(grad[2] + grad[8]);

  if (!radians) {
    mstreal f = 180.0/M_PI;
    d *= f;
    for (int i = 0; i < grad.size(); i++) grad[i] *= f;
  }
  return an;
}

// see derivation in http://grigoryanlab.org/docs/dynamics_derivatives.pdf
template <class T>
mstreal CartesianGeometry::dihedral(const CartesianPoint& atom1, const CartesianPoint& atom2, const CartesianPoint& atom3, const CartesianPoint& atom4, T& grad, bool radians) {
  mstreal x21 = atom2.getX() - atom1.getX();
  mstreal y21 = atom2.getY() - atom1.getY();
  mstreal z21 = atom2.getZ() - atom1.getZ();
  mstreal x32 = atom3.getX() - atom2.getX();
  mstreal y32 = atom3.getY() - atom2.getY();
  mstreal z32 = atom3.getZ() - atom2.getZ();
  mstreal x43 = atom4.getX() - atom3.getX();
  mstreal y43 = atom4.getY() - atom3.getY();
  mstreal z43 = atom4.getZ() - atom3.getZ();
  mstreal x31 = atom3.getX() - atom1.getX();
  mstreal y31 = atom3.getY() - atom1.getY();
  mstreal z31 = atom3.getZ() - atom1.getZ();
  mstreal x42 = atom4.getX() - atom2.getX();
  mstreal y42 = atom4.getY() - atom2.getY();
  mstreal z42 = atom4.getZ() - atom2.getZ();

  CartesianPoint N1(z21*y32 - y21*z32, x21*z32 - z21*x32, y21*x32 - x21*y32);
  CartesianPoint N2(y43*z32 - z43*y32, z43*x32 - x43*z32, x43*y32 - y43*x32);
  CartesianPoint angleGrad(9);
  mstreal th = CartesianGeometry::angle(N1, CartesianPoint(0, 0, 0), N2, angleGrad, radians);

  grad[0]  = -angleGrad[1]*z32 + angleGrad[2]*y32;
  grad[1]  =  angleGrad[0]*z32 - angleGrad[2]*x32;
  grad[2]  = -angleGrad[0]*y32 + angleGrad[1]*x32;

  grad[3]  =  angleGrad[1]*z31 - angleGrad[2]*y31 - angleGrad[7]*z43 + angleGrad[8]*y43;
  grad[4]  = -angleGrad[0]*z31 + angleGrad[2]*x31 + angleGrad[6]*z43 - angleGrad[8]*x43;
  grad[5]  =  angleGrad[0]*y31 - angleGrad[1]*x31 - angleGrad[6]*y43 + angleGrad[7]*x43;

  grad[6]  = -angleGrad[1]*z21 + angleGrad[2]*y21 + angleGrad[7]*z42 - angleGrad[8]*y42;
  grad[7]  =  angleGrad[0]*z21 - angleGrad[2]*x21 - angleGrad[6]*z42 + angleGrad[8]*x42;
  grad[8]  = -angleGrad[0]*y21 + angleGrad[1]*x21 + angleGrad[6]*y42 - angleGrad[7]*x42;

  grad[9]  = -angleGrad[7]*z32 + angleGrad[8]*y32;
  grad[10] =  angleGrad[6]*z32 - angleGrad[8]*x32;
  grad[11] = -angleGrad[6]*y32 + angleGrad[7]*x32;

  if (N1 * CartesianPoint(x43, y43, z43) > 0) {
    for (int i = 0; i < grad.size(); i++) grad[i] = -grad[i];
    th = -th;
  }

  return th;
}

bool CartesianGeometry::testPrimitiveGradients(bool radians) {
  long int x = std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
  srand(x);

  // test a bunch of times
  CartesianPoint bondGrad(6, 0.0), angleGrad(9, 0.0), diheGrad(12, 0.0);
  CartesianPoint bondGradFD(6, 0.0), angleGradFD(9, 0.0), diheGradFD(12, 0.0);
  for (int k = 0; k < 100; k++) {
    mstreal L = 1 + MstUtils::randUnit()*10; // length scale
    mstreal del = 0.0001; // finite-difference coordinate step
    mstreal tol = 0.01;   // tollerated norm difference between analytical and FD

    // make some random points
    vector<CartesianPoint> P;
    for (int i = 0; i < 4; i++) {
      P.push_back(CartesianPoint(L*MstUtils::randUnit(), L*MstUtils::randUnit(), L*MstUtils::randUnit()));
    }

    // test bond derivatives
    CartesianGeometry::distance(P[0], P[1], bondGrad);
    bondGradFD[0] = P[1].distance(P[0] + CartesianPoint(del, 0, 0)) - P[1].distance(P[0] - CartesianPoint(del, 0, 0));
    bondGradFD[1] = P[1].distance(P[0] + CartesianPoint(0, del, 0)) - P[1].distance(P[0] - CartesianPoint(0, del, 0));
    bondGradFD[2] = P[1].distance(P[0] + CartesianPoint(0, 0, del)) - P[1].distance(P[0] - CartesianPoint(0, 0, del));
    bondGradFD[3] = P[0].distance(P[1] + CartesianPoint(del, 0, 0)) - P[0].distance(P[1] - CartesianPoint(del, 0, 0));
    bondGradFD[4] = P[0].distance(P[1] + CartesianPoint(0, del, 0)) - P[0].distance(P[1] - CartesianPoint(0, del, 0));
    bondGradFD[5] = P[0].distance(P[1] + CartesianPoint(0, 0, del)) - P[0].distance(P[1] - CartesianPoint(0, 0, del));
    bondGradFD /= 2*del;
    if ((bondGrad - bondGradFD).norm() > tol) {
      cout << "test FAILED for distance gradient, point set: " << MstUtils::vecToString(P, " | ") << endl;
      cout << "analytical: " << MstUtils::vecToString(bondGrad) << endl;
      cout << "numerical: " << MstUtils::vecToString(bondGradFD) << endl;
      return false;
    }

    // test angle derivatives
    CartesianGeometry::angle(P[0], P[1], P[2], angleGrad, radians);
    angleGradFD[0] = CartesianGeometry::angle(P[0] + CartesianPoint(del, 0, 0), P[1], P[2], radians) - CartesianGeometry::angle(P[0] - CartesianPoint(del, 0, 0), P[1], P[2], radians);
    angleGradFD[1] = CartesianGeometry::angle(P[0] + CartesianPoint(0, del, 0), P[1], P[2], radians) - CartesianGeometry::angle(P[0] - CartesianPoint(0, del, 0), P[1], P[2], radians);
    angleGradFD[2] = CartesianGeometry::angle(P[0] + CartesianPoint(0, 0, del), P[1], P[2], radians) - CartesianGeometry::angle(P[0] - CartesianPoint(0, 0, del), P[1], P[2], radians);
    angleGradFD[3] = CartesianGeometry::angle(P[0], P[1] + CartesianPoint(del, 0, 0), P[2], radians) - CartesianGeometry::angle(P[0], P[1] - CartesianPoint(del, 0, 0), P[2], radians);
    angleGradFD[4] = CartesianGeometry::angle(P[0], P[1] + CartesianPoint(0, del, 0), P[2], radians) - CartesianGeometry::angle(P[0], P[1] - CartesianPoint(0, del, 0), P[2], radians);
    angleGradFD[5] = CartesianGeometry::angle(P[0], P[1] + CartesianPoint(0, 0, del), P[2], radians) - CartesianGeometry::angle(P[0], P[1] - CartesianPoint(0, 0, del), P[2], radians);
    angleGradFD[6] = CartesianGeometry::angle(P[0], P[1], P[2] + CartesianPoint(del, 0, 0), radians) - CartesianGeometry::angle(P[0], P[1], P[2] - CartesianPoint(del, 0, 0), radians);
    angleGradFD[7] = CartesianGeometry::angle(P[0], P[1], P[2] + CartesianPoint(0, del, 0), radians) - CartesianGeometry::angle(P[0], P[1], P[2] - CartesianPoint(0, del, 0), radians);
    angleGradFD[8] = CartesianGeometry::angle(P[0], P[1], P[2] + CartesianPoint(0, 0, del), radians) - CartesianGeometry::angle(P[0], P[1], P[2] - CartesianPoint(0, 0, del), radians);

    angleGradFD /= 2*del;
    if ((angleGrad - angleGradFD).norm() > tol) {
      cout << "test FAILED for angle gradient, point set: " << MstUtils::vecToString(P, " | ") << endl;
      cout << "analytical: " << MstUtils::vecToString(angleGrad) << endl;
      cout << "numerical: " << MstUtils::vecToString(angleGradFD) << endl;
      return false;
    }

    // test dihedral derivatives
    CartesianGeometry::dihedral(P[0], P[1], P[2], P[3], diheGrad, radians);
    diheGradFD[0] =  CartesianGeometry::dihedral(P[0] + CartesianPoint(del, 0, 0), P[1], P[2], P[3], radians) - CartesianGeometry::dihedral(P[0] - CartesianPoint(del, 0, 0), P[1], P[2], P[3], radians);
    diheGradFD[1] =  CartesianGeometry::dihedral(P[0] + CartesianPoint(0, del, 0), P[1], P[2], P[3], radians) - CartesianGeometry::dihedral(P[0] - CartesianPoint(0, del, 0), P[1], P[2], P[3], radians);
    diheGradFD[2] =  CartesianGeometry::dihedral(P[0] + CartesianPoint(0, 0, del), P[1], P[2], P[3], radians) - CartesianGeometry::dihedral(P[0] - CartesianPoint(0, 0, del), P[1], P[2], P[3], radians);
    diheGradFD[3] =  CartesianGeometry::dihedral(P[0], P[1] + CartesianPoint(del, 0, 0), P[2], P[3], radians) - CartesianGeometry::dihedral(P[0], P[1] - CartesianPoint(del, 0, 0), P[2], P[3], radians);
    diheGradFD[4] =  CartesianGeometry::dihedral(P[0], P[1] + CartesianPoint(0, del, 0), P[2], P[3], radians) - CartesianGeometry::dihedral(P[0], P[1] - CartesianPoint(0, del, 0), P[2], P[3], radians);
    diheGradFD[5] =  CartesianGeometry::dihedral(P[0], P[1] + CartesianPoint(0, 0, del), P[2], P[3], radians) - CartesianGeometry::dihedral(P[0], P[1] - CartesianPoint(0, 0, del), P[2], P[3], radians);
    diheGradFD[6] =  CartesianGeometry::dihedral(P[0], P[1], P[2] + CartesianPoint(del, 0, 0), P[3], radians) - CartesianGeometry::dihedral(P[0], P[1], P[2] - CartesianPoint(del, 0, 0), P[3], radians);
    diheGradFD[7] =  CartesianGeometry::dihedral(P[0], P[1], P[2] + CartesianPoint(0, del, 0), P[3], radians) - CartesianGeometry::dihedral(P[0], P[1], P[2] - CartesianPoint(0, del, 0), P[3], radians);
    diheGradFD[8] =  CartesianGeometry::dihedral(P[0], P[1], P[2] + CartesianPoint(0, 0, del), P[3], radians) - CartesianGeometry::dihedral(P[0], P[1], P[2] - CartesianPoint(0, 0, del), P[3], radians);
    diheGradFD[9] =  CartesianGeometry::dihedral(P[0], P[1], P[2], P[3] + CartesianPoint(del, 0, 0), radians) - CartesianGeometry::dihedral(P[0], P[1], P[2], P[3] - CartesianPoint(del, 0, 0), radians);
    diheGradFD[10] = CartesianGeometry::dihedral(P[0], P[1], P[2], P[3] + CartesianPoint(0, del, 0), radians) - CartesianGeometry::dihedral(P[0], P[1], P[2], P[3] - CartesianPoint(0, del, 0), radians);
    diheGradFD[11] = CartesianGeometry::dihedral(P[0], P[1], P[2], P[3] + CartesianPoint(0, 0, del), radians) - CartesianGeometry::dihedral(P[0], P[1], P[2], P[3] - CartesianPoint(0, 0, del), radians);
    diheGradFD /= 2*del;
    if ((diheGrad - diheGradFD).norm() > tol) {
      cout << "test FAILED for dihedral gradient, point set: " << MstUtils::vecToString(P, " | ") << endl;
      cout << "analytical: " << MstUtils::vecToString(diheGrad) << endl;
      cout << "numerical: " << MstUtils::vecToString(diheGradFD) << endl;
      cout << "difference: " << (diheGrad - diheGradFD) << endl;
      cout << "norm difference: " << (diheGrad - diheGradFD).norm() << endl;
      return false;
    }
  }
  return true;
}

// forward declarations of template functions
template mstreal CartesianGeometry::distance<vector<mstreal> >(const CartesianPoint& atom1, const CartesianPoint& atom2, vector<mstreal>& grad);
template mstreal CartesianGeometry::angle<vector<mstreal> >(const CartesianPoint& atom1, const CartesianPoint& atom2, const CartesianPoint& atom3, vector<mstreal>& grad, bool radians);
template mstreal CartesianGeometry::dihedral<vector<mstreal> >(const CartesianPoint& atom1, const CartesianPoint& atom2, const CartesianPoint& atom3, const CartesianPoint& atom4, vector<mstreal>& grad, bool radians);

/* --------- selector --------------- */

selector::selector(const Structure& S) {
  atoms = S.getAtoms();
  residues = S.getResidues();
  atomResidues.resize(atoms.size());
  atomChains.resize(atoms.size());
  for (int i = 0; i < atoms.size(); i++) {
    atomResidues[i] = atoms[i]->getParent();
    if (atomResidues[i] != NULL) atomChains[i] = atomResidues[i]->getParent();
    if ((atomResidues[i] == NULL) || (atomChains[i] == NULL)) MstUtils::error("internally inconsistent Structure given", "selector::selector");
  }
}

AtomPointerVector selector::select(string selStr) {
  expressionTree* tree = buildExpressionTree(selStr);
  AtomPointerVector sel;
  select(tree, sel);
  delete tree;
  set<Atom*> selAtoms(sel.begin(), sel.end());
  AtomPointerVector sortedSel;
  for (int i = 0; i < atoms.size(); i++) {
    if (selAtoms.find(atoms[i]) != selAtoms.end()) sortedSel.push_back(atoms[i]);
  }
  return sortedSel;
}

vector<Residue*> selector::selectRes(string selStr) {
  AtomPointerVector sel = select(selStr);
  set<Residue*> selRes;
  for (int i = 0; i < sel.size(); i++) {
    selRes.insert(sel[i]->getParent());
  }
  vector<Residue*> selResSorted;
  for (int i = 0; i < residues.size(); i++) {
    if (selRes.find(residues[i]) != selRes.end()) selResSorted.push_back(residues[i]);
  }
  return selResSorted;
}

void selector::select(expressionTree* tree, AtomPointerVector& sel) {
  if (tree->numChildren() == 0) {
    // this is a terminal node, so just do the selection
    for (int i = 0; i < atoms.size(); i++) {
      switch(tree->getProperty()) {
        case (expressionTree::selProperty::ALL):
          sel.push_back(atoms[i]);
          break;
        case (expressionTree::selProperty::RESID):
          if ((tree->hasNumSet() && tree->inNumSet(atomResidues[i]->getNum())) ||
              (!tree->hasNumSet() && tree->hasNums() && (atomResidues[i]->getNum() >= tree->getNumByIdx(0)) && (atomResidues[i]->getNum() <= tree->getNumByIdx(1))) ||
              (!tree->hasNumSet() && !tree->hasNums() && (atomResidues[i]->getNum() == tree->getNum()))) {
            sel.push_back(atoms[i]);
          }
          break;
        case (expressionTree::selProperty::RESNAME):
          if (atomResidues[i]->getName() == tree->getString()) sel.push_back(atoms[i]);
          break;
        case (expressionTree::selProperty::ICODE):
          if (string(1, atomResidues[i]->getIcode()) == tree->getString()) sel.push_back(atoms[i]);
          break;
        case (expressionTree::selProperty::CHAIN):
          if (atomChains[i]->getID() == tree->getString()) sel.push_back(atoms[i]);
          break;
        case (expressionTree::selProperty::SEGID):
          if (atomChains[i]->getSegID() == tree->getString()) sel.push_back(atoms[i]);
          break;
        case (expressionTree::selProperty::NAME):
          if (atoms[i]->getName() == tree->getString()) sel.push_back(atoms[i]);
          break;
        default:
          MstUtils::error("unknown selectable property " + MstUtils::toString(tree->getProperty()), "selector::select");
      }
    }
  } else {
    AtomPointerVector selA, selB;
    switch(tree->getLogicalOperator()) {
      case (expressionTree::logicalOp::AND):
        if (tree->numChildren() != 2)
          MstUtils::error("poorly parsed expression: expected two operands for AND", "selector::select(expressionTree* )");
        select(tree->getChild(0), selA);
        select(tree->getChild(1), selB);
        sel = intersect(selA, selB);
        break;
      case (expressionTree::logicalOp::OR):
        if (tree->numChildren() != 2)
          MstUtils::error("poorly parsed expression: expected two operands for OR", "selector::select(expressionTree* )");
        select(tree->getChild(0), selA);
        select(tree->getChild(1), selB);
        sel = combine(selA, selB);
        break;
      case (expressionTree::logicalOp::NOT):
        if (tree->numChildren() != 1)
          MstUtils::error("poorly parsed expression: expected one operand for NOT", "selector::select(expressionTree* )");
        select(tree->getChild(0), selA);
        sel = invert(selA);
        break;
      case (expressionTree::logicalOp::BYRES):
        if (tree->numChildren() != 1)
          MstUtils::error("poorly parsed expression: expected one operand for BYRES", "selector::select(expressionTree* )");
        select(tree->getChild(0), selA);
        sel = byRes(selA);
        break;
      case (expressionTree::logicalOp::BYCHAIN):
        if (tree->numChildren() != 1)
          MstUtils::error("poorly parsed expression: expected one operand for BYCHAIN", "selector::select(expressionTree* )");
        select(tree->getChild(0), selA);
        sel = byChain(selA);
        break;
      case (expressionTree::logicalOp::IS):
        if (tree->numChildren() != 1)
          MstUtils::error("poorly parsed expression: expected one operand for IS", "selector::select(expressionTree* )");
        select(tree->getChild(0), sel);
        break;
      case (expressionTree::logicalOp::AROUND):
        if (tree->numChildren() != 1)
          MstUtils::error("poorly parsed expression: expected one selection operand for AROUND", "selector::select(expressionTree* )");
        select(tree->getChild(0), selA);
        sel = around(selA, tree->getVal());
        break;
      default:
        MstUtils::error("unknown selectable property " + MstUtils::toString(tree->getProperty()), "selector::select");
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
    tree->setLogicalOperator(expressionTree::logicalOp::IS);
    tree->setProperty(expressionTree::selProperty::RESID);
    if (MstUtils::isInt(str)) {
      tree->setNum(MstUtils::toInt(str));
    } else {
      vector<string> parts = MstUtils::split(str, "+");
      for (const string& part : parts) {
        vector<string> range = MstUtils::split(part, "-");
        int lower, upper;
        if (range.size() == 1) {
          if (!MstUtils::isInt(part)) MstUtils::error("bad selection, expected an integer, not " + part, "selector::buildExpressionTree(string)");
          lower = MstUtils::toInt(part);
          upper = lower;
        } else if (range.size() == 2) {
          if (!MstUtils::isInt(range[0]) || !MstUtils::isInt(range[1])) {
            MstUtils::error("bad selection, a range should be two integers separated by a dash, not " + part, "selector::buildExpressionTree(string)");
          }
          lower = MstUtils::toInt(range[0]);
          upper = MstUtils::toInt(range[1]);
        } else MstUtils::error("bad selection, each token separated by '+' should be an integer or a range of integers, not " + part, "selector::buildExpressionTree(string)");
        tree->addToNumSet(MstUtils::range(lower, upper + 1));
      }
    }
  } else if (MstUtils::stringsEqual(token, "all")) {
    tree->setLogicalOperator(expressionTree::logicalOp::IS);
    tree->setProperty(expressionTree::selProperty::ALL);
  } else if (MstUtils::stringsEqual(token, "resname")) {
    string str = getNextSelectionToken(selStr);
    tree->setLogicalOperator(expressionTree::logicalOp::IS);
    tree->setProperty(expressionTree::selProperty::RESNAME);
    tree->setString(str);
  } else if (MstUtils::stringsEqual(token, "icode")) {
    string str = getNextSelectionToken(selStr);
    tree->setLogicalOperator(expressionTree::logicalOp::IS);
    tree->setProperty(expressionTree::selProperty::ICODE);
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
    MstUtils::assertCond(extra.empty(), "poorly formed selection expression: extra stuff after an 'around' operator; use parentheses as necessary.");
    return root;
  } else {
    MstUtils::error("bad selection, unrecognized connector keyword '" + connector + "' after token '" + token + "'", "selector::buildExpressionTree(string)");
  }

  root->addChild(buildExpressionTree(selStr));
  return root;
}

string selector::getNextSelectionToken(string& selStr) {
  selStr = MstUtils::trim(selStr);
  if (selStr.empty()) return "";

  if (selStr[0] == '(') { // if parenthetical, find matching paren
    // find matching closing paren
    int n = 1; int i;
    char quoteChar = '\0';
    for (i = 1; i < selStr.size(); i++) {
      if (quoteChar == '\0') {
        if (selStr[i] == '(') n++;
        else if (selStr[i] == ')') n--;
        else if ((selStr[i] == '\'') && (selStr[i-1] != '\\')) quoteChar = '\'';
        else if ((selStr[i] == '"') && (selStr[i-1] != '\\')) quoteChar = '"';
      } else if ((selStr[i] == quoteChar) && (selStr[i-1] != '\\')) quoteChar = '\0';
      if (n == 0) break;
    }
    if (n != 0) MstUtils::error("ill-formed selection expression '" + selStr + "'", "Structure::getNextSelectionToken");
    string token = selStr.substr(0, i+1); // keep trailing ()
    selStr = selStr.substr(i+1);
    return token;
  }

  // otherwise, just get the next space-delimited token
  string token = MstUtils::nextQuoteAwareToken(selStr);
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

vector<mstreal> RMSDCalculator::bestRMSD(const vector<vector<Atom*>> &_alignSet, const vector<Atom*> &_ref, int start, int fromEnd) {
  vector<mstreal> rmsdVector(_alignSet.size());
  fill(rmsdVector.begin(), rmsdVector.end(), -1.0);
  RMSDCalculator rc;
  for (int i = start; i < _alignSet.size() - fromEnd; i++) {
    rmsdVector[i] = rc.bestRMSD(_alignSet[i], _ref);
  }
  return rmsdVector;
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

void RMSDCalculator::applyLastTransformation(vector<Atom*>& _moveable) {
  mstreal x[3], x1[3];
  for(int k = 0; k < _moveable.size(); k++) {
      x[0] = _moveable[k]->getX();
      x[1] = _moveable[k]->getY();
      x[2] = _moveable[k]->getZ();
      x1[0] = t[0]+u[0][0]*x[0]+u[0][1]*x[1]+u[0][2]*x[2];
      x1[1] = t[1]+u[1][0]*x[0]+u[1][1]*x[1]+u[1][2]*x[2];
      x1[2] = t[2]+u[2][0]*x[0]+u[2][1]*x[1]+u[2][2]*x[2];
      _moveable[k]->setCoor(x1[0], x1[1], x1[2]);
  }
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
    mstreal ax, ay, az, rx, ry, rz;
    for (m = 0; m < n; m++) {
      ax = _align[m]->getX() - xc[0];
      ay = _align[m]->getY() - xc[1];
      az = _align[m]->getZ() - xc[2];
      rx = _ref[m]->getX() - yc[0];
      ry = _ref[m]->getY() - yc[1];
      rz = _ref[m]->getZ() - yc[2];
      e0 += ax * ax + rx * rx;
      e0 += ay * ay + ry * ry;
      e0 += az * az + rz * rz;
      r[0][0] += rx * ax;
      r[0][1] += rx * ay;
      r[0][2] += rx * az;
      r[1][0] += ry * ax;
      r[1][1] += ry * ay;
      r[1][2] += ry * az;
      r[2][0] += rz * ax;
      r[2][1] += rz * ay;
      r[2][2] += rz * az;
    }

    //compute determinat of matrix r
    det = r[0][0] * ( r[1][1]*r[2][2] - r[1][2]*r[2][1] )       \
        - r[0][1] * ( r[1][0]*r[2][2] - r[1][2]*r[2][0] )       \
        + r[0][2] * ( r[1][0]*r[2][1] - r[1][1]*r[2][0] );
    sigma = det;

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

mstreal RMSDCalculator::rmsdCutoff(const vector<vector<int> >& I, mstreal rmsdMax, mstreal L0) {
  int N = 0;
  mstreal c = 0;

  // disjoint segments are counted as independent, so their correlation
  // with respect to each other is zero
  for (int i = 0; i < I.size(); i++) {
    for (int j = 0; j < I[i].size(); j++) {
      for (int k = j + 1; k < I[i].size(); k++) {
        c = c + exp(-abs(I[i][j] - I[i][k])/L0);
      }
    }
    N += I[i].size();
  }
  double df = N*(1 - (2.0/(N*(N-1)))*c);

  return rmsdMax/sqrt(N/df);
}

mstreal RMSDCalculator::rmsdCutoff(const Structure& S, mstreal rmsdMax, mstreal L0) {
  vector<int> L(S.chainSize());
  for (int i = 0; i < S.chainSize(); i++) L[i] = S[i].residueSize();
  return RMSDCalculator::rmsdCutoff(L, rmsdMax, L0);
}

mstreal RMSDCalculator::rmsdCutoff(const vector<int>& J, const Structure& S, mstreal rmsdMax, mstreal L0) {
  vector<vector<int> > I;
  map<Chain*, vector<int> > residuesFromChain;
  for (int i = 0; i < J.size(); i++) residuesFromChain[S.getResidue(J[i]).getChain()].push_back(J[i]);
  vector<Chain*> keys = MstUtils::keys(residuesFromChain);
  I.resize(keys.size());
  for (int i = 0; i < keys.size(); i++) I[i] = residuesFromChain[keys[i]];
  return RMSDCalculator::rmsdCutoff(I, rmsdMax, L0);
}

// QCP algorithm from: http://onlinelibrary.wiley.com/doi/10.1002/jcc.21439/epdf
template <class T>
mstreal RMSDCalculator::qcpRMSD(const T& A, const T& B, bool setTransform, bool setResiduals) {
  int i, j;
  int N = A.size();
  if (N != B.size()) MstUtils::error("structures are of different length", "CartesianGeometry::qcpRMSD");

  //compute centers for vector sets x, y
  mstreal cA[3], cB[3];
  for (i = 0; i < 3; i++) { cA[i] = 0.0; cB[i] = 0.0; }
  for (i = 0; i < N; i++) {
    cA[0] += A[i]->getX();
    cA[1] += A[i]->getY();
    cA[2] += A[i]->getZ();
    cB[0] += B[i]->getX();
    cB[1] += B[i]->getY();
    cB[2] += B[i]->getZ();
  }
  for (i = 0; i < 3; i++){
    cA[i] = cA[i]/N;
    cB[i] = cB[i]/N;
  }

  // compute the correlation matrix S and the inner products of the two structures
  mstreal S[3][3];
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) S[i][j] = 0.0;
  }
  mstreal ax, ay, az, bx, by, bz;
  mstreal GA = 0, GB = 0;
  for (i = 0; i < N; i++) {
    ax = A[i]->getX() - cA[0];
    ay = A[i]->getY() - cA[1];
    az = A[i]->getZ() - cA[2];
    bx = B[i]->getX() - cB[0];
    by = B[i]->getY() - cB[1];
    bz = B[i]->getZ() - cB[2];
    GA += ax*ax + ay*ay + az*az;
    GB += bx*bx + by*by + bz*bz;
    S[0][0] += bx * ax;
    S[0][1] += bx * ay;
    S[0][2] += bx * az;
    S[1][0] += by * ax;
    S[1][1] += by * ay;
    S[1][2] += by * az;
    S[2][0] += bz * ax;
    S[2][1] += bz * ay;
    S[2][2] += bz * az;
  }

  // square of S
  mstreal S2[3][3];
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) S2[i][j] = S[i][j]*S[i][j];
  }

  // calculate characteristic polynomial coefficients
  // NOTE: though there are repeated terms in F, G, H, and I, it actually turns
  // out to be better to let the compiler optimize this rathre than declare temp variables
  mstreal C2 = -2*(S2[0][0] + S2[0][1] + S2[0][2] + S2[1][0] + S2[1][1] + S2[1][2] + S2[2][0] + S2[2][1] + S2[2][2]);
  mstreal C1 = 8*(S[0][0]*S[1][2]*S[2][1] + S[1][1]*S[2][0]*S[0][2] + S[2][2]*S[0][1]*S[1][0] -
                  S[0][0]*S[1][1]*S[2][2] - S[1][2]*S[2][0]*S[0][1] - S[2][1]*S[1][0]*S[0][2]);
  mstreal D = (S2[0][1] + S2[0][2] - S2[1][0] - S2[2][0]); D = D*D;
  mstreal E1 = -S2[0][0] + S2[1][1] + S2[2][2] + S2[1][2] + S2[2][1];
  mstreal E2 = 2*(S[1][1]*S[2][2] - S[1][2]*S[2][1]);
  mstreal E = (E1 - E2) * (E1 + E2);
  mstreal F = (-(S[0][2] + S[2][0])*(S[1][2] - S[2][1]) + (S[0][1] - S[1][0])*(S[0][0] - S[1][1] - S[2][2])) *
              (-(S[0][2] - S[2][0])*(S[1][2] + S[2][1]) + (S[0][1] - S[1][0])*(S[0][0] - S[1][1] + S[2][2]));
  mstreal G = (-(S[0][2] + S[2][0])*(S[1][2] + S[2][1]) - (S[0][1] + S[1][0])*(S[0][0] + S[1][1] - S[2][2])) *
              (-(S[0][2] - S[2][0])*(S[1][2] - S[2][1]) - (S[0][1] + S[1][0])*(S[0][0] + S[1][1] + S[2][2]));
  mstreal H = ( (S[0][1] + S[1][0])*(S[1][2] + S[2][1]) + (S[0][2] + S[2][0])*(S[0][0] - S[1][1] + S[2][2])) *
              (-(S[0][1] - S[1][0])*(S[1][2] - S[2][1]) + (S[0][2] + S[2][0])*(S[0][0] + S[1][1] + S[2][2]));
  mstreal I = ( (S[0][1] + S[1][0])*(S[1][2] - S[2][1]) + (S[0][2] - S[2][0])*(S[0][0] - S[1][1] - S[2][2])) *
              (-(S[0][1] - S[1][0])*(S[1][2] + S[2][1]) + (S[0][2] - S[2][0])*(S[0][0] + S[1][1] - S[2][2]));
  mstreal C0 = D + E + F + G + H + I;

  // now iterate Newton–Raphson method to find max eigenvalue
  mstreal tol = 10E-11;
  mstreal L = (GA + GB)/2; // this initial guess is key (see http://journals.iucr.org/a/issues/2005/04/00/sh5029/sh5029.pdf)
  mstreal Lold, L2, L3, L4, C22 = 2*C2;
  do {
    Lold = L;
    L2 = L*L;
    L3 = L2*L;
    L4 = L3*L;
    L = L - (L4 + C2*L2 + C1*L + C0)/(4*L3 + C22*L + C1);
  } while (fabs(L - Lold) > tol*L);

  // compute optimal rotation matrix
  if (setTransform) {
    mstreal K[4][4];
    K[0][0] =  S[0][0] + S[1][1] + S[2][2] - L;
    K[1][1] =  S[0][0] - S[1][1] - S[2][2] - L;
    K[2][2] = -S[0][0] + S[1][1] - S[2][2] - L;
    K[3][3] = -S[0][0] - S[1][1] + S[2][2] - L;
    K[0][1] = K[1][0] = S[1][2] - S[2][1];
    K[0][2] = K[2][0] = S[2][0] - S[0][2];
    K[0][3] = K[3][0] = S[0][1] - S[1][0];
    K[1][2] = K[2][1] = S[0][1] + S[1][0];
    K[1][3] = K[3][1] = S[0][2] + S[2][0];
    K[2][3] = K[3][2] = S[1][2] + S[2][1];

    mstreal evectol = 10E-11;
    mstreal q0, q1, q2, q3; // a column from the adjoint matrix
    q0 = K[1][1]*K[2][2]*K[3][3] - K[1][1]*K[2][3]*K[3][2] - K[1][2]*K[2][1]*K[3][3] + K[1][2]*K[2][3]*K[3][1] + K[1][3]*K[2][1]*K[3][2] - K[1][3]*K[2][2]*K[3][1];
    q1 = K[1][0]*K[2][3]*K[3][2] - K[1][0]*K[2][2]*K[3][3] + K[1][2]*K[2][0]*K[3][3] - K[1][2]*K[2][3]*K[3][0] - K[1][3]*K[2][0]*K[3][2] + K[1][3]*K[2][2]*K[3][0];
    q2 = K[1][0]*K[2][1]*K[3][3] - K[1][0]*K[2][3]*K[3][1] - K[1][1]*K[2][0]*K[3][3] + K[1][1]*K[2][3]*K[3][0] + K[1][3]*K[2][0]*K[3][1] - K[1][3]*K[2][1]*K[3][0];
    q3 = K[1][0]*K[2][2]*K[3][1] - K[1][0]*K[2][1]*K[3][2] + K[1][1]*K[2][0]*K[3][2] - K[1][1]*K[2][2]*K[3][0] - K[1][2]*K[2][0]*K[3][1] + K[1][2]*K[2][1]*K[3][0];
    mstreal qsqr = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
    if (qsqr < evectol) {
      q0 = K[0][1]*K[2][3]*K[3][2] - K[0][1]*K[2][2]*K[3][3] + K[0][2]*K[2][1]*K[3][3] - K[0][2]*K[2][3]*K[3][1] - K[0][3]*K[2][1]*K[3][2] + K[0][3]*K[2][2]*K[3][1];
      q1 = K[0][0]*K[2][2]*K[3][3] - K[0][0]*K[2][3]*K[3][2] - K[0][2]*K[2][0]*K[3][3] + K[0][2]*K[2][3]*K[3][0] + K[0][3]*K[2][0]*K[3][2] - K[0][3]*K[2][2]*K[3][0];
      q2 = K[0][0]*K[2][3]*K[3][1] - K[0][0]*K[2][1]*K[3][3] + K[0][1]*K[2][0]*K[3][3] - K[0][1]*K[2][3]*K[3][0] - K[0][3]*K[2][0]*K[3][1] + K[0][3]*K[2][1]*K[3][0];
      q3 = K[0][0]*K[2][1]*K[3][2] - K[0][0]*K[2][2]*K[3][1] - K[0][1]*K[2][0]*K[3][2] + K[0][1]*K[2][2]*K[3][0] + K[0][2]*K[2][0]*K[3][1] - K[0][2]*K[2][1]*K[3][0];
      qsqr = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
      if (qsqr < evectol) {
        q0 = K[0][1]*K[1][2]*K[3][3] - K[0][1]*K[1][3]*K[3][2] - K[0][2]*K[1][1]*K[3][3] + K[0][2]*K[1][3]*K[3][1] + K[0][3]*K[1][1]*K[3][2] - K[0][3]*K[1][2]*K[3][1];
        q1 = K[0][0]*K[1][3]*K[3][2] - K[0][0]*K[1][2]*K[3][3] + K[0][2]*K[1][0]*K[3][3] - K[0][2]*K[1][3]*K[3][0] - K[0][3]*K[1][0]*K[3][2] + K[0][3]*K[1][2]*K[3][0];
        q2 = K[0][0]*K[1][1]*K[3][3] - K[0][0]*K[1][3]*K[3][1] - K[0][1]*K[1][0]*K[3][3] + K[0][1]*K[1][3]*K[3][0] + K[0][3]*K[1][0]*K[3][1] - K[0][3]*K[1][1]*K[3][0];
        q3 = K[0][0]*K[1][2]*K[3][1] - K[0][0]*K[1][1]*K[3][2] + K[0][1]*K[1][0]*K[3][2] - K[0][1]*K[1][2]*K[3][0] - K[0][2]*K[1][0]*K[3][1] + K[0][2]*K[1][1]*K[3][0];
        qsqr = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
        if (qsqr < evectol) {
          q0 = K[0][1]*K[1][3]*K[2][2] - K[0][1]*K[1][2]*K[2][3] + K[0][2]*K[1][1]*K[2][3] - K[0][2]*K[1][3]*K[2][1] - K[0][3]*K[1][1]*K[2][2] + K[0][3]*K[1][2]*K[2][1];
          q1 = K[0][0]*K[1][2]*K[2][3] - K[0][0]*K[1][3]*K[2][2] - K[0][2]*K[1][0]*K[2][3] + K[0][2]*K[1][3]*K[2][0] + K[0][3]*K[1][0]*K[2][2] - K[0][3]*K[1][2]*K[2][0];
          q2 = K[0][0]*K[1][3]*K[2][1] - K[0][0]*K[1][1]*K[2][3] + K[0][1]*K[1][0]*K[2][3] - K[0][1]*K[1][3]*K[2][0] - K[0][3]*K[1][0]*K[2][1] + K[0][3]*K[1][1]*K[2][0];
          q3 = K[0][0]*K[1][1]*K[2][2] - K[0][0]*K[1][2]*K[2][1] - K[0][1]*K[1][0]*K[2][2] + K[0][1]*K[1][2]*K[2][0] + K[0][2]*K[1][0]*K[2][1] - K[0][2]*K[1][1]*K[2][0];
          qsqr = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
        }
      }
    }

    // this is from http://onlinelibrary.wiley.com/doi/10.1002/jcc.20110/abstract (eq. 33)
    if (qsqr < evectol) {
      // if eigenvector still too small, then the optimal rotation is likely just the identity matrix
      u[0][0] = 1; u[0][1] = 0; u[0][2] = 0;
      u[1][0] = 0; u[1][1] = 1; u[1][2] = 0;
      u[2][0] = 0; u[2][1] = 0; u[2][2] = 1;
    } else {
      mstreal normq = sqrt(qsqr);
      q0 /= normq; q1 /= normq; q2 /= normq; q3 /= normq;
      u[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
      u[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
      u[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
      u[0][1] = 2*(q1*q2 - q0*q3); u[0][2] = 2*(q1*q3 + q0*q2);
      u[1][0] = 2*(q1*q2 + q0*q3); u[2][0] = 2*(q1*q3 - q0*q2);
      u[1][2] = 2*(q2*q3 - q0*q1);
      u[2][1] = 2*(q2*q3 + q0*q1);
    }

    /* To superimpose A onto B, we would:
     * 1. translate A by -cA to bring it to the origin
     * 2. apply the optimal rotatoin matrix
     * 3. translate the resulting coordinates by +cB to superimpose with B
     * So, if the original 3xN coordinate matrix is MA, then the operqtion is:
     * rot*(MA - cA) + cB = rot*MA + (cB - rot*cA) = rot*MA + trans
     * where trans = cB - rot*cA */
    for (i = 0; i < 3; i++) {
      t[i] = cB[i] - (u[i][0]*cA[0] + u[i][1]*cA[1] + u[i][2]*cA[2]);
    }

    // compute residuals after superposition
    if (setResiduals) {
      residuals.resize(N, vector<mstreal>(3, 0.0));
      mstreal xB, yB, zB;
      for (i = 0; i < N; i++) {
        xB = B[i]->getX() - cB[0];
        yB = B[i]->getY() - cB[1];
        zB = B[i]->getZ() - cB[2];
        residuals[i][0] = (A[i]->getX() - cA[0]) - (u[0][0]*xB + u[0][1]*yB + u[0][2]*zB);
        residuals[i][1] = (A[i]->getY() - cA[1]) - (u[1][0]*xB + u[1][1]*yB + u[1][2]*zB);
        residuals[i][2] = (A[i]->getZ() - cA[2]) - (u[2][0]*xB + u[2][1]*yB + u[2][2]*zB);
      }
    }
  }

   // when RMSD is almost exactly zero, the argument of the sqrt can end up very
   // small, but negative; so apply fabs to avoid nan answers
  return sqrt(fabs(GA + GB - 2*L)/N);
}

template <class T>
mstreal RMSDCalculator::qcpRMSDGrad(const T& A, const T& B, vector<mstreal>& grad) {
  mstreal rmsd = qcpRMSD(A, B, true, true); // set optimal residuals
  // if RMSD is zero, the gradient should be zero also (RMSD = 0 is the global minimum)
  mstreal f = (MstUtils::closeEnough(rmsd, 0.0) ? 0.0 : 1.0/(A.size() * rmsd));
  int k = 0;
  for (int i = 0; i < residuals.size(); i++) {
    for (int j = 0; j < 3; j++) {
      grad[k] = f * residuals[i][j];
      k++;
    }
  }
  return rmsd;
}

bool RMSDCalculator::testQCP(bool testGrad) {
  long int x = std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
  srand(x);
  RMSDCalculator rc;

  // test a bunch of times
  long int Tkabsch = 0, Tqcp = 0;
  bool failed = false; int N = 100000;
  for (int k = 0; k < N; k++) {
    int N = MstUtils::randInt(100, 10); // number of atoms
    mstreal L = (1 + MstUtils::randUnit())*50; // length scale

    // create random atoms
    AtomPointerVector A(N), B(N);
    for (int i = 0; i < N; i++) {
      A[i] = new Atom(0, "X", MstUtils::randUnit()*L, MstUtils::randUnit()*L, MstUtils::randUnit()*L, 0, 0, false);
      B[i] = new Atom(0, "X", MstUtils::randUnit()*L, MstUtils::randUnit()*L, MstUtils::randUnit()*L, 0, 0, false);
    }

    // superimpose different ways
    bool ord = (MstUtils::randUnit() > 0.5);
    mstreal rmsd, rmsdQCP;
    for (int i = 0; i < 2; i++) {
      if (ord == (bool) i) {
        Tkabsch -= std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
        rmsd = rc.bestRMSD(A, B);
        Tkabsch += std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
      } else {
        Tqcp -= std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
        rmsdQCP = rc.qcpRMSD(A, B, true);
        Tqcp += std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
      }
    }

    // test
    if (fabs(rmsd - rmsdQCP) > 10E-8) {
      failed = true;
      cout << "test FAILED for RMSD calculation, point sets:\n" << A << endl << B << endl;
      cout << "Kabsch: " << rmsd << endl;
      cout << "QCP: " << rmsdQCP << endl;
      cout << "difference: " << (rmsd - rmsdQCP) << endl;
    }

    // test gradient calculation
    if (testGrad) {
      CartesianPoint rmsdGrad(3*A.size(), 0), rmsdGradNum(3*A.size(), 0);
      mstreal del = 0.0001*L;
      rc.qcpRMSDGrad(A, B, rmsdGrad);
      int k = 0;
      for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < 3; j++) {
          mstreal old = (*A[i])[j];
          (*A[i])[j] = old + del; mstreal p = rc.bestRMSD(A, B);
          (*A[i])[j] = old - del; mstreal m = rc.bestRMSD(A, B);
          rmsdGradNum[k] = (p - m)/(2*del);
          (*A[i])[j] = old;
          k++;
        }
      }
      if ((rmsdGrad - rmsdGradNum).norm() > 10E-3) {
        failed = true;
        cout << "test FAILED for RMSD gradient calculation, point sets:\n" << A << endl << B << endl;
        cout << "QCP\tnumeric\tdifference" << endl;
        for (int i = 0; i < rmsdGrad.size(); i++) cout << rmsdGrad[i] << "\t" << rmsdGradNum[i] << "\t" << rmsdGrad[i] - rmsdGradNum[i] << endl;
        cout << "difference norm: " << (rmsdGrad - rmsdGradNum).norm() << endl;
      }
    }

    A.deletePointers();
    B.deletePointers();
    if (failed) return false;
  }
  cout << "Kabsch: " << (10E9/(Tkabsch/N)) << " per second" << endl;
  cout << "Qcp:    " << (10E9/(Tqcp/N)) << " per second" << endl;
  cout << "Kabsch/Qcp = " << Tqcp*1.0/Tkabsch << endl;
  return true;
}

// forward declarations of template functions
template mstreal RMSDCalculator::qcpRMSD<vector<Atom*> >(const vector<Atom*>& A, const vector<Atom*>& B, bool setTransform, bool setResiduals);
template mstreal RMSDCalculator::qcpRMSDGrad<vector<Atom*> >(const vector<Atom*>& A, const vector<Atom*>& B, vector<mstreal>& grad);
template mstreal RMSDCalculator::qcpRMSD<AtomPointerVector>(const AtomPointerVector& A, const AtomPointerVector& B, bool setTransform, bool setResiduals);
template mstreal RMSDCalculator::qcpRMSDGrad<AtomPointerVector>(const AtomPointerVector& A, const AtomPointerVector& B, vector<mstreal>& grad);
template mstreal RMSDCalculator::qcpRMSD<vector<CartesianPoint*> >(const vector<CartesianPoint*>& A, const vector<CartesianPoint*>& B, bool setTransform, bool setResiduals);
template mstreal RMSDCalculator::qcpRMSDGrad<vector<CartesianPoint*> >(const vector<CartesianPoint*>& A, const vector<CartesianPoint*>& B, vector<mstreal>& grad);


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

ProximitySearch::ProximitySearch(const vector<CartesianPoint>& _points, mstreal _characteristicDistance, bool _addPoints, vector<int>* tags, mstreal pad) {
  calculateExtent(_points);
  if (xlo == xhi) { xlo -= _characteristicDistance/2; xhi += _characteristicDistance/2; }
  if (ylo == yhi) { ylo -= _characteristicDistance/2; yhi += _characteristicDistance/2; }
  if (zlo == zhi) { zlo -= _characteristicDistance/2; zhi += _characteristicDistance/2; }
  xlo -= pad; ylo -= pad; zlo -= pad;
  xhi += pad; yhi += pad; zhi += pad;
  int _N = int(ceil(max(max((xhi - xlo), (yhi - ylo)), (zhi - zlo))/_characteristicDistance));
  reinitBuckets(_N);
  setBinWidths();
  if (_addPoints) {
    for (int i = 0; i < _points.size(); i++) {
      addPoint(_points[i], (tags == NULL) ? i : (*tags)[i]);
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
  pointTags = ps.pointTags;
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

void ProximitySearch::calculateExtent(const vector<CartesianPoint>& _points, mstreal& _xlo, mstreal& _ylo, mstreal& _zlo, mstreal& _xhi, mstreal& _yhi, mstreal& _zhi) {
  if (_points.size() == 0) { cout << "Error in ProximitySearch::calculateExtent() -- empty cartesian point vector passed!\n"; exit(-1); }
  _xlo = _xhi = _points[0].getX();
  _ylo = _yhi = _points[0].getY();
  _zlo = _zhi = _points[0].getZ();
  for (int i = 0; i < _points.size(); i++) {
    if (_xlo > _points[i].getX()) _xlo = _points[i].getX();
    if (_ylo > _points[i].getY()) _ylo = _points[i].getY();
    if (_zlo > _points[i].getZ()) _zlo = _points[i].getZ();
    if (_xhi < _points[i].getX()) _xhi = _points[i].getX();
    if (_yhi < _points[i].getY()) _yhi = _points[i].getY();
    if (_zhi < _points[i].getZ()) _zhi = _points[i].getZ();
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
  *i = int((px - xlo)/xbw + 0.5);
  *j = int((py - ylo)/ybw + 0.5);
  *k = int((pz - zlo)/zbw + 0.5);
}

void ProximitySearch::limitIndex(int *ind) {
  if (*ind < 0) *ind = 0;
  if (*ind > N-1) *ind = N-1;
}

mstreal ProximitySearch::limitX(mstreal x) {
  if (x < xlo) return xlo;
  if (x > xhi) return xhi;
  return x;
}

mstreal ProximitySearch::limitY(mstreal y) {
  if (y < ylo) return ylo;
  if (y > yhi) return yhi;
  return y;
}

mstreal ProximitySearch::limitZ(mstreal z) {
  if (z < zlo) return zlo;
  if (z > zhi) return zhi;
  return z;
}

bool ProximitySearch::pointsWithin(const CartesianPoint& c, mstreal dmin, mstreal dmax, vector<int>* list, bool byTag) {
  mstreal cx = c.getX(); mstreal cy = c.getY(); mstreal cz = c.getZ();
  // first check if the point is outside of the bounding box of the point cloud by a sufficient amount
  if ((cx < xlo - dmax) || (cy < ylo - dmax) || (cz < zlo - dmax) || (cx > xhi + dmax) || (cy > yhi + dmax) || (cz > zhi + dmax)) return false;

  mstreal d2, dmin2, dmax2;
  int ci, cj, ck;
  int iOutLo, jOutLo, kOutLo, iOutHi, jOutHi, kOutHi; // external box (no point in looking beyond it, points there are too far)
  int iInLo, jInLo, kInLo, iInHi, jInHi, kInHi;       // internal box (no point in looking within it, points there are too close)
  pointBucket(limitX(cx), limitY(cy), limitZ(cz), &ci, &cj, &ck);
  pointBucket(limitX(cx - dmax), limitY(cy - dmax), limitZ(cz - dmax), &iOutLo, &jOutLo, &kOutLo);
  pointBucket(limitX(cx + dmax), limitY(cy + dmax), limitZ(cz + dmax), &iOutHi, &jOutHi, &kOutHi);
  if (dmin > 0) {
    mstreal sr3 = sqrt(3);
    pointBucket(cx - dmin/sr3, cy - dmin/sr3, cz - dmin/sr3, &iInLo, &jInLo, &kInLo);
    pointBucket(cx + dmin/sr3, cy + dmin/sr3, cz + dmin/sr3, &iInHi, &jInHi, &kInHi);
  } else {
    iInLo = iInHi = ci;
    jInLo = jInHi = cj;
    kInLo = kInHi = ck;
  }
  // limitIndex(&iInLo); limitIndex(&iInHi); limitIndex(&jInLo); limitIndex(&jInHi); limitIndex(&kInLo); limitIndex(&kInHi);
  // limitIndex(&iOutLo); limitIndex(&iOutHi); limitIndex(&jOutLo); limitIndex(&jOutHi); limitIndex(&kOutLo); limitIndex(&kOutHi);

  // search only within the boxes where points of interest can be, in principle
  if (list != NULL) list->clear();
  bool found = false, insi, ins;
  bool yesno = (list == NULL);
  dmin2 = dmin*dmin; dmax2 = dmax*dmax;
  int i, j, k, ii, pi;
  for (i = iOutLo; i <= iOutHi; i++) {
    insi = (i > iInLo) && (i < iInHi);
    vector<vector<vector<int> > >& Bi = buckets[i];
    for (j = jOutLo; j <= jOutHi; j++) {
      ins = insi && (j > jInLo) && (j < jInHi);
      vector<vector<int> >& Bij = Bi[j];
      for (k = kOutLo; k <= kOutHi; k++) {
        vector<int>& Bijk = Bij[k];
        // check all points in bucket i, j, k
        for (ii = 0; ii < Bijk.size(); ii++) {
          pi = Bijk[ii];
          d2 = c.distance2nc(pointList[pi]);
          if ((d2 >= dmin2) && (d2 <= dmax2)) {
            if (yesno) return true;
            list->push_back(byTag ? pointTags[pi] : pi);
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
vector<vector<int>> Clusterer::greedyCluster(const vector<vector<Atom*>>& units, mstreal rmsdCut, int Nmax, mstreal coverage, int maxClusters, bool verbose) {
  vector<vector<int>> clusters;
  if (units.empty()) return clusters;
  set<int> remIndices;
  for (int i = 0; i < units.size(); i++) remIndices.insert(i);
  if (remIndices.size() <= Nmax) return Clusterer::greedyClusterBruteForce(units, remIndices, rmsdCut);
  
  cout << "There are " << units.size() << " total points to cluster, will continue until " << floor(units.size()*coverage) << " (" << coverage << ") or " << maxClusters << " are covered" << endl;

  // create some dummy storage vectors
  int L = units[0].size();
  AtomPointerVector mean(L, NULL), copy(L, NULL);
  for (int i = 0; i < L; i++) { mean[i] = new Atom(); copy[i] = new Atom(); }

  // determine how many points can be left over (to be distributed to already formed clusters) given the specified coverage fraction
  int total = remIndices.size();
  int numToLeave = total * (1 - coverage);

  // iterate to find a new cluster each time
  while ((remIndices.size() > Nmax) && (remIndices.size() > numToLeave) && ((maxClusters <= 0) || ((maxClusters > 0) && (clusters.size() < maxClusters)))) {
    // sub-sample Nmax elements
    set<int> subSample = Clusterer::randomSubsample(remIndices, Nmax);

    // get the top cluster from these and use its centroid
    vector<int> topClustSub = Clusterer::greedyClusterBruteForce(units, subSample, rmsdCut, 1)[0];
    vector<int> topClust = Clusterer::elementsWithin(units, remIndices, units[topClustSub[0]], rmsdCut);
    if (verbose) cout << "picked initial cluster with " << topClust.size() << " points..." << endl;

    // now try to improve the centroid by moving it closer to the average
    while (1) {
      if (verbose) cout << "\timproving..." << endl;
      mean.copyCoordinates(units[topClust[0]]);
      for (int i = 1; i < topClust.size(); i++) {
        copy.copyCoordinates(units[topClust[i]]);
        if (optimAlign) rCalc.align(copy, units[topClust[0]], copy);
        mean *= (mstreal) i; mean += copy; mean /= (mstreal) (i + 1);
      }
      vector<int> topClustNew = Clusterer::elementsWithin(units, remIndices, mean, rmsdCut);
      topClustNew = Clusterer::elementsWithin(units, remIndices, units[topClustNew[0]], rmsdCut);
      if (topClustNew.size() <= topClust.size()) break;
      topClust = topClustNew;
      if (verbose) cout << "\timproved to " << topClust.size() << " points" << endl;
    }

    // // now try to grow the cluster by improving the centroid, if there are compute cycles left
    // int Ntry = int(round(Nmax * Nmax * 1.0 / remIndices.size()));
    // for (int i = 1; i <= MstUtils::min(Ntry, (int) topClust.size()-1); i++) {
    //   vector<int> newTopClust = Clusterer::elementsWithin(units, remIndices, units[topClust[i]], rmsdCut);
    //   if (newTopClust.size() > topClust.size()) topClust = newTopClust;
    // }

    // keep whatever cluster end up with, exclude its elements
    clusters.push_back(topClust);
    for (int i = 0; i < topClust.size(); i++) remIndices.erase(topClust[i]);
    if (verbose) cout << remIndices.size() << " points remaining" << endl;
  }

  // brute force through the rest
  if ((maxClusters > 0) && (clusters.size() < maxClusters)) {
    int remainingClusters = maxClusters - clusters.size();
    vector<vector<int>> remClusters = greedyClusterBruteForce(units, remIndices, rmsdCut, remainingClusters);
    clusters.insert(clusters.end(), remClusters.begin(), remClusters.end());
  } else if ((maxClusters <= 0) && (remIndices.size() > numToLeave)) {
    set<int> indsToCluster;
    auto curInd = remIndices.begin();
    for (int i = 0; i < (remIndices.size() - numToLeave); i++) {
      auto it = remIndices.begin();
      indsToCluster.insert(*it);
      remIndices.erase(it);
    }
    vector<vector<int>> remClusters = greedyClusterBruteForce(units, indsToCluster, rmsdCut);
    clusters.insert(clusters.end(), remClusters.begin(), remClusters.end());
  }

  // clean up
  mean.deletePointers();
  copy.deletePointers();

  return clusters;
}

vector<vector<int> > Clusterer::greedyClusterBruteForce(const vector<vector<Atom*> >& units, set<int> remIndices, mstreal rmsdCut, int nClusts) {
  vector<vector<int> > clusters;
  while (remIndices.size() != 0) {
    // pick the best current centroid
    vector<int> bestClust;
    for (auto it = remIndices.begin(); it != remIndices.end(); ++it) {
      vector<int> clust = elementsWithin(units, remIndices, units[*it], rmsdCut);
      if (clust.size() > bestClust.size()) {
        bestClust = clust;
      }
    }
    // add its corresponding cluster
    clusters.push_back(bestClust);
    for (int i = 0; i < bestClust.size(); i++) remIndices.erase(bestClust[i]);

    // where we asked for at most the top some number of clusters?
    if ((nClusts > 0) && (clusters.size() >= nClusts)) break;
  }
  return clusters;
}

vector<int> Clusterer::elementsWithin(const vector<vector<Atom*> >& units, set<int>& remIndices, const vector<Atom*>& fromUnit, mstreal rmsdCut) {
  vector<int> neigh; vector<mstreal> rmsds;
  for (auto it = remIndices.begin(); it != remIndices.end(); ++it) {
    mstreal r = optimAlign ? rCalc.bestRMSD(fromUnit, units[*it]) : rCalc.rmsd(fromUnit, units[*it]);
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

void MstUtils::fileToArray(const string& _filename, vector<string>& lines) {
  fstream inp;
  MstUtils::openFile(inp, _filename, ios_base::in, "MstUtils::fileToArray");
  string line;
  while (true) {
    getline(inp, line);
    // if eof set upon trying to read the line, and the line is empty, then it
    // was not really a valid line (but just the last eof that was not yet read)
      if (!inp.eof() || !line.empty()) lines.push_back(line);
    if (inp.eof()) break;
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
    ret = str.substr(0, i);
    str = str.substr(min((size_t) i+1, str.length()-1));
    return ret;
  }
  ret = str.substr(0, 1);
  str = str.substr(min((size_t) 1, str.length()));
  return ret;
}

string MstUtils::nextQuoteAwareToken(string& str, string quoteMarks, string delimiters, bool skipTrailingDelims) {
  string ret;
  if (!delimiters.empty()) {
    if (skipTrailingDelims) str = trim(str, delimiters);
    size_t firstQuote = string::npos;
    for (size_t i = 0; i < str.size(); i++) {
      for (int j = 0; j < quoteMarks.size(); j++) {
        if ((str[i] == quoteMarks[j]) && ((i == 0) || (str[i-1] != '\\'))) { firstQuote = i; break; }
      }
      if (firstQuote != string::npos) break;
    }
    size_t firstDelim = string::npos;
    for (int i = 0; i < delimiters.size(); i++) firstDelim = min(firstDelim, str.find_first_of(delimiters[i]));
    if (firstQuote < firstDelim) {
      char quoteChar = str[firstQuote];
      size_t nextQuote = string::npos;
      for (size_t i = firstQuote+1; i < str.size(); i++) {
        if ((str[i] == quoteChar) && (str[i-1] != '\\')) { nextQuote = i; break; }
      }
      if (nextQuote == string::npos) {
        MstUtils::error("malformed quotation (missing the trailing quote mark " + string(1, quoteChar) + ")", "MstUtils::nextQuoteAwareToken");
      }
      string beforeQuote = str.substr(0, firstQuote);
      string insideQuote = str.substr(firstQuote+1, nextQuote-firstQuote-1);
      string trailingChars = str.substr(nextQuote+1);
      string afterQuote = MstUtils::nextToken(trailingChars, delimiters, false);
      ret = beforeQuote + insideQuote + afterQuote;
      str = trailingChars;
    } else {
      string beforeQuote = str.substr(0, firstQuote);
      string quoteStart = str.substr(min(firstQuote, str.size()));
      ret = MstUtils::nextToken(beforeQuote, delimiters, false);
      str = beforeQuote + quoteStart;
    }
    ret = unescape(ret);
  } else {
    if (!str.empty()) {
      bool foundDelim = false;
      for (int i = 0; i < quoteMarks.size(); i++) {
        if (str[0] == quoteMarks[i]) {
          char quoteChar = quoteMarks[i];
          int j = 1;
          while ((j < str.size()) && !((str[j] == quoteChar) && (str[j-1] != '\\'))) j++;
          if (j == str.size()) {
            MstUtils::error("malformed quotation (missing the trailing quote mark " + string(1, quoteChar) + ")", "MstUtils::nextQuoteAwareToken");
          }
          ret = str.substr(1, j-1);
          str = str.substr(j+1);
          foundDelim = true;
          break;
        }
      }
      if (!foundDelim) {
        if (str[0] == '\\') str = unescape(str.substr(0, 2)) + str.substr(2);
        ret = str.substr(0, 1);
        str = str.substr(1);
      }
    } else ret = "";
  }
  return ret;
}

string MstUtils::escape(string str, string special) {
  for (int i = 0; i < str.size()-1; i++) {
    for (int j = 0; j < special.size(); j++) {
      if (str[i] == special[j]) {
        str.replace(i, 1, "\\" + string(1, str[i]));
        i++; break;
      }
    }
  }
  if (str[str.size()-1] == '\\') {
    MstUtils::error("malformed escape sequence (any \\ must be followed by another character, with another \\ producing a literal backslash)", "MstUtils::escape");
  }
  return str;
}

string MstUtils::unescape(string str) {
  string unescapedStr = str;
  int i = 0;
  while (i < unescapedStr.size()) {
    if (unescapedStr[i] == '\\') {
      if (i+1 < unescapedStr.size()) unescapedStr = unescapedStr.substr(0, i) + unescapedStr.substr(i+1);
      else MstUtils::error("malformed escape sequence (any \\ must be followed by another character, with another \\ producing a literal backslash)", "MstUtils::unescape");
    }
    i++;
  }
  return unescapedStr;
}

vector<string> MstUtils::split(const string& str, string delimiters, bool skipTrailingDelims) {
  string strCopy = str;
  vector<string> tokens;
  while (!strCopy.empty()) tokens.push_back(nextToken(strCopy, delimiters, skipTrailingDelims));
  return tokens;
}

string MstUtils::removeComment(const string& str, string commStart) {
  for (int i = 0; i < str.length() - commStart.length() + 1; i++) {
    if (str.substr(i, commStart.length()).compare(commStart) == 0) {
      return str.substr(0, i);
    }
  }
  return str;
}

vector<MST::mstreal> MstUtils::splitToReal(const string& str, string delimiters, bool skipTrailingDelims, bool strict) {
  vector<string> tokens = MstUtils::split(str, delimiters, skipTrailingDelims);
  vector<MST::mstreal> vec(tokens.size());
  for (int i = 0; i < tokens.size(); i++) vec[i] = MstUtils::toReal(tokens[i], strict);
  return vec;
}

vector<int> MstUtils::splitToInt(const string& str, string delimiters, bool skipTrailingDelims, bool strict) {
  vector<string> tokens = MstUtils::split(str, delimiters, skipTrailingDelims);
  vector<int> vec(tokens.size());
  for (int i = 0; i < tokens.size(); i++) vec[i] = MstUtils::toInt(tokens[i], strict);
  return vec;
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

string MstUtils::trim(const string& str, string delimiters) {
  int i = str.find_first_not_of(delimiters);
  if (i == string::npos) return "";
  int j = str.find_last_not_of(delimiters);
  return str.substr(i, j - i + 1);
}
vector<string> MstUtils::trim(const vector<string>& strings, string delimiters) {
  vector<string> ret = strings;
  for (int i = 0; i < strings.size(); i++) ret[i] = MstUtils::trim(strings[i], delimiters);
  return ret;
}


void MstUtils::warn(const string& message, string from) {
  string head = from.empty() ? "Warning: " : "Warning in " + from + ": ";
  cerr << head << wrapText(message, 100, 0, head.length()) << endl;
}

void MstUtils::error(const string& message, string from, int code) {
  string head = from.empty() ? "Error: " : "Error in " + from + ": ";
  cerr << head << wrapText(message, 100, 0, head.length()) << endl;

  // print backtrace
  void *array[100];
  size_t size = backtrace(array, 100);
  backtrace_symbols_fd(array, size, STDERR_FILENO);

  throw(code);
  // exit(code);
}

MST::mstreal MstUtils::randNormal(MST::mstreal mu, MST::mstreal sig) {
  normal_distribution<MST::mstreal> dist(mu, sig);
  return dist(MstUtils::randEngine());
}

MST::mstreal MstUtils::randUnit(MST::mstreal mi, MST::mstreal ma) {
  uniform_real_distribution<MST::mstreal> dist(mi, ma);
  return dist(MstUtils::randEngine());
}

int MstUtils::randInt(int lower, int upper) {
  uniform_int_distribution<int> dist(lower, upper);
  return dist(MstUtils::randEngine());
}

void MstUtils::errorHandler(int sig) {
  fprintf(stderr, "Error: signal %d:\n", sig);
  MstUtils::error("printing trace:");
}

void MstUtils::setSignalHandlers() {
  signal(SIGABRT, MstUtils::errorHandler);
  signal(SIGFPE, MstUtils::errorHandler);
  signal(SIGILL, MstUtils::errorHandler);
  signal(SIGINT, MstUtils::errorHandler);
  signal(SIGSEGV, MstUtils::errorHandler);
  signal(SIGTERM, MstUtils::errorHandler);
}

void MstUtils::assertCond(bool condition, string message, string from, int exitCode) {
  if(!condition) {
    MstUtils::error(message, from, exitCode);
  }
}

string MstUtils::wrapText(const string& message, int width, int leftSkip, int startingOffset) {
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

int MstUtils::toInt(const string& num, bool strict) {
  int ret = 0;
  try { ret = stoi(num); }
  catch (...) {
    if (strict) MstUtils::error("failed to convert '" + num + "' to integer", "MstUtils::toInt");
  }
  return ret;
  // int ret = 0;
  // if ((sscanf(num.c_str(), "%d", &ret) != 1) && strict) MstUtils::error("failed to convert '" + num + "' to integer", "MstUtils::toInt");
  // return ret;
}

bool MstUtils::isInt(const string& num) {
  string::size_type k;
  try {
    stoi(num, &k);
  }
  catch (...) { // no matter what exception, fail
    return false;
  }
  bool period = false;
  for (int i = k; i < num.size(); i++) {
    if ((num[i] == '.') && !period) { period = true; continue; } // one period is okay (e.g., "12." is an int)
    if (!isspace(num[i])) return false;
  }
  return true;
}

bool MstUtils::isReal(const string& num) {
  double ret;
  return (sscanf(num.c_str(), "%lf", &ret) == 1);
}

MST::mstreal MstUtils::toReal(const string& num, bool strict) {
  double ret = 0.0;
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

vector<pair<int, int> > MstUtils::splitTasks(int numTasks, int numJobs) {
  int numPerJob = numTasks / numJobs;
  int remTasks = numTasks % numJobs;
  vector<pair<int, int> > division(numJobs);
  int k = 0;
  for (int i = 0; i < numJobs; i++) {
    division[i].first = k;
    division[i].second = k + numPerJob - 1;
    if (i < remTasks) division[i].second++;
    k = division[i].second + 1;
  }
  if (k != numTasks) MstUtils::error("something went very wrong!", "MstUtils::splitTasks(int, int)");
  return division;
}
