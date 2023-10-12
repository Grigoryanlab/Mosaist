#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <typeinfo>
#include <stdio.h>
#include <map>
#include <iomanip>
#include <unistd.h>
#include <ctime>

#include "msttypes.h"
#include "mstrotlib.h"
#include "mstcondeg.h"
#include "mstoptions.h"
#include "mstsystem.h"

using namespace std;
using namespace MST;

// ---- utility functions
void proteinOnly(System& S, System& So, vector<string>& legalNames);

class options {
  public:
    options() {
      dcut = 25.0;
      clashDist = 2.0; contDist = 3.0; rotLibFile = ""; beblFile = ""; rotOutFile = ""; rotLevel = "";
      verbose = renumPDB = phi_psi = omega = printFileNames = seq_const = false;
      aaProp["ALA"] = 7.73; aaProp["CYS"] = 1.84; aaProp["ASP"] = 5.82; aaProp["GLU"] = 6.61; aaProp["PHE"] = 4.05;
      aaProp["GLY"] = 7.11; aaProp["HIS"] = 2.35; aaProp["HSD"] = 2.35; aaProp["ILE"] = 5.66; aaProp["LYS"] = 6.27;
      aaProp["LEU"] = 8.83; aaProp["MET"] = 2.08; aaProp["ASN"] = 4.50; aaProp["PRO"] = 4.52; aaProp["GLN"] = 3.94;
      aaProp["ARG"] = 5.03; aaProp["SER"] = 6.13; aaProp["THR"] = 5.53; aaProp["VAL"] = 6.91; aaProp["TRP"] = 1.51; aaProp["TYR"] = 3.54;
    }
    vector<string> pdbfs, omapfs, opdbfs;
    bool verbose, renumPDB, phi_psi, printFileNames, omega, freeB, seq_const;
    string focus, rotLibFile, beblFile, rotOutFile, rotLevel;
    double dcut, clashDist, contDist;
    map<string, double> aaProp; // amino-acid propensities (in percent)
};

string option(string opt, string mes, int w, int p1, int p2) {
  // first print the name of the option
  string text(p1, ' ');
  text += opt;
  if (p2 > text.size()) text += string(p2 - text.size(), ' ');

  // next print the description text
  int i = 0, k, L = text.size(), n;
  while (i < mes.size()) {
    k = mes.find_first_of(" ", i);
    if (k == string::npos) k = mes.size();
    n = k - i;
    if ((L + n >= w) && (L > 0)) { text += "\n" + string(p2, ' '); L = p2; }
    text += mes.substr(i, n) + " ";
    L += n + 1;
    i = k+1;
  }
  return text;
}


// ---- Main Program
int main(int argc, char *argv[]) {
  fstream rotOut, *rotOutPtr = NULL;

  // process input arguments
  options iopts;
  MstOptions op;
  op.setTitle("Identifies inter-positional contacts and environment information from input PDB file(s). Options:");
  op.addOption("p", "input PDB file.");
  op.addOption("pL", "a file with a list of PDB files. Either --p or --pL must be given.");
  op.addOption("rLib", "a path to an MSL-formatter rotamer library, with WEIGHTS information. If a file, will be treated as a backbone-independent library. If a directory, will look for files EBL.out and BEBL.out in it to read a backbone-dependent library.", true);
  op.addOption("o", "optional: output file name for writing contacts. If not given, will write to standard output.");
  op.addOption("oL", "optional: a file with a list of contact file names (one per input PDB structure).");
  op.addOption("opdb", "optional: output post-processed PDB file (useful for keeping track of how all PDB weirdnesses got parsed).");
  op.addOption("opdbL", "optional: a file with a list of file names for post-processed PDBs, one per input structure.");
  op.addOption("freeB", "optional: if given and output of post-processed PDB file(s) was requested, will replace B-factors with freedom values scaled by 100 (-1 for residues not considered).");
  op.addOption("rout", "optional: name of a file into which to place PDB-formated coordinates of rotamers that ended up surviving at each considered position.");
  op.addOption("sel", "optional: a selection string for limiting the residues for which properties will be computed.");
  op.addOption("pp", "optional: print phi/psi angles for each residue (will print next to all positional scores score).");
  op.addOption("omg", "optional: print omega angles for each residue (will print next to all positional scores score). Omega for the current position is defined as the -CA, C, N, CA dihedral angle.");
  op.addOption("seq_const","optional: recompute the CD and INT using sequence constraints.");
  op.addOption("verb", "optional: generate lots of detailed output (i.e., print the details of which rotamer pairs are in contact).");
  op.addOption("pf", "if flag specified, will print the name of the PDB file being analyzed next to all positional scores. This is especially convenient when a list of PDB file is specified as input and the output goes to a single file.");
  op.addOption("ren", "if flag specified, will renumber the structure before output. Useful for keeping track of residues in the output list of contacts if the input PDB file is strangely numbered.");

  op.setOptions(argc, argv);
  
  MstTimer timer; timer.start();

  MstUtils::assertCond(op.isGiven("p") || op.isGiven("pL"), "either --p or --pL must be specified!");
  if (op.isGiven("p")) iopts.pdbfs.push_back(op.getString("p"));
  if (op.isGiven("pL")) MstUtils::fileToArray(op.getString("pL"), iopts.pdbfs);
  if (op.isGiven("o")) iopts.omapfs.push_back(op.getString("o"));
  if (op.isGiven("oL")) MstUtils::fileToArray(op.getString("oL"), iopts.omapfs);
  if (op.isGiven("opdb")) iopts.opdbfs.push_back(op.getString("opdb"));
  if (op.isGiven("opdbL")) iopts.opdbfs.push_back(op.getString("opdbL"));
  if (op.isGiven("sel")) iopts.focus = string(op.getString("sel"));
  if (op.isGiven("rLib")) iopts.rotLibFile = op.getString("rLib");
  if (op.isGiven("rout")) iopts.rotOutFile = op.getString("rout");
  iopts.verbose = op.isGiven("verb");
  iopts.phi_psi = op.isGiven("pp");
  iopts.omega = op.isGiven("omg");
  iopts.renumPDB = op.isGiven("ren");
  iopts.printFileNames = op.isGiven("pf");
  iopts.freeB = op.isGiven("freeB");
  iopts.seq_const = op.isGiven("seq_const");

  // error checking
  // make sure lists are of the proper size
  if (iopts.omapfs.size() > 1) MstUtils::assertCond(iopts.omapfs.size() == iopts.pdbfs.size(), "the number of input PDB files and output map files does not agree");
  if (iopts.opdbfs.size() > 1) MstUtils::assertCond(iopts.opdbfs.size() == iopts.pdbfs.size(), "the number of input PDB files and output PDB files does not agree");
  if (iopts.pdbfs.size() > 1) {
    if (iopts.omapfs.size() == 1) {
      iopts.omapfs.resize(iopts.pdbfs.size());
      string base = MstSys::pathBase(iopts.omapfs[0]);
      for (int i = 0; i < iopts.pdbfs.size(); i++) { iopts.omapfs[i] = base + ".f" + MstUtils::toString(i+1) + ".cont"; }
    }
    if (iopts.opdbfs.size() == 1) {
      iopts.opdbfs.resize(iopts.pdbfs.size());
      string base = MstSys::pathBase(iopts.opdbfs[0]);
      for (int i = 0; i < iopts.pdbfs.size(); i++) { iopts.opdbfs[i] = base + ".f" + MstUtils::toString(i+1) + ".pdb"; }
    }
  }

  // legal residue names that are considered "protein" here
  vector<string> legalNames;
  legalNames.push_back("ALA"); legalNames.push_back("CYS"); legalNames.push_back("ASP"); legalNames.push_back("GLU"); legalNames.push_back("PHE"); legalNames.push_back("GLY");
  legalNames.push_back("HIS"); legalNames.push_back("ILE"); legalNames.push_back("LYS"); legalNames.push_back("LEU"); legalNames.push_back("MET"); legalNames.push_back("ASN");
  legalNames.push_back("PRO"); legalNames.push_back("GLN"); legalNames.push_back("ARG"); legalNames.push_back("SER"); legalNames.push_back("THR"); legalNames.push_back("VAL");
  legalNames.push_back("TRP"); legalNames.push_back("TYR"); legalNames.push_back("HSD"); legalNames.push_back("HSE"); legalNames.push_back("HSC"); legalNames.push_back("HSP");
  legalNames.push_back("MSE");
  legalNames.push_back("CSO"); legalNames.push_back("HIP"); legalNames.push_back("SEC"); legalNames.push_back("SEP"); legalNames.push_back("TPO"); legalNames.push_back("PTR");

  // pre-read rotamer library once
  RotamerLibrary RL(iopts.rotLibFile);

  // go through all PDB files
  for (int si = 0; si < iopts.pdbfs.size(); si++) {
    Structure So(iopts.pdbfs[si]);                       // original input PDB structure
    Structure S = So;
    // Structure S; proteinOnly(S, So, legalNames);
    if (iopts.renumPDB) S.renumber();
    ConFind C(&RL, S);
    if (!iopts.rotOutFile.empty()) C.openLogFile(iopts.rotOutFile, si > 0);
    if (iopts.freeB) {
      AtomPointerVector atoms = S.getAtoms();
      for (int ai = 0; ai < atoms.size(); ai++) atoms[ai]->setB(-1);
    }

    // optionally select the relevant residue subset
    vector<Residue*> allRes;
    if (!iopts.focus.empty()) {
      selector sel(S);
      allRes = sel.selectRes(iopts.focus);
    } else {
      allRes = S.getResidues();
    }

    // open output file and write header
    fstream of;
    streambuf * buf;
    if (!iopts.omapfs.empty()) {
      MstUtils::openFile(of, iopts.omapfs[si], fstream::out);
      buf = of.rdbuf();
    } else {
      cout << iopts.pdbfs[si] << endl;
      buf = cout.rdbuf();
    }
    ostream out(buf);

    // print degrees
    contactList L;
    C.getContacts(allRes, 0, &L);
    vector<pair<Residue*, Residue*> > list = L.getOrderedContacts();
    for (int k = 0; k < list.size(); k++) {
      Residue* resA = list[k].first;
      Residue* resB = list[k].second;
      out << "contact\t" << resA->getChainID() << "," << resA->getNum() << "\t" << resB->getChainID() << "," << resB->getNum();
      out << "\t" << std::setprecision(6) << std::fixed << L.degree(resA, resB);
      out << "\t" << resA->getName() << "\t" << resB->getName();
      if (iopts.printFileNames) out << "\t" << iopts.pdbfs[si];
      out << endl;
    }

    // print crowdedness
    for (int k = 0; k < allRes.size(); k++) {
      Residue* res = allRes[k];
      out << "crwdnes\t" << res->getChainID() << "," << res->getNum() << "\t";
      out << std::setprecision(6) << std::fixed << C.getCrowdedness(res) << "\t";
      if (iopts.phi_psi) out << res->getPhi() << "\t" << res->getPsi() << "\t";
      if (iopts.omega) out << res->getOmega() << "\t";
      out << res->getName();
      if (iopts.printFileNames) out << "\t" << iopts.pdbfs[si];
      out << endl;
    }

    // print freedoms
    vector<mstreal> freedoms = C.getFreedom(allRes);
    for (int k = 0; k < allRes.size(); k++) {
      Residue* res = allRes[k];
      out << "freedom\t" << res->getChainID() << "," << res->getNum() << "\t";
      out << std::setprecision(6) << std::fixed << freedoms[k] << "\t";
      if (iopts.phi_psi) out << res->getPhi() << "\t" << res->getPsi() << "\t";
      if (iopts.omega) out << res->getOmega() << "\t";
      out << res->getName();
      if (iopts.printFileNames) out << "\t" << iopts.pdbfs[si];
      out << endl;
      if (iopts.freeB) {
        for (int ai = 0; ai < res->atomSize(); ai++) { (*res)[ai].setB(100*freedoms[k]); }
      }
    }
    
    // print interference
    contactList intL;
    C.getInterference(allRes, 0, &intL);
    vector<pair<Residue*, Residue*> > intList = intL.getOrderedContacts();
    for (int k = 0; k < intList.size(); k++) {
      Residue* resA = intList[k].first;
      Residue* resB = intList[k].second;
      out << "interference\t" << resA->getChainID() << "," << resA->getNum() << "\t" << resB->getChainID() << "," << resB->getNum();
      out << "\t" << std::setprecision(6) << std::fixed << intL.degree(resA, resB);
      out << "\t" << resA->getName() << "\t" << resB->getName();
      if (iopts.printFileNames) out << "\t" << iopts.pdbfs[si];
      out << endl;
    }
    
    if (iopts.seq_const) {
      // print contact degrees with sequence constraints
      contactList cL;
      C.getConstrainedContacts(allRes, 0, &cL);
      for (int k = 0; k < cL.size(); k++) {
        Residue* resA = cL.residueA(k);
        Residue* resB = cL.residueB(k);
        out << "seq_const_contact\t" << resA->getChainID() << "," << resA->getNum() << "\t" << resB->getChainID() << "," << resB->getNum();
        out << "\t" << std::setprecision(6) << std::fixed << cL.degree(k);
        out << "\t" << resA->getName() << "\t" << resB->getName();
        out << "\t" << *cL.alphabetA(k).begin() << "\t" << "XXX";
        if (iopts.printFileNames) out << "\t" << iopts.pdbfs[si];
        out << endl;
      }
    }

    // print sequence
    out << "SEQUENCE: ";
    for (int k = 0; k < allRes.size(); k++) {
      out << allRes[k]->getName() << " ";
    }
    out << endl;

    // close output file
    if (!iopts.omapfs.empty()) of.close();

    // write out the parsed region of interest
    if (!iopts.opdbfs.empty()) S.writePDB(iopts.opdbfs[si]);

    if (!iopts.rotOutFile.empty()) C.closeLogFile();
  }
  
  cout << "In total, took " << timer.getDuration() << " to complete" << endl;
}


void proteinOnly(System& S, System& So, vector<string>& legalNames) {
  AtomPointerVector A;
  for (int i = 0; i < So.chainSize(); i++) {
    Chain& chain = So[i];
    for (int j = 0; j < chain.residueSize(); j++) {
      Residue& res = chain[j];
      for (int k = 0; k < legalNames.size(); k++) {
        if (res.isNamed(legalNames[k])) {
          for (int ai = 0; ai < res.atomSize(); ai++) A.push_back(&(res[ai]));
          break;
        }
      }
    }
  }
  S.addAtoms(&A);
}
