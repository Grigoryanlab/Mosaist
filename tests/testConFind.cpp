#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <typeinfo>
#include <stdio.h>
#include <getopt.h>
#include <map>
#include <iomanip>
#include <unistd.h>
#include <ctime>

#include "msttypes.h"
#include "mstrotlib.h"
#include "mstcondeg.h"

using namespace std;
using namespace MST;

// ---- Utility Functions
void proteinOnly(System& S, System& So, vector<string>& legalNames);

class options {
  public:
    options() {
      dcut = 25.0;
      clashDist = 2.0; contDist = 3.0; rotLibFile = ""; beblFile = ""; rotOutFile = ""; rotLevel = "";
      verbose = renumPDB = phi_psi = omega = printFileNames = false; calcContacts = true;
      aaProp["ALA"] = 7.73; aaProp["CYS"] = 1.84; aaProp["ASP"] = 5.82; aaProp["GLU"] = 6.61; aaProp["PHE"] = 4.05;
      aaProp["GLY"] = 7.11; aaProp["HIS"] = 2.35; aaProp["HSD"] = 2.35; aaProp["ILE"] = 5.66; aaProp["LYS"] = 6.27;
      aaProp["LEU"] = 8.83; aaProp["MET"] = 2.08; aaProp["ASN"] = 4.50; aaProp["PRO"] = 4.52; aaProp["GLN"] = 3.94;
      aaProp["ARG"] = 5.03; aaProp["SER"] = 6.13; aaProp["THR"] = 5.53; aaProp["VAL"] = 6.91; aaProp["TRP"] = 1.51; aaProp["TYR"] = 3.54;
    }
    vector<string> pdbfs, omapfs, opdbfs;
    bool verbose, renumPDB, phi_psi, printFileNames, omega, calcContacts;
    string selection, focus, rotLibFile, beblFile, rotOutFile, rotLevel;
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

void usage() {
  int w = 80, p1 = 3, p2 = p1+8;
  cout << endl << option("", "Identifies inter-positional contacts and environment information from input PDB file(s). Options:", w, 0, 0) << endl;
  cout << option("--p", "input PDB file.", w, p1, p2) << endl;
  cout << option("--pL", "a file with a list of PDB files. Either --p or --pL must be given.", w, p1, p2) << endl;
  cout << option("--o", "output file name for writing contacts. If not given, will write to standard output.", w, p1, p2) << endl;
  cout << option("--oL", "a file with a list of contact file names (one per input PDB structure).", w, p1, p2) << endl;
  cout << option("--opdb", "optional: output post-processed PDB file (useful for keeping track of how all PDB weirdnesses got parsed).", w, p1, p2) << endl;
  cout << option("--opdbL", "optional: a file with a list of file names for post-processed PDBs, one per input structure.", w, p1, p2) << endl;
  cout << option("--rLib", "a path to an MSL-formatter rotamer library, with WEIGHTS information. If a file, will be treated as a backbone-independent library. If a directory, will look for files EBL.out and BEBL.out in it to read a backbone-dependent library.", w, p1, p2) << endl;
  cout << option("--rout", "name of a file into which to place PDB-formated coordinates of rotamers that ended up surviving at each considered position.", w, p1, p2) << endl;
  cout << option("--psel", "optional: pre-selection string to apply before doing anything (only the selected part of structure will be considered). Will select a residue if its CA atom is included in the given selection. E.g., 'NAME CA WITHIN 25 OF CHAIN A'.", w, p1, p2) << endl;
  cout << option("--sel", "optional: selection string for defining which residues to compute properties for. This is different from --psel. For example, one could pre-select a certain chain with --psel, but only be concerned with data about specific resdues within that chain, by defining --sel. The entire chain would be preserved in the structure, so that (for example) backbone rotameric collisions are properly detected. However, time would not be wasted on analyzing residues not in the --sel set. E.g., 'RESID 20-30'. Again, a residue is included if its CA is in the selection. NOTE: typically, values for some of the residues in --sel will be wrong (e.g., sum contact degree or crowdedness; because relevant neighboring residues will be missing from --sel). Thus, in practice --sel will need to include both residues of interest AND important surrounding residues, for which values will be wrong, but they will assure that values for the set of interest are correct.", w, p1, p2) << endl;
  cout << option("--pp", "optional: print phi/psi angles for each residue (will print next to all positional scores score).", w, p1, p2) << endl;
  cout << option("--omg", "optional: print omega angles for each residue (will print next to all positional scores score). Omega for the current position is defined as the -CA, C, N, CA dihedral angle.", w, p1, p2) << endl;
  cout << option("--verb", "optional: generate lots of detailed output (i.e., print the details of which rotamer pairs are in contact).", w, p1, p2) << endl;
  cout << option("--pf", "optional: if flag specified, will print the name of the PDB file being analyzed next to all positional scores. This is especially convenient when a list of PDB file is specified as input and the output goes to a single file.", w, p1, p2) << endl;
  cout << option("--ren", "optional: if flag specified, will renumber the structure before output. Useful for keeping track of residues in the output list of contacts if the input PDB file is strangely numbered.", w, p1, p2) << endl;
  cout << option("--nc", "optional: if flag specified, contact information will not be calculated/printed and only self information will.", w, p1, p2) << endl;
}

void parseCommandLine(int argc, char** argv, options& iopts) {
  map<string, bool> spec;

  while (1) {
    int oind = 0;
    static struct option opts[] = {
      {"p", 1, 0, 1},
      {"o", 1, 0, 2},
      {"opdb", 1, 0, 7},
      {"pL", 1, 0, 10},
      {"oL", 1, 0, 11},
      {"opdbL", 1, 0, 12},
      {"psel", 1, 0, 18},
      {"sel", 1, 0, 24},
      {"rLib", 1, 0, 19},
      {"pp", 0, 0, 23},
      {"omg", 0, 0, 26},
      {"verb", 0, 0, 20},
      {"rout", 1, 0, 21},
      {"ren", 0, 0, 22},
      {"pf", 0, 0, 25},
      {"nc", 0, 0, 27},
      {0, 0, 0, 0}
    };

    int c = getopt_long (argc, argv, "", opts, &oind);
    if (c == -1) break;

    switch (c) {
      case 1:
        if (iopts.pdbfs.size() > 0) { usage(); exit(-1); }
        iopts.pdbfs.push_back(optarg);
        spec[string(opts[oind].name)] = true;
        break;

      case 2:
        if (iopts.omapfs.size() > 0) { usage(); exit(-1); }
        iopts.omapfs.push_back(optarg);
        spec[string(opts[oind].name)] = true;
        break;

      case 7:
        if (iopts.opdbfs.size() > 0) { usage(); exit(-1); }
        iopts.opdbfs.push_back(optarg);
        spec[string(opts[oind].name)] = true;
        break;

      case 10:
        if (iopts.pdbfs.size() > 0) { usage(); exit(-1); }
        MstUtils::fileToArray(optarg, iopts.pdbfs);
        spec[string(opts[oind].name)] = true;
        break;

      case 11:
        if (iopts.omapfs.size() > 0) { usage(); exit(-1); }
        MstUtils::fileToArray(optarg, iopts.omapfs);
        spec[string(opts[oind].name)] = true;
        break;

      case 12:
        if (iopts.opdbfs.size() > 0) { usage(); exit(-1); }
        MstUtils::fileToArray(optarg, iopts.opdbfs);
        spec[string(opts[oind].name)] = true;
        break;

      case 18:
        iopts.selection = string(optarg);
        spec[string(opts[oind].name)] = true;
        break;

      case 24:
        iopts.focus = string(optarg);
        spec[string(opts[oind].name)] = true;
        break;

      case 19:
        if (MstUtils::isDir(optarg)) {
          iopts.rotLibFile = string(optarg) + "/EBL.out";
          iopts.beblFile = string(optarg) + "/BEBL.out";
        } else {
          iopts.rotLibFile = string(optarg);
        }
        spec[string(opts[oind].name)] = true;
        break;

      case 20:
        iopts.verbose = true;
        spec[string(opts[oind].name)] = true;
        break;

      case 23:
        iopts.phi_psi = true;
        spec[string(opts[oind].name)] = true;
        break;

      case 26:
        iopts.omega = true;
        spec[string(opts[oind].name)] = true;
        break;

      case 21:
        iopts.rotOutFile = string(optarg);
        spec[string(opts[oind].name)] = true;
        break;

      case 22:
        iopts.renumPDB = true;
        spec[string(opts[oind].name)] = true;
        break;

      case 25:
        iopts.printFileNames = true;
        spec[string(opts[oind].name)] = true;
        break;

      case 27:
        iopts.calcContacts = false;
        spec[string(opts[oind].name)] = true;
        break;

      case '?':
        break;

      default:
        printf ("?? getopt returned character code 0%o ??\n", c);
    }
  }

  // make sure all required options have been specified
  if (!(((spec.find(string("p")) != spec.end()) || (spec.find(string("pL")) != spec.end())) && (spec.find(string("rLib")) != spec.end()))) {
    cout << "Not all required options specified!\n"; usage(); exit(-1);
  }

  // error checking
  // make sure lists are of the proper size
  if (iopts.omapfs.size() > 1) MstUtils::assert(iopts.omapfs.size() == iopts.pdbfs.size(), "the number of input PDB files and output map files does not agree");
  if (iopts.opdbfs.size() > 1) MstUtils::assert(iopts.opdbfs.size() == iopts.pdbfs.size(), "the number of input PDB files and output PDB files does not agree");
  if (iopts.pdbfs.size() > 1) {
    if (iopts.omapfs.size() == 1) {
      iopts.omapfs.resize(iopts.pdbfs.size());
      string base = MstUtils::pathBase(iopts.omapfs[0]);
      for (int i = 0; i < iopts.pdbfs.size(); i++) { iopts.omapfs[i] = base + ".f" + MstUtils::toString(i+1) + ".cont"; }
    }
    if (iopts.opdbfs.size() == 1) {
      iopts.opdbfs.resize(iopts.pdbfs.size());
      string base = MstUtils::pathBase(iopts.opdbfs[0]);
      for (int i = 0; i < iopts.pdbfs.size(); i++) { iopts.opdbfs[i] = base + ".f" + MstUtils::toString(i+1) + ".pdb"; }
    }
  }
}

// ---- Main Program
int main(int argc, char *argv[]) {
  int ii, jj; double d;
  fstream rotOut, *rotOutPtr = NULL;

  // process input arguments
  options iopts;
  parseCommandLine(argc, argv, iopts);
  if (!iopts.rotOutFile.empty()) { MstUtils::openFile(rotOut, iopts.rotOutFile, fstream::out); rotOutPtr = &rotOut; }

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
    Structure S;                                         // just the region of the original structure corresponding to the map
    proteinOnly(S, So, legalNames);
    if (iopts.renumPDB) S.renumber();
    ConFind C(&RL, S); // this structure will be used for computing contact degrees

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

    // --- compute contact degrees
    contactList L = C.getContacts(S);
    for (int k = 0; k < L.size(); k++) {
      cout << *(L.residueA(k)) << " - " << *(L.residueB(k)) << " " << L.degree(k) << endl;
    }
    //
    //   // --- write contact degree information
    //   vector<double> sumContDeg(resIndex.size(), 0);
    //   for (int i = 0; i < conts.size(); i++) {
    //     ii = conts[i].resi;
    //     jj = conts[i].resj;
    //     d = conts[i].degree;
    //     sumContDeg[ii] += d;
    //     sumContDeg[jj] += d;
    //     out << contactString(S, resIndex[ii], resIndex[jj], d) << endl;
    //     if (iopts.verbose) { printf("%s", conts[i].info.c_str()); }
    //   }
    //
    //   // -- write out sum contact degrees
    //   for (int i = 0; i < sumContDeg.size(); i++) {
    //     out << "sumcond\t" << S.getPosition(resIndex[i]).getPositionId() << "\t" << std::setprecision(6) << std::fixed << sumContDeg[i];
    //     if (iopts.phi_psi) out << "\t" << pp[i].first << "\t" << pp[i].second;
    //     if (iopts.omega) out << "\t" << pp[i].third;
    //     out << "\t" << (S.getPosition(resIndex[i])).getResidueName();
    //     if (iopts.printFileNames) out << "\t" << iopts.pdbfs[si];
    //     out << endl;
    //   }
    //
    //   // -- write out the freedom parameter
    //   for (int i = 0; i < freedom.size(); i++) {
    //     out << "freedom\t" << S.getPosition(resIndex[i]).getPositionId() << "\t" << std::setprecision(6) << std::fixed << freedom[i];
    //     if (iopts.phi_psi) out << "\t" << pp[i].first << "\t" << pp[i].second;
    //     if (iopts.omega) out << "\t" << pp[i].third;
    //     out << "\t" << (S.getPosition(resIndex[i])).getResidueName();
    //     if (iopts.printFileNames) out << "\t" << iopts.pdbfs[si];
    //     out << endl;
    //   }
    // }
    //
    // // -- write out permanent contacts
    // for (int i = 0; i < permanentContacts.size(); i++) {
    //   for (set<int>::iterator it = permanentContacts[i].begin(); it != permanentContacts[i].end(); ++it) {
    //     out << contactString(S, resIndex[i], resIndex[*it], -1, true) << endl;
    //   }
    // }
    //
    // // -- write out the crowdedness parameter
    // for (int i = 0; i < fractionPruned.size(); i++) {
    //   out << "crwdnes\t" << S.getPosition(resIndex[i]).getPositionId() << "\t" << std::setprecision(6) << std::fixed << fractionPruned[i];
    //   if (iopts.phi_psi) out << "\t" << pp[i].first << "\t" << pp[i].second;
    //   if (iopts.omega) out << "\t" << pp[i].third;
    //   out << "\t" << (S.getPosition(resIndex[i])).getResidueName();
    //   if (iopts.printFileNames) out << "\t" << iopts.pdbfs[si];
    //   out << endl;
    // }
    //
    // // -- write out free volume parameter
    // for (int i = 0; i < freeVolume.size(); i++) {
    //   out << "freevol\t" << S.getPosition(resIndex[i]).getPositionId() << "\t" << std::setprecision(6) << std::fixed << freeVolume[i];
    //   if (iopts.phi_psi) out << "\t" << pp[i].first << "\t" << pp[i].second;
    //   if (iopts.omega) out << "\t" << pp[i].third;
    //   out << "\t" << (S.getPosition(resIndex[i])).getResidueName();
    //   if (iopts.printFileNames) out << "\t" << iopts.pdbfs[si];
    //   out << endl;
    // }
    //
    // // --- free near-neighbor structures if was dealing with contact probability maps
    // for (int i = 0; i < rotamers.size(); i++) {
    //   for (int j = 0; j < rotamers[i].size(); j++) {
    //     delete(rotamers[i][j].gridSC());
    //     delete(rotamers[i][j].gridBB());
    //   }
    // }
    //
    // // write sequence information
    // out << "SEQUENCE:";
    // for (int i = 0; i < resIndex.size(); i++) {
    //   Position &p = S.getPosition(resIndex[i]);
    //   out << " " << p.getResidueName();
    // }
    // out << endl;

    // close output file
    if (!iopts.omapfs.empty()) of.close();

    // write out the parsed region of interest
    if (!iopts.opdbfs.empty()) S.writePDB(iopts.opdbfs[si]);

  }
  if (!iopts.rotOutFile.empty()) rotOut.close();

}


// ---- utility functions definitions
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
