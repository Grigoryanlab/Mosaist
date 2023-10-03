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
      verbose = renumPDB = phi_psi = omega = printFileNames = false;
      aaProp["ALA"] = 7.73; aaProp["CYS"] = 1.84; aaProp["ASP"] = 5.82; aaProp["GLU"] = 6.61; aaProp["PHE"] = 4.05;
      aaProp["GLY"] = 7.11; aaProp["HIS"] = 2.35; aaProp["HSD"] = 2.35; aaProp["ILE"] = 5.66; aaProp["LYS"] = 6.27;
      aaProp["LEU"] = 8.83; aaProp["MET"] = 2.08; aaProp["ASN"] = 4.50; aaProp["PRO"] = 4.52; aaProp["GLN"] = 3.94;
      aaProp["ARG"] = 5.03; aaProp["SER"] = 6.13; aaProp["THR"] = 5.53; aaProp["VAL"] = 6.91; aaProp["TRP"] = 1.51; aaProp["TYR"] = 3.54;
    }
    vector<string> pdbfs, omapfs, opdbfs;
    bool verbose, renumPDB, phi_psi, printFileNames, omega;
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
  cout << option("--sel", "optional: a selection string for limiting the residues for which properties will be computed.", w, p1, p2) << endl;
  cout << option("--pp", "optional: print phi/psi angles for each residue (will print next to all positional scores score).", w, p1, p2) << endl;
  cout << option("--omg", "optional: print omega angles for each residue (will print next to all positional scores score). Omega for the current position is defined as the -CA, C, N, CA dihedral angle.", w, p1, p2) << endl;
  cout << option("--verb", "optional: generate lots of detailed output (i.e., print the details of which rotamer pairs are in contact).", w, p1, p2) << endl;
  cout << option("--pf", "optional: if flag specified, will print the name of the PDB file being analyzed next to all positional scores. This is especially convenient when a list of PDB file is specified as input and the output goes to a single file.", w, p1, p2) << endl;
  cout << option("--ren", "optional: if flag specified, will renumber the structure before output. Useful for keeping track of residues in the output list of contacts if the input PDB file is strangely numbered.", w, p1, p2) << endl;
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
      {"sel", 1, 0, 24},
      {"rLib", 1, 0, 19},
      {"pp", 0, 0, 23},
      {"omg", 0, 0, 26},
      {"verb", 0, 0, 20},
      {"rout", 1, 0, 21},
      {"ren", 0, 0, 22},
      {"pf", 0, 0, 25},
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

      case 24:
        iopts.focus = string(optarg);
        spec[string(opts[oind].name)] = true;
        break;

      case 19:
        if (MstSys::isDir(optarg)) {
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
}

int maxIndex(map<string, int>& aaToIndex) {
  int mi = 0;
  for (map<string, int>::iterator it = aaToIndex.begin(); it != aaToIndex.end(); ++it) {
    if (it->second > mi) mi = it->second;
  }
  return mi;
}

void addName(vector<string>& legalNames, map<string, int>& aaToIndex, vector<string> newNames) {
  int newIdx = maxIndex(aaToIndex) + 1;
  for (int i = 0; i < newNames.size(); i++) {
    legalNames.push_back(newNames[i]);
    aaToIndex[newNames[i]] = newIdx;
  }
}

// ---- Main Program
int main(int argc, char *argv[]) {
  mstreal lcpMin = 0.5;       // the lowest low collision probability (CP) cutoff value to consider
  mstreal lcpMax = 1.5;       // the highest low CP cutoff value to consider
  mstreal hcpMin = 1.5;       // the lowest high CP cutoff value to consider
  mstreal hcpMax = 2.5;       // the highest high CP cutoff value to consider
  int N = 21;              // number of grid points between min and max in each parameter
  int Nb = 50;             // number of freedom bins

  vector<mstreal> lcpVals(N), hcpVals(N);
  vector<vector<vector<mstreal > > > freedoms; // freedom value at each residue for each parameter combo
  vector<int> AA;                           // aa at each position
  freedoms.resize(N);
  for (int i = 0; i < N; i++) freedoms[i].resize(N);
  for (int i = 0; i < N; i++) {
    lcpVals[i] = lcpMin + (lcpMax - lcpMin)/(N-1)*i;
    hcpVals[i] = hcpMin + (hcpMax - hcpMin)/(N-1)*i;
  }

  // process input arguments
  options iopts;
  parseCommandLine(argc, argv, iopts);

  // legal residue names that are considered "protein" here
  vector<string> legalNames;
  map<string, int> aaToIndex;
  addName(legalNames, aaToIndex, {"ALA"}); addName(legalNames, aaToIndex, {"CYS", "CSO", "SEC"});
  addName(legalNames, aaToIndex, {"ASP"}); addName(legalNames, aaToIndex, {"GLU"}); addName(legalNames, aaToIndex, {"PHE"});
  addName(legalNames, aaToIndex, {"GLY"}); addName(legalNames, aaToIndex, {"HIS", "HIP", "HSE", "HSC", "HSP", "HSD"});
  addName(legalNames, aaToIndex, {"ILE"}); addName(legalNames, aaToIndex, {"LYS"}); addName(legalNames, aaToIndex, {"LEU"});
  addName(legalNames, aaToIndex, {"MET", "MSE"}); addName(legalNames, aaToIndex, {"ASN"}); addName(legalNames, aaToIndex, {"PRO"});
  addName(legalNames, aaToIndex, {"GLN"}); addName(legalNames, aaToIndex, {"ARG"}); addName(legalNames, aaToIndex, {"SER", "SEP"});
  addName(legalNames, aaToIndex, {"THR", "TPO"}); addName(legalNames, aaToIndex, {"VAL"}); addName(legalNames, aaToIndex, {"TRP"});
  addName(legalNames, aaToIndex, {"TYR", "PTR"});
  int NAATypes = maxIndex(aaToIndex) + 1;

  // pre-read rotamer library once
  RotamerLibrary* RL = new RotamerLibrary(iopts.rotLibFile);

  // go through all PDB files
  for (int si = 0; si < iopts.pdbfs.size(); si++) {
    cout << iopts.pdbfs[si] << "..." << endl;
    Structure So(iopts.pdbfs[si]);                       // original input PDB structure
    Structure S;                                         // just the region of the original structure corresponding to the map
    proteinOnly(S, So, legalNames);
    if (iopts.renumPDB) S.renumber();
    ConFind C(RL, S);

    // optionally select the relevant residue subset
    vector<Residue*> allRes;
    if (!iopts.focus.empty()) {
      selector sel(S);
      allRes = sel.selectRes(iopts.focus);
    } else {
      allRes = S.getResidues();
    }
    C.getContacts(allRes);

    // consider all parameter combinations
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        C.setFreedomParams(lcpVals[i], hcpVals[j], 3);
        C.clearFreedom();
        vector<mstreal> structFreedoms = C.getFreedom(allRes);
        freedoms[i][j].insert(freedoms[i][j].end(), structFreedoms.begin(), structFreedoms.end());
      }
    }

    // save amino-acid sequence
    for (int k = 0; k < allRes.size(); k++) {
      AA.push_back(aaToIndex[allRes[k]->getName()]);
    }
  }
  delete RL;

  // now find the combination of parameters that maximizes the expected
  // freedom-based sequence recovery and find the
  mstreal maxSeqRec = 0; int bestI = -1; int bestJ = -1;
  int pc = 1;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      vector<mstreal>& F = freedoms[i][j];
      vector<vector<int> > aaDist(Nb, vector<int>(NAATypes, 0));
      for (int k = 0; k < F.size(); k++) {
        int bi = int(MstUtils::min((mstreal) F[k], (mstreal) (1 - 1.0/(2*Nb))*Nb)); // since freedom is always [0, 1]
        aaDist[bi][AA[k]]++;
      }
      mstreal expSeqRec = 0; int Nall = 0;
      for (int bi = 0; bi < Nb; bi++) {
        int Ntot = 0;
        for (int aai = 0; aai < NAATypes; aai++) Ntot += aaDist[bi][aai];
        Nall += Ntot;
        mstreal binExpSeqRec = 0;
        for (int aai = 0; aai < NAATypes; aai++) {
          binExpSeqRec += ((aaDist[bi][aai] + pc) * 1.0 / (Ntot + pc*NAATypes)) * ((aaDist[bi][aai] + pc) * 1.0 / (Ntot + pc*NAATypes));
        }
        expSeqRec += binExpSeqRec * Ntot;
      }
      expSeqRec /= Nall;
      cout << "--> expected sequence recovery is " << expSeqRec << " for parameter combo " << lcpVals[i] << " x " << hcpVals[j] << endl;
      if (expSeqRec > maxSeqRec) {
        maxSeqRec = expSeqRec;
        bestI = i; bestJ = j;
      }
    }
  }
  cout << "BEST: expected sequence recovery is " << maxSeqRec << " for parameter combo " << lcpVals[bestI] << " x " << hcpVals[bestJ] << endl;
  cout << "Histogram:" << endl;
  vector<mstreal>& F = freedoms[bestI][bestJ];
  vector<int> hist(Nb, 0);
  for (int k = 0; k < F.size(); k++) {
    int bi = int(MstUtils::min((mstreal) F[k], (mstreal) (1 - 1.0/(2*Nb))*Nb)); // since freedom is always [0, 1]
    hist[bi]++;
  }
  for (int bi = 0; bi < Nb; bi++) cout << hist[bi] << endl;
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
