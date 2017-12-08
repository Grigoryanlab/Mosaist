#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>
#include <typeinfo>
#include <stdio.h>

#include "msttypes.h"
#include "mstfasst.h"
#include "mstoptions.h"

using namespace std;
using namespace MST;

int main(int argc, char** argv) {
  MstOptions op;
  op.setTitle("Find all TERM instances within a given structure. Options:");
  op.addOption("p", "input PDB file.", true);
  op.addOption("terms", "a file with a list of PDB files, each corresponding to a single TERM centroid.", true);
  op.addOption("o", "output file name, which will include the original structure as well as all of the identified TERM matches as independent objects.", true);
  op.addOption("max", "RMSD cutoff plateau (default is 1.1).");
  op.addOption("L", "RMSD cutoff persistance length (default is 15).");
  op.setOptions(argc, argv);
  MstUtils::assert(op.isReal("max") && (op.getReal("max") > 0), "--max must be a positive number");
  MstUtils::assert(op.isInt("L") && (op.getInt("L") > 0), "--L must be a positive integer");
  mstreal max = op.getReal("max", 1.1);
  int L0 = op.getInt("L", 15);

  Structure T(op.getString("p"));
  FASST S;
  S.addTarget(T);
  vector<string> termPaths = MstUtils::fileToArray(op.getString("terms"));
  fstream ofs; MstUtils::openFile(ofs, op.getString("o"), fstream::out);
  Structure(T).writePDB(ofs);
  for (int i = 0; i < termPaths.size(); i++) {
    Structure term(termPaths[i]);
    cout << termPaths[i] << "..." << endl;
    S.setRMSDCutoff(RMSDCalculator::rmsdCutoff(term, max, L0));
    S.setQuery(term);
    S.search();
    fasstSolutionSet matches = S.getMatches();
    for (auto it = matches.begin(); it != matches.end(); ++it, ++i) {
      Structure match = S.getMatchStructure(*it);
      match.writePDB(ofs);
    }
  }
  ofs.close();

}
