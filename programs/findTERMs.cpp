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
  op.addOption("minSegs", "minimal number of segments in the TERMs to consider.");
  op.addOption("maxSegs", "maximal number of segments in the TERMs to consider");
  op.setOptions(argc, argv);
  if (op.isGiven("max")) MstUtils::assertCond(op.isReal("max") && (op.getReal("max") > 0), "--max must be a positive number");
  if (op.isGiven("L")) MstUtils::assertCond(op.isInt("L") && (op.getInt("L") > 0), "--L must be a positive integer");
  if (op.isGiven("minSegs")) MstUtils::assertCond(op.isInt("minSegs") && (op.getInt("minSegs") > 0), "--minSegs must be a positive integer");
  if (op.isGiven("maxSegs")) MstUtils::assertCond(op.isInt("maxSegs") && (op.getInt("maxSegs") > 0), "--maxSegs must be a positive integer");
  mstreal max = op.getReal("max", 1.1);
  int L0 = op.getInt("L", 15);
  int maxSegs = op.getInt("maxSegs", 10E6);
  int minSegs = op.getInt("minSegs", 0);

  Structure T(op.getString("p"));
  FASST S;
  RMSDCalculator rc;
  S.addTarget(T);
  vector<string> termPaths = MstUtils::fileToArray(op.getString("terms"));
  fstream ofs; MstUtils::openFile(ofs, op.getString("o"), fstream::out);
  T.writePDB(ofs);
  int Nm = 0;
  for (int i = 0; i < termPaths.size(); i++) {
    Structure term(termPaths[i]);
    S.stripSidechains(term);
    term = term.reassignChainsByConnectivity();
    if ((term.chainSize() > maxSegs) || (term.chainSize() < minSegs)) continue;
    S.setRMSDCutoff(RMSDCalculator::rmsdCutoff(term, max, L0));
    S.setQuery(term);
    S.search();
    fasstSolutionSet matches = S.getMatches();
    for (auto it = matches.begin(); it != matches.end(); ++it) {
      Structure match = S.getMatchStructure(*it, false, FASST::matchType::REGION, false);
      Structure termAligned(term);
      rc.align(termAligned.getAtoms(), match.getAtoms(), termAligned);
      ofs << "REM TERM " << termPaths[i] << endl;
      termAligned.writePDB(ofs);
    }
    Nm += matches.size();
    if (i % 100 == 0) cout << i+1 << "/" << termPaths.size() << ", " << Nm << " matches" << endl;
  }
  ofs.close();
}
