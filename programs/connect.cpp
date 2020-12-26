#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>
#include <typeinfo>
#include <stdio.h>

#include "msttypes.h"
#include "mstfasst.h"
#include "mstcondeg.h"
#include "mstfuser.h"
#include "mstcondeg.h"
#include "mstoptions.h"
#include "mstmagic.h"
#include "mstlinalg.h"

using namespace std;
using namespace MST;

Structure standardize(const Structure& input) {
  Structure output;
  RotamerLibrary::extractProtein(output, input);
  if (!RotamerLibrary::hasFullBackbone(output)) MstUtils::error("inputs must have a fully defined backbone");
  return output;
}

int main(int argc, char** argv) {
  MstOptions op;
  op.setTitle("Connects two structures (of equal length) via a smooth transition in conformation space, with each point being well populated with examples in the database. Options:");
  op.addOption("beg", "first structure (start of the trajectory).", true);
  op.addOption("end", "last structure (end of the trajectory).", true);
  op.addOption("b", "binary FASST database to use.", true);
  op.addOption("o", "output base name.", true);
  op.addOption("tol", "RMSD tolerance to consider the trajectory as having arrived at the last structure.");
  op.setOptions(argc, argv);
  RMSDCalculator rc;

  // read structures and standardize/remove side-chains
  Structure B = standardize(Structure(op.getString("beg")));
  Structure E = standardize(Structure(op.getString("end")));
  Structure C = B;

  // read FASST database
  FASST F;
  F.readDatabase(op.getString("b"), 1);
  F.setRedundancyCut(0.5);

  // start the trajectory output
  int c = 0;
  fstream out; MstUtils::openFile(out, op.getString("o") + ".traj.pdb", ios::out);
  out << "MODEL " << c << endl; B.writePDB(out); out << "ENDMDL" << endl;

  // find the closest thing in the database to the end structure (this defines the convergence tolerance)
  F.setQuery(E);
  F.setRMSDCutoff(999.0);
  F.setMaxNumMatches(1);
  F.setMinNumMatches(1);
  fasstSolutionSet matchList = F.search();
  mstreal tol = matchList.bestRMSD();
  cout << "closest match to the final structure in the database has RMSD of " << tol << endl;

  // the connect loop
  int numMatches = 100;
  cout << "[" << numMatches << "-th match RMSD] [RMSD to final]" << endl;
  while (rc.bestRMSD(RotamerLibrary::getBackbone(C), RotamerLibrary::getBackbone(E)) > op.getReal("tol", 1.1*tol)) {
    // search for the closest numMatches to the current conformation
    F.setQuery(C);
    F.setRMSDCutoff(999.0);
    F.setMaxNumMatches(numMatches);
    F.setMinNumMatches(numMatches);
    fasstSolutionSet matchList = F.search();

    // identify the one closest to the final point
    vector<mstreal> rmsds(matchList.size());
    for (int i = 0; i < matchList.size(); i++) {
      rmsds[i] = rc.bestRMSD(RotamerLibrary::getBackbone(E), F.getMatchStructure(matchList[i]).getAtoms());
    }
    int mi; MstUtils::min(rmsds, -1, -1, &mi);
    cout << matchList.worstRMSD() << " " << rmsds[mi] << endl;
    C = F.getMatchStructure(matchList[mi]);

    out << "MODEL " << c << endl;
    C.writePDB(out); out << "ENDMDL" << endl;
    c++;
  }

  // end trajectory output
  out << "MODEL " << c << endl;
  rc.align(RotamerLibrary::getBackbone(E), RotamerLibrary::getBackbone(C), E);
  E.writePDB(out); out << "ENDMDL" << endl;
  out.close();

  return 0;
}
