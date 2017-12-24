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

using namespace std;
using namespace MST;

void selectAround(const vector<Residue*>& cenRes, int pm, vector<Residue*>& fragRes, vector<int>& fragResIdx) {
  Structure* S = cenRes[0]->getChain()->getParent();
  vector<bool> selected(S->residueSize(), false);
  for (int i = 0; i < cenRes.size(); i++) {
    Residue& res = *(cenRes[i]);
    Chain& C = *(res.getChain());
    int ri = res.getResidueIndex();
    for (int k = ri - pm; k <= ri + pm; k++) {
      if ((k < 0) || (k >= C.residueSize())) continue;
      selected[k] = true;
    }
  }
  for (int k = 0; k < selected.size(); k++) {
    if (selected[k]) {
      fragRes.push_back(&(S->getResidue(k)));
      fragResIdx.push_back(k);
    }
  }
}

void addMatches(FASST& F, vector<Residue*>& fragRes, vector<Structure*>& allMatches, vector<vector<Residue*> >& resTopo) {
  Structure frag(fragRes);
  F.setQuery(frag);
  F.setRMSDCutoff(RMSDCalculator::rmsdCutoff(frag));
  F.setMaxNumMatches(10);
  F.setMinNumMatches(2);
  F.search();
  fasstSolutionSet matches = F.getMatches();
  cout << "found " << matches.size() << " matches" << endl;
  for (auto it = matches.begin(); it != matches.end(); ++it) {
    allMatches.push_back(new Structure(F.getMatchStructure(*it, false, FASST::matchType::REGION)));
    Structure& match = *(allMatches.back());
    for (int k = 0; k < fragRes.size(); k++) {
      resTopo[fragRes[k]->getResidueIndex()].push_back(&(match.getResidue(k)));
    }
  }

}

int main(int argc, char** argv) {
  MstOptions op;
  op.setTitle("Starting from some arbitrary conformation of a chain, iteratively build a TERM-based compact structure. Options:");
  op.addOption("p", "starting conformation PDB file.", true);
  op.addOption("rLib", "a path to an MST-formatter rotamer library.", true);
  op.addOption("d", "a database file with a list of PDB files.");
  op.addOption("b", "a binary database file. If both --d and --b are given, will overwrite this file with a corresponding binary database.");
  op.addOption("o", "output base name.", true);
  op.addOption("rad", "compactness radius. Default will be based on protein length.");
  op.addOption("cyc", "number of iteration cycles (10 by default).");
  op.setOptions(argc, argv);
  RMSDCalculator rc;
  Structure I(op.getString("p"));
  FASST F;
  F.setMemorySaveMode(true);
  if (op.isGiven("d")) {
    F.addTargets(MstUtils::fileToArray(op.getString("d")));
    if (op.isGiven("b")) {
      F.writeDatabase(op.getString("b"));
    }
  } else if (op.isGiven("b")) {
    F.readDatabase(op.getString("b"));
  } else {
    MstUtils::error("either --b or --d must be given!");
  }
  F.pruneRedundancy(0.5);
  RotamerLibrary RL(op.getString("rLib"));
  int pmSelf = 2, pmPair = 1;
  int Ni = 1000;
  fusionScores scores;

  // TERMify loop
  Structure S = I;
  fstream out; MstUtils::openFile(out, op.getString("o") + ".traj.pdb", ios_base::out);
  out << "MODEL " << 0 << endl; S.writePDB(out); out << "ENDMDL" << endl;
  for (int c = 0; c < op.getInt("cyc", 10); c++) {
    out << "MODEL " << c+1 << endl;
    cout << "Cycle " << c+1 << "..." << endl;
    /* --- Decorate the current conformation with TERMs --- */
    // first self TERMs
    cout << "Searching for self TERMs..." << endl;
    vector<vector<Residue*> > resTopo(I.residueSize());
    for (int ri = 0; ri < S.residueSize(); ri++) resTopo[ri].push_back(&(S.getResidue(ri)));
    vector<Structure*> allMatches;
    for (int ci = 0; ci < S.chainSize(); ci++) {
      Chain& C = S[ci];
      for (int ri = 0; ri < C.residueSize(); ri++) {
        vector<Residue*> fragRes; vector<int> fragResIdx;
        selectAround({&C[ri]}, pmSelf, fragRes, fragResIdx);
        Structure frag(fragRes);
        cout << "TERM around " << C[ri] << " ";
        addMatches(F, fragRes, allMatches, resTopo);
      }
    }

    // then pair TERMs
    cout << "Searching for pair TERMs..." << endl;
    ConFind C(&RL, S);
    contactList L = C.getContacts(S, 0.01);
    vector<pair<Residue*, Residue*> > list = L.getOrderedContacts();
    for (int k = 0; k < list.size(); k++) {
      Residue* resA = list[k].first;
      Residue* resB = list[k].second;
      vector<Residue*> fragRes; vector<int> fragResIdx;
      selectAround({resA, resB}, pmPair, fragRes, fragResIdx);
      cout << "TERM around " << *resA << " x " << *resB << " ";
      addMatches(F, fragRes, allMatches, resTopo);
    }

// fstream ofs;
// MstUtils::openFile(ofs, "/tmp/allmatches.pdb", ios_base::out);
// S.writePDB(ofs);
// for (int i = 0; i < allMatches.size(); i++) allMatches[i]->writePDB(ofs);
// ofs.close();

    /* --- Fuse --- */
    fusionParams opts; opts.setNumIters(Ni); opts.setVerbose(false);
    opts.setGradientDescentFlag(true);
    opts.setRepFC(1);
    opts.setCompFC(0.1);
    mstreal compactnessRadius = op.getReal("rad", pow(I.residueSize() * 1.0, 1.0/3)*5.0);
    opts.setCompRad(compactnessRadius);
    cout << "will be trying to combine structure to a radius of " << compactnessRadius << endl;
    Structure fused = Fuser::fuse(resTopo, scores, vector<int>(), opts);
    cout << "\t" << scores << endl;
    // opts.setStartingStructure(fused);
    // opts.setGradientDescentFlag(false);
    // fused = Fuser::fuse(resTopo, scores, vector<int>(), opts);
    // cout << "\t" << scores << endl;

    /* --- write intermediate result and clean up--- */
    fused.writePDB(out); out << "ENDMDL" << endl;
    S = fused;
    for (int mi = 0; mi < allMatches.size(); mi++) delete(allMatches[mi]);
  }
  out.close();
  S.writePDB(op.getString("o") + ".fin.pdb");
}
