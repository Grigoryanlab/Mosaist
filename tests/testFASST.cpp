#include "msttypes.h"
#include "mstoptions.h"
#include "mstfasst.h"
#include <chrono>

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Implements the FASST (FAst Structure Search Algorithm). Options:");
  op.addOption("q", "query PDB file.", true);
  op.addOption("d", "a database file with a list of PDB files.");
  op.addOption("b", "a binary database file. If both --d and --b are given, will overwrite this file with a corresponding binary database.");
  op.addOption("r", "RMSD cutoff (takes the size-dependent cutoff by default).");
  op.addOption("red", "set redundancy cutoff level in percent (default is 100, so no redundancy filtering).");
  op.addOption("redProp", "set redundancy property name. If defined, will assume the FASST database encodes this relational property and will define redundancy via it.");
  op.addOption("min", "min number of matches.");
  op.addOption("max", "max number of matches.");
  op.addOption("strOut", "dump structures into this directory.");
  op.addOption("pp", "store phi/psi properties in the database, if creating a new one from PDB files.");
  op.setOptions(argc, argv);
  if (op.isGiven("redProp")) MstUtils::assert(!op.getString("redProp").empty(), "--redProp must specify a property name");

  FASST S;
  cout << "Reading the database..." << endl;
  auto begin = chrono::high_resolution_clock::now();
  Structure query(op.getString("q"));
  S.setQuery(query);
  S.setMemorySaveMode(true);
  if (op.isGiven("d")) {
    vector<string> pdbFiles = MstUtils::fileToArray(op.getString("d"));
    for (int i = 0; i < pdbFiles.size(); i++) {
      Structure P(pdbFiles[i]);
      S.addTarget(P);
      // compute and add some properties
      if (op.isGiven("pp")) {
        vector<Residue*> residues = P.getResidues();
        vector<mstreal> phi(residues.size()), psi(residues.size());
        for (int ri = 0; ri < residues.size(); ri++) {
          phi[ri] = residues[ri]->getPhi(false);
          psi[ri] = residues[ri]->getPsi(false);
        }
        S.addResidueProperties(S.numTargets() - 1, "phi", phi);
        S.addResidueProperties(S.numTargets() - 1, "psi", psi);
      }
    }
    if (op.isGiven("b")) {
      S.writeDatabase(op.getString("b"));
    }
  } else if (op.isGiven("b")) {
    S.readDatabase(op.getString("b"));
  } else {
    MstUtils::error("either --b or --d must be given!");
  }
  if (op.isGiven("r")) { S.setRMSDCutoff(op.getReal("r")); }
  else {
    cout << "setting RMSD cutoff to " << RMSDCalculator::rmsdCutoff(query) << endl;
    S.setRMSDCutoff(RMSDCalculator::rmsdCutoff(query));
  }
  S.setMaxNumMatches(op.getInt("max", -1));
  S.setMinNumMatches(op.getInt("min", -1));
  // S.setMaxGap(0, 1, 6); S.setMinGap(0, 1, 0);
  S.setRedundancyCut(op.getReal("red", 100.0)/100.0);
  if (op.isGiven("redProp")) S.setRedundancyProperty(op.getString("redProp"));
  auto end = chrono::high_resolution_clock::now();
  cout << "DB reading took " << chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " ms" << endl;
  cout << "Searching..." << endl;
  begin = chrono::high_resolution_clock::now();
  S.search();
  end = chrono::high_resolution_clock::now();
  cout << "Search took " << chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " ms" << endl;
  cout << "found " << S.numMatches() << " matches:" << endl;
  fasstSolutionSet matches = S.getMatches(); int i = 0;
  vector<vector<mstreal> > phi, psi;
  for (auto it = matches.begin(); it != matches.end(); ++it, ++i) {
    cout << S.toString(*it) << endl;
    cout << *it << endl;
    if (op.isGiven("pp")) {
      if (S.isPropertyDefined("phi")) cout << "\tphi: " << MstUtils::vecToString(S.getResidueProperties(*it, "phi")) << endl;
      if (S.isPropertyDefined("psi")) cout << "\tpsi: " << MstUtils::vecToString(S.getResidueProperties(*it, "psi")) << endl;
    }
    if (op.isGiven("strOut")) {
      // Structure match = S.getMatchStructure(*it, true, FASST::matchType::FULL);
      Structure match = S.getMatchStructure(*it);
      match.writePDB(op.getString("strOut") + "/match" + MstUtils::toString(i) + ".pdb");
    }
  }
  for (auto it = matches.begin(); it != matches.end(); ++it, ++i) {
    Sequence seq = S.getMatchSequence(*it);
    cout << seq << endl;
  }
}
