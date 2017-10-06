#include "msttypes.h"
#include "mstoptions.h"
#include "mstfasst.h"
#include <chrono>

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Implements the FASSA (FAst Structure Search Algorithm). Options:");
  op.addOption("q", "query PDB file.", true);
  op.addOption("d", "a database file with a list of PDB files.");
  op.addOption("b", "a binary database file. If both --d and --b are given, will overwrite this file with a corresponding binary database.");
  op.addOption("r", "RMSD cutoff.", true);
  op.addOption("min", "min number of matches.");
  op.addOption("max", "max number of matches.");
  op.setOptions(argc, argv);

  FASST S;
  cout << "Reading the database..." << endl;
  auto begin = chrono::high_resolution_clock::now();
  S.setQuery(op.getString("q"));
  S.setMemorySaveMode(true);
  if (op.isGiven("d")) {
    S.addTargets(MstUtils::fileToArray(op.getString("d")));
    if (op.isGiven("b")) {
      S.writeDatabase(op.getString("b"));
    }
  } else if (op.isGiven("b")) {
    S.readDatabase(op.getString("b"));
  } else {
    MstUtils::error("either --b or --d must be given!");
  }
  S.setRMSDCutoff(op.getReal("r"));
  S.setMaxNumMatches(op.getInt("max", -1));
  S.setMinNumMatches(op.getInt("min", -1));
  // S.setMaxGap(0, 1, 6); S.setMinGap(0, 1, 1);
  auto end = chrono::high_resolution_clock::now();
  cout << "DB reading took " << chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " ms" << endl;
  cout << "Searching..." << endl;
  begin = chrono::high_resolution_clock::now();
  S.search();
S.setQuery(op.getString("q"));
S.search();
  end = chrono::high_resolution_clock::now();
  cout << "Search took " << chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " ms" << endl;
  cout << "found " << S.numMatches() << " matches:" << endl;
  set<fasstSolution> matches = S.getMatches(); int i = 0;
  for (auto it = matches.begin(); it != matches.end(); ++it, ++i) {
    cout << S.toString(*it) << endl;
    Structure match = S.getMatchStructure(*it, true, FASST::matchType::FULL);
//    match.writePDB("/tmp/match" + MstUtils::toString(i) + ".pdb");
  }
}
