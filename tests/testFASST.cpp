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
  cout << "Building the database..." << endl;
  auto beginDB = chrono::high_resolution_clock::now();
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
  S.setMaxGap(0, 1, 10); S.setMinGap(0, 1, 0);
  auto endDB = chrono::high_resolution_clock::now();
  cout << "Took " << chrono::duration_cast<std::chrono::milliseconds>(endDB-beginDB).count() << " ms" << endl;
  cout << "Searching..." << endl;
  if (false) {
    int N = 10;
    for (int i = 0; i < N; i++) {
      mstreal d = i*20.0/(N-1) + (N-i-1)*5.0/(N-1);
      cout << "Searching with grid spacing " << d << "..." << endl;
      S.setGridSpacing(d);
      S.search();
    }
  } else {
    S.search();
  }
  cout << "found " << S.numMatches() << " matches:" << endl;
  set<fasstSolution> matches = S.getMatches(); int i = 0;
  for (auto it = matches.begin(); it != matches.end(); ++it, ++i) {
    cout << *it << endl;
    Structure match = S.getMatchStructure(*it, true, FASST::matchType::WITHGAPS);
    match.writePDB("/tmp/match" + MstUtils::toString(i) + ".pdb");
  }
}
