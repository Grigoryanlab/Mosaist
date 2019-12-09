#include "msttypes.h"
#include "msttermanal.h"
#include "mstfasst.h"
#include "mstoptions.h"
#include <chrono>

using namespace std;
using namespace MST;

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Tests MST::TERMANAL. Options:");
  op.addOption("t", "path to MST test directory.", true);
  op.addOption("b", "path to a binary FASST database.", true);
  op.addOption("s", "test smoothing feature.", true);
  op.setOptions(argc, argv);

  MstUtils::setSignalHandlers();
  Structure motif(op.getString("t") + "/heptad.0388_0001.pdb");
  //Structure motif(op.getString("t") + "/heptad.0388_0014.pdb"); motif.getResidue(7).setNum(2);
  FASST F; F.readDatabase(op.getString("b"));
  TERMANAL T(&F);

  auto begin = chrono::high_resolution_clock::now();
  cout << "scoring..." << endl;
  vector<Residue*> central({&motif.getResidue(1)});
  //vector<Residue*> central({&motif.getResidue(7)});
  map<Residue*, vector<Residue*>> overlapSets;
  bool verbose = true;
  mstreal ss;
  if (op.getBool("s")) {
    Structure overlap(op.getString("t") + "/heptad.0388_0014.pdb"); overlap.getResidue(7).setNum(2);
    //Structure overlap(op.getString("t") + "/heptad.0388_0001.pdb");
    overlapSets = T.getOverlappingTerms(central, {&overlap});
    ss = T.structureScore(motif, central, overlapSets, verbose);
  } else ss = T.structureScore(motif, central, overlapSets, verbose);
  cout << "structure score = " << ss << endl;
  auto end = chrono::high_resolution_clock::now();
  cout << "scoring took " << chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " ms" << endl;
}
