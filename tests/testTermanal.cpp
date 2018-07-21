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
  op.setOptions(argc, argv);

  MstUtils::setSignalHandlers();
  Structure motif(op.getString("t") + "/heptad.0388_0001.pdb");
  FASST F; F.readDatabase(op.getString("b"));
  TERMANAL T(&F);
  auto begin = chrono::high_resolution_clock::now();
  cout << "scoring..." << endl;
  mstreal ss = T.structureScore(motif, {&(motif.getResidue(8))}, true);
  auto end = chrono::high_resolution_clock::now();
  cout << "scoring took " << chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " ms" << endl;
  cout << "structure score = " << ss << endl;
}
