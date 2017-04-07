#include "msttypes.h"
#include "mstfuser.h"
#include "mstoptim.h"

using namespace std;
using namespace MST;

int main(int argc, char *argv[]) {
  if (argc < 3) {
    MstUtils::error("Usage: ./testFuser [testfiles/] [out.pdb]", "main");
  }
  string testdir(argv[1]), opdbf(argv[2]);
  Structure unitA(testdir + "heptad.0388_0001.pdb");
  Structure unitB(testdir + "heptad.0388_0007.pdb");
  Structure bridge(testdir + "heptad.0388_0014.pdb");

  // will fuse chain B of unitA with chain A of unitB
  Chain& chainA = unitA[1];
  Chain& chainB = unitB[0];
  Chain& chainBridge = bridge[0];
  int overlapN = 2, overlapC = 2; // overlap lengths on the N- and C-termini of the bridge
  int L = chainA.residueSize() + chainB.residueSize() + chainBridge.residueSize() - overlapN - overlapC; // length of fused structure
  vector<vector<Residue*> > resTopo(L);

  // first chain
  for (int i = 0; i < chainA.residueSize(); i++) {
    resTopo[i].push_back(&(chainA[i]));
  }

  // bridge
  for (int i = 0; i < chainBridge.residueSize(); i++) {
    resTopo[i + chainA.residueSize() - overlapN].push_back(&(chainBridge[i]));
  }

  // second chain
  for (int i = 0; i < chainB.residueSize(); i++) {
    resTopo[i + chainA.residueSize() + chainBridge.residueSize() - overlapN - overlapC].push_back(&(chainB[i]));
  }

  // fix just the parts of unit chains that are not overlapped by the bridge
  vector<int> fixed;
  for (int i = 0; i < chainA.residueSize() - overlapN; i++) fixed.push_back(i);
  for (int i = L - 1; i >= L - (chainB.residueSize() - overlapC); i--) fixed.push_back(i);
  cout << "Will fix residues: " << MstUtils::vecToString(fixed) << endl;
  cout << "Leaving " << L - fixed.size() << " mobile" << endl;

  // Structure fused = Fuser::fuse(resTopo, fixed, 10, true);
  Structure fused = Fuser::fuse(resTopo, fixed, 10, true);
  fused.writePDB(opdbf);
}
