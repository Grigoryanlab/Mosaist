#include "msttypes.h"
#include "mstfuser.h"
#include "mstoptim.h"

using namespace std;
using namespace MST;

int main(int argc, char *argv[]) {
  if (argc < 3) {
    MstUtils::error("Usage: ./testFuser [testfiles/] [out-base]", "main");
  }
  string testdir(argv[1]), outBase(argv[2]);
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

  fusionParams opts; opts.setNumIters(1000); opts.setVerbose(false);
  Structure fused = Fuser::fuse(resTopo, fixed, opts);
  fused.writePDB(outBase + ".fused.pdb");


  // --- now try a different, more challenging case
  // assembly with three terms - in chains A, B, and C
	Structure s(testdir + "/fuserinput.pdb");
  Structure termA = Structure(s[0]).reassignChainsByConnectivity();
  Structure termB = Structure(s[1]).reassignChainsByConnectivity();
  Structure termC = Structure(s[2]).reassignChainsByConnectivity();

	// build residue topology
	vector<vector<Residue* > > residues(32);
	for (int i = 0; i < 10; i++) { // chain A of term1
		residues[i].push_back(&(termA[0][i]));
	}

	for (int i = 0; i < 10; i++) { // chain B of term1
		residues[i+11].push_back(&(termA[1][i]));
	}

	for (int i = 0; i < 10; i++) { // chain A of term2
		residues[i+11].push_back(&(termB[0][i]));
	}

	for (int i = 0; i < 10; i++) { // chain B of term2
		residues[i+22].push_back(&(termB[1][i]));
	}

	for (int i = 0; i < 5; i++) { // chain A of term3
		residues[i+21].push_back(&(termC[0][i]));
	}

	for (int i = 0; i < 5; i++) { // chain B of term3
		residues[i+10].push_back(&(termC[1][i]));
	}

	// print out residues
	for (int i = 0; i < residues.size(); i++) {
		cout << i;
		for (int j = 0; j < residues[i].size(); j++) {
			cout << " " << *residues[i][j];
		}
		cout << endl;
	}

  opts.setNumCycles(2);
	Structure fusedStruct = Fuser::fuse(residues, vector<int>(1, 0), opts);
	fusedStruct.writePDB(outBase + ".fused1.pdb");
}
