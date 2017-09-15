#include "msttypes.h"
#include "mstfuser.h"
#include "mstoptim.h"

using namespace std;
using namespace MST;

int main(int argc, char *argv[]) {
  if (argc < 3) {
    MstUtils::error("Takes a PDB file with overlapping segments, which must imply a single chain, and fuses into a single chain. Usage:\n   ./testAutofuser [e.g., testfiles/abb-unfused.pdb] [output.pdb]\n   ./testAutofuser [e.g., pieces*.pdb] [output.pdb]", "main");
  }
  string outFile(argv[argc-1]);
  vector<Residue*> toFuse;
  vector<Structure> inputs;
  for (int i = 1; i < argc-1; i++) {
    inputs.push_back(Structure(argv[i]));
    vector<Residue*> res = inputs.back().getResidues();
    toFuse.insert(toFuse.end(), res.begin(), res.end());
  }
  Structure fused = Fuser::autofuse(toFuse, 1, 100, 1, true);
	fused.writePDB(outFile);
}
