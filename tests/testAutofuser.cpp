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
  for (int i = 1; i < argc-1; i++) {
    Structure unfused(argv[i]);
    vector<Residue*> res = unfused.getResidues();
    for (int j = 0; j < unfused.residueSize(); j++) {
      toFuse.push_back(new Residue(*res[j]));
    }
  }
  Structure fused = Fuser::autofuse(toFuse, 1);
	fused.writePDB(outFile);
  for (int j = 0; j < toFuse.size(); j++) delete(toFuse[j]);
}
