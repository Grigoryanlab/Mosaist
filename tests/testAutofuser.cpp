#include "msttypes.h"
#include "mstfuser.h"
#include "mstoptim.h"

using namespace std;
using namespace MST;

int main(int argc, char *argv[]) {
  if (argc < 3) {
    MstUtils::error("Takes a PDB file with overlapping segments, which must imply a single chain, and fuses into a single chain.\nUsage: ./testAutofuser [e.g., testfiles/abb-unfused.pdb] [output.pdb]", "main");
  }
  Structure unfused(argv[1]);
  string outFile(argv[2]);
  Structure fused = Fuser::autofuse(unfused.getResidues());
	fused.writePDB(outFile);
}
