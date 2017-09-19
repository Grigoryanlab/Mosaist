#include "msttypes.h"
#include "mstsystem.h"

using namespace MST;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    MstUtils::error("Renumbers residues and (if necessary) renames chains in input PDB file(s). Outputs with suffix _renum.", "main");
  }
  for (int i = 1; i < argc; i++) {
    Structure S(argv[i]); // this already renames chains as necessary
    for (int ci = 0; ci < S.chainSize(); ci++) {
      Chain& C = S[ci];
      for (int ri = 0; ri < C.residueSize(); ri++) {
        C[ri].setNum(ri+1);
      }
    }
    S.writePDB(MstSys::pathBase(argv[i]) + "_renum.pdb");
  }
}
