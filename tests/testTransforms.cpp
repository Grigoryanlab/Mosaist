#include "msttypes.h"
#include "msttransforms.h"

using namespace MST;

int main(int argc, char** argv) {

  if (argc < 3) {
    MstUtils::error("Usage: ./test [PDB file] [output file base]", "main");
  }
  string pdbFile(argv[1]);
  string outBase(argv[2]);

  // read structure
  Structure S(pdbFile);

  // do some manipulations to it
  Transform T = TransformFactory::rotateAroundX(90);
  T.apply(S);
  S.writePDB(outBase + ".X90.pdb");
  T.apply(S);
  S.writePDB(outBase + ".X180.pdb");
  T.apply(S);
  S.writePDB(outBase + ".X270.pdb");
  T.apply(S);
  S.writePDB(outBase + ".X360.pdb");

  printf("TEST DONE\n");
}
