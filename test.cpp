#include "mstlib.h"

int main(int argc, char** argv) {
  if (argc < 3) {
    MstUtils::error("Usage: ./test [PDB file] [output file base]", "main");
  }
  string pdbFile(argv[1]);
  string outBase(argv[2]);
  Structure S(pdbFile);
  S.writePDB(outBase + ".out.pdb");
}
