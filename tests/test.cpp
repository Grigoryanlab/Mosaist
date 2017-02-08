#include "msttypes.h"

using namespace MST;

int main(int argc, char** argv) {

  if (argc < 3) {
    MstUtils::error("Usage: ./test [PDB file] [output file base]", "main");
  }
  string pdbFile(argv[1]);
  string outBase(argv[2]);

  // quick read and write
  printf("quick read/write test...\n");
  Structure S(pdbFile);
  S.writePDB(outBase + ".out.pdb");

  // memory test (sort of)
  int n = 100;
  printf("reading %d Structure objects into memory...\n", n);
  vector<Structure> many(n);
  for (int i = 0; i < n; i++) {
    many[i].readPDB(pdbFile);
  }
  printf("writing %d Structure objects to file...\n", n);
  for (int i = 0; i < n; i++) {
    many[i].writePDB(outBase + ".out.pdb");
  }
  printf("deleting %d Structure objects from memory...\n", n);
  many.clear();
  printf("TEST DONE\n");

  // phi/psi test
  for (int i = 0; i < S.chainSize(); i++) {
    Chain& chain = S[i];
    for (int j = 0; j < chain.residueSize(); j++) {
      Residue& res = chain[j];
      cout << res << ", phi = " << res.getPhi(false) << ", psi = " << res.getPsi(false) << endl;
    }
  }
}
