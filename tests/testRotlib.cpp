#include "msttypes.h"
#include "msttransforms.h"
#include "mstrotlib.h"

using namespace MST;

int main(int argc, char** argv) {

  if (argc < 4) {
    MstUtils::error("Usage: ./test [PDB file] [rotlib.bin file] [output file base]", "main");
  }
  string pdbFile(argv[1]);
  string rotLibFile(argv[2]);
  string outBase(argv[3]);

  // read structure
  Structure S(pdbFile);

  // read rotamer library
  RotamerLibrary R(rotLibFile);
  vector<string> aas = R.availableAminoAcids();

  // place some rotamers
  for (int i = 0; i < 10; i++) {
    Chain& chain = S[MstUtils::randInt(S.chainSize())];
    Residue& res = chain[MstUtils::randInt(chain.residueSize())];
    string aa = aas[MstUtils::randInt(aas.size())];
    int ri = MstUtils::randInt(R.numberOfRotamers(aa, res.getPhi(), res.getPsi()));
    cout << "At position " << res << " placing rotamer " << ri << " of " << aa << endl;
    R.placeRotamer(res, aa, ri);
  }

  printf("TEST DONE\n");
}
