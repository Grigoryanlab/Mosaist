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
  printf("placing a bunch of rotamers randomly...\n");
  for (int i = 0; i < 100; i++) {
    Chain& chain = S[MstUtils::randInt(S.chainSize())];
    Residue& res = chain[MstUtils::randInt(chain.residueSize())];
    string aa = aas[MstUtils::randInt(aas.size())];
    int ri = MstUtils::randInt(R.numberOfRotamers(aa, res.getPhi(), res.getPsi()));
    mstreal prob = R.rotamerProbability(aa, ri, res.getPhi(), res.getPsi());
    cout << "At position " << res << " placing rotamer " << ri << " of " << aa << ", with probability " << prob << endl;
    R.placeRotamer(res, aa, ri);
  }
  S.writePDB(outBase + ".pdb", "renumber");

  // place every rotamer at every position (as a test)
  printf("now placing every rotamer at every position...\n");
  for (int i = 0; i < S.chainSize(); i++) {
    Chain& chain = S[i];
    for (int j = 0; j < chain.residueSize(); j++) {
      Residue& res = chain[j];
      for (int k = 0; k < aas.size(); k++) {
        string aa = aas[k];
        int nr = R.numberOfRotamers(aa, res.getPhi(), res.getPsi());
        for (int r = 0; r < nr; r++) {
          R.placeRotamer(res, aa, r);
        }
      }
    }
  }
  printf("done!\n");

  printf("TEST DONE\n");
}
