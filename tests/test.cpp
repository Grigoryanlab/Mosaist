#include "msttypes.h"
#include "mstsystem.h"

using namespace MST;

int main(int argc, char** argv) {

  if (argc < 3) {
    MstUtils::error("Usage: ./test [PDB file] [output file base]", "main");
  }
  string pdbFile(argv[1]);
  string outBase(argv[2]);

  // cout << "getting shared lock...\n" << MstSys::getNetLock("gevorg", true) << endl;
  // cout << "working..." << endl; sleep(10);
  // cout << "releasing the lock...\n" << MstSys::releaseNetLock("gevorg") << endl;
  // cout << "getting exlusive lock...\n" << MstSys::getNetLock("gevorg") << endl;
  // cout << "working..." << endl; sleep(10);
  // cout << "releasing the lock...\n" << MstSys::releaseNetLock("gevorg") << endl;
  // exit(0);

  // quick read and write
  printf("quick read/write test...\n");
  Structure S(pdbFile);
  S.writePDB(outBase + ".out.pdb");

  // reassign by connectivity
  S = S.reassignChainsByConnectivity();
  S.writePDB(outBase + ".conn.out.pdb");

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

  // chains test
  Structure s1(pdbFile);
  Chain* c = s1.appendChain("Z");
  Structure s2(pdbFile);
  for (int i = 0; i < s2.residueSize(); i++) {
      c->appendResidue(new Residue(s2.getResidue(i)));
  }
  s1.writePDB(outBase + ".ch1.out.pdb");
  s1.deleteChain(c);
  s1.writePDB(outBase + ".ch2.out.pdb");

  // write to binary data test
  cout << "writing to/reading from binary file..." << endl;
  Structure M(pdbFile);
  M.writeData(outBase + ".bin");
  Structure Mread;
  Mread.readData(outBase + ".bin");
  Mread.writePDB(outBase + ".bin.pdb");
}
