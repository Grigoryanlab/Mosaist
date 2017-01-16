#include "msttypes.h"

using namespace MST;

int main(int argc, char** argv) {

  double angle;
  double base = 360.0;
  angle = -12.5; printf("fmod(%f, %f) = %f\n", angle, base, fmod(angle, base));
  angle = 12.5; printf("fmod(%f, %f) = %f\n", angle, base, fmod(angle, base));
  angle = -180.3; printf("fmod(%f, %f) = %f\n", angle, base, fmod(angle, base));
  angle = 180.3; printf("fmod(%f, %f) = %f\n", angle, base, fmod(angle, base));
  angle = -365.2; printf("fmod(%f, %f) = %f\n", angle, base, fmod(angle, base));
  angle = 365.2; printf("fmod(%f, %f) = %f\n", angle, base, fmod(angle, base));

  angle = -12.5; printf("mod(%f, %f) = %f\n", angle, base, MstUtils::mod(angle, base));
  angle = 12.5; printf("mod(%f, %f) = %f\n", angle, base, MstUtils::mod(angle, base));
  angle = -180.3; printf("mod(%f, %f) = %f\n", angle, base, MstUtils::mod(angle, base));
  angle = 180.3; printf("mod(%f, %f) = %f\n", angle, base, MstUtils::mod(angle, base));
  angle = -365.2; printf("mod(%f, %f) = %f\n", angle, base, MstUtils::mod(angle, base));
  angle = 365.2; printf("mod(%f, %f) = %f\n", angle, base, MstUtils::mod(angle, base));

  if (argc < 3) {
    MstUtils::error("Usage: ./test [PDB file] [output file base]", "main");
  }
  string pdbFile(argv[1]);
  string outBase(argv[2]);

  // quick read and write
  printf("quick read/write test...\n");
  Structure S(pdbFile);
  S.writePDB(outBase + ".out.pdb");

  // memory test
  int n = 1000;
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
}
