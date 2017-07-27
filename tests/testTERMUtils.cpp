#include "msttypes.h"
#include "mstmagic.h"
#include <algorithm>

using namespace MST;

int main(int argc, char** argv) {
  mstreal pad = 100;    // search for matching segments in a pad x pad x pad box around the center of the central unit
  int N = 1000;      // examine top this many most promissing points to expand into regions
  int n = 20;        // write this many best-scoring regions
  mstreal dcut = 2.0;   // inter-atomic distance cutoff for counting neighbors
  int ri;

  if (argc < 2) {
    MstUtils::error("Usage: ./test [file with a list of PDB file paths]", "main");
  }
  string pdbListFile(argv[1]);
  vector<string> pdbFiles = MstUtils::fileToArray(pdbListFile);
  if (pdbFiles.size() <= 1) MstUtils::error("expected at least two PDB files in list", "main");

  // read in the central TERM and all the attached ones
  Structure C(pdbFiles[0]);
  vector<Structure*> TERMs(pdbFiles.size());
  for (int k = 1; k < pdbFiles.size(); k++) {
    TERMs[k] = new Structure(pdbFiles[k]);
  }

  // do the magic
  printf("starting to look for the most designable fragments...\n");
  TERMUtils::mostDesignableFragments(C, TERMs, 10, NULL, NULL, "region");
  printf("done looking for the most designable fragments...\n");

  printf("done!\n");

  // free
  for (int k = 0; k < TERMs.size(); k++) delete(TERMs[k]);
}
