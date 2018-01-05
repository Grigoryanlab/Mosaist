#include "msttypes.h"
#include "mstrotlib.h"
#include "mstcondeg.h"

using namespace MST;

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Compute the TERM-induced amino-acid substitution matrix. Options:");
  op.addOption("d", "a FASST database.", true);
  op.setOptions(argc, argv);
  RMSDCalculator rc;
  FASST F;

  F.setMemorySaveMode(true);
  F.readDatabase(op.getString());

  while (true) {
    int ti = MstUtils::randInt(0, F.numTargets() - 1);
    Structure S(F.getTarget(ti));
    
  }
}
