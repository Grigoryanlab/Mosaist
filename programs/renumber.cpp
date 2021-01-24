#include "msttypes.h"
#include "mstsystem.h"
#include "mstoptions.h"

using namespace MST;

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Renumber/clean a PDB file. Options:");
  op.addOption("p", "input PDB file.", true);
  op.addOption("s", "split chains by connectivity.");
  op.addOption("nr", "do not actually do the renumbering.");
  op.addOption("o", "output PDB file name. If not specified, suffix _renum will be added to the input name.");
  op.setOptions(argc, argv);

  Structure S(op.getString("p"));
  if (op.isGiven("s")) S = S.reassignChainsByConnectivity();
  if (!op.isGiven("nr")) {
    for (int ci = 0; ci < S.chainSize(); ci++) {
      Chain& C = S[ci];
      for (int ri = 0; ri < C.residueSize(); ri++) {
        C[ri].setNum(ri+1);
      }
    }
  }
  if (op.isGiven("o")) S.writePDB(op.getString("o"));
  else S.writePDB(MstSys::pathBase(op.getString("p")) + "_renum.pdb");
}
