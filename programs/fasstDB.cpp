#include "msttypes.h"
#include "mstoptions.h"
#include "mstfasst.h"
#include "mstrotlib.h"
#include <chrono>

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Creates a FASST database from input PDB files. Options:");
  op.addOption("p", "a file with a list of PDB files.", true);
  op.addOption("o", "output database file name.", true);
  op.addOption("m", "memory save flag (will store backbone only).");
  op.addOption("c", "clean up PDB files, so that only protein residues with enough of a backbone to support rotamer building survive.");
  op.addOption("pp", "store phi/psi properties in the database, if creating a new one from PDB files.");
  op.setOptions(argc, argv);

  FASST S;
  cout << "Reading structures..." << endl;
  S.setMemorySaveMode(op.isGiven("m"));
  vector<string> pdbFiles = MstUtils::fileToArray(op.getString("p"));
  for (int i = 0; i < pdbFiles.size(); i++) {
    Structure P(pdbFiles[i]);
    if (op.isGiven("c")) {
      Structure C; RotamerLibrary::extractProtein(C, P);
      if (P.residueSize() != C.residueSize()) {
        cout << pdbFiles[i] << ", had " << P.residueSize() << " residues, and " << C.residueSize() << " residues after cleaning..." << endl;
      }
      P = C;
    }
    S.addTarget(P);
    // compute and add some properties
    if (op.isGiven("pp")) {
      vector<Residue*> residues = P.getResidues();
      vector<mstreal> phi(residues.size()), psi(residues.size());
      for (int ri = 0; ri < residues.size(); ri++) {
        phi[ri] = residues[ri]->getPhi(false);
        psi[ri] = residues[ri]->getPsi(false);
      }
      S.addResidueProperties(S.numTargets() - 1, "phi", phi);
      S.addResidueProperties(S.numTargets() - 1, "psi", psi);
    }
  }
  S.writeDatabase(op.getString("o"));
}
