#include "msttypes.h"
#include "dtermen.h"
#include "mstoptions.h"
#include <chrono>

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Tests the dTERMen class. Options:");
  op.addOption("t", "test files directory.", true);
  op.setOptions(argc, argv);

  Structure S(op.getString("t") + "/2ZTA.pdb");

  dTERMen D(op.getString("t") + "/dtermen.conf");
  cout << "-----" << endl;
  cout << "omega energy for ALA with omega 0 is: " << D.bbOmegaEner(0, "ALA") << endl;
  cout << "omega energy for ALA with omega 200 is: " << D.bbOmegaEner(200, "ALA") << endl;
  cout << "omega energy for ALA with omega -5 is: " << D.bbOmegaEner(-5, "ALA") << endl;
  cout << "omega energy for LEU with omega -180 is: " << D.bbOmegaEner(-180, "ALA") << endl;
  cout << "-----" << endl;
  cout << "env energy for GLU with freedom 0.2 is: " << D.envEner(0.2, "GLU") << endl;
  cout << "env energy for GLU with freedom 1.0 is: " << D.envEner(1.0, "GLU") << endl;
  cout << "-----" << endl;
  cout << "phi/psi energy for ALA in [-45, -65] is: " << D.bbPhiPsiEner(-45, -65, "ALA") << endl;
  cout << "phi/psi energy for GLY in [-45, -65] is: " << D.bbPhiPsiEner(-45, -65, "GLY") << endl;
  cout << "phi/psi energy for PRO in [-45, -65] is: " << D.bbPhiPsiEner(-45, -65, "PRO") << endl;
  cout << "phi/psi energy for LYS in [-45, -65] is: " << D.bbPhiPsiEner(-45, -65, "LYS") << endl;

  Residue* R = &(S.getResidue(MstUtils::randInt(S.residueSize())));
  vector<mstreal> selfE = D.selfEnergies(R);
  cout << "self energies at site " << *R << ":" << endl;
  for (int i = 0; i < D.globalAlphabetSize(); i++) {
    cout << D.indexToResName(i) << " => " << selfE[i] << endl;
  }
}
