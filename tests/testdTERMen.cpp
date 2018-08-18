#include "msttypes.h"
#include "dtermen.h"
#include "mstoptions.h"
#include <chrono>

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Tests the dTERMen class. Options:");
  op.addOption("p", "template PDB file.", true);
  op.addOption("t", "test files directory.", true);
  op.addOption("o", "output base.", true);
  op.setOptions(argc, argv);

  Structure S(op.getString("p"));
  dTERMen D(op.getString("t") + "/dtermen.conf");

  EnergyTable E = D.buildEnergyTable(S.getResidues());
  vector<int> bestSol = E.mc(100, 1000000, 1.0, 0.01);
  mstreal lowE = E.scoreSolution(bestSol);
  cout << "lowest energy found is " << lowE << endl;
  cout << "lowest-energy sequence: " << (E.solutionToSequence(bestSol)).toString() << endl;
  cout << "mean energy is " << E.meanEnergy() << endl;
  cout << "estimated energy standard deviation is " << E.energyStdEst() << endl;
  E.writeToFile(op.getString("o") + ".etab");
}

// cout << "-----" << endl;
// cout << "omega energy for ALA with omega 0 is: " << D.bbOmegaEner(0, "ALA") << endl;
// cout << "omega energy for ALA with omega 200 is: " << D.bbOmegaEner(200, "ALA") << endl;
// cout << "omega energy for ALA with omega -5 is: " << D.bbOmegaEner(-5, "ALA") << endl;
// cout << "omega energy for LEU with omega -180 is: " << D.bbOmegaEner(-180, "ALA") << endl;
// cout << "-----" << endl;
// cout << "env energy for GLU with freedom 0.2 is: " << D.envEner(0.2, "GLU") << endl;
// cout << "env energy for GLU with freedom 1.0 is: " << D.envEner(1.0, "GLU") << endl;
// cout << "-----" << endl;
// cout << "phi/psi energy for ALA in [-45, -65] is: " << D.bbPhiPsiEner(-45, -65, "ALA") << endl;
// cout << "phi/psi energy for GLY in [-45, -65] is: " << D.bbPhiPsiEner(-45, -65, "GLY") << endl;
// cout << "phi/psi energy for PRO in [-45, -65] is: " << D.bbPhiPsiEner(-45, -65, "PRO") << endl;
// cout << "phi/psi energy for LYS in [-45, -65] is: " << D.bbPhiPsiEner(-45, -65, "LYS") << endl;
//
// ConFind C(D.getRotamerLibrary(), S);
// Residue* Ri = &(S.getResidue(MstUtils::randInt(S.residueSize())));
// vector<Residue*> cRes = C.getContactingResidues(Ri, 0.01);
// Residue* Rj = cRes[MstUtils::randInt(cRes.size())];
// int dt = std::chrono::time_point_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
// vector<vector<mstreal> > pairE = D.pairEnergies(Ri, Rj);
// dt = std::chrono::time_point_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count() - dt;
// cout << "computing pair energies took: " << dt << " seconds" << endl;
// cout << "pair energies for " << *Ri << " x " << *Rj << " are:" << endl;
// for (int i = 0; i < D.globalAlphabetSize(); i++) cout << " " << D.indexToResName(i);
// cout << endl;
// for (int i = 0; i < D.globalAlphabetSize(); i++) {
//   for (int j = 0; j < D.globalAlphabetSize(); j++) {
//     cout << pairE[i][j] << " ";
//   }
//   cout << endl;
// }
//
// // test pair energies
// Residue* R = &(S.getResidue(MstUtils::randInt(S.residueSize())));
// dt = std::chrono::time_point_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
// vector<mstreal> selfE = D.selfEnergies(R, true);
// dt = std::chrono::time_point_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count() - dt;
// cout << "computing self energies took: " << dt << " seconds" << endl;
// cout << "self energies at site " << *R << ":" << endl;
// for (int i = 0; i < D.globalAlphabetSize(); i++) {
//   cout << D.indexToResName(i) << " => " << selfE[i] << endl;
// }
