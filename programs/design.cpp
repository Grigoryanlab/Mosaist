#include "msttypes.h"
#include "dtermen.h"
#include "mstoptions.h"
#include "mstsystem.h"
#include "mstrotlib.h"
#include <chrono>

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("dTERMen design. Options:");
  op.addOption("p", "template PDB file.", true);
  op.addOption("c", "dTERMen configuration file.", true);
  op.addOption("s", "selection of variable positions that should be designed. By default, all positions are designed.");
  op.addOption("ctx", "selection of fixed positions that should be treated as the fixed specificity context. If specified, will produce a specificity gap table.");
  op.addOption("o", "output base.", true);
  op.setOptions(argc, argv);

  Structure So(op.getString("p")), S;
  RotamerLibrary::extractProtein(S, So);
  vector<Residue*> variable, specContext;
  selector sel(S);
  if (op.isGiven("s")) {
    variable = sel.selectRes(op.getString("s"));
    cout << "selected " << variable.size() << " residues as variable according to selection '" << op.getString("s") << "':" << endl;
    for (int i = 0; i < variable.size(); i++) cout << "\t" << *(variable[i]) << endl;
  } else {
     variable = S.getResidues();
     cout << "selected all " << variable.size() << " residues as variable..." << endl;
  }
  if (op.isGiven("ctx")) {
    specContext = sel.selectRes(op.getString("ctx"));
    cout << "selected " << specContext.size() << " residues as variable according to selection '" << op.getString("ctx") << "':" << endl;
    for (int i = 0; i < specContext.size(); i++) cout << "\t" << *(specContext[i]) << endl;
    if (!MstUtils::setintersect(specContext, variable).empty()) MstUtils::error("specificity context selection cannot contain variable positions!");
  }
  EnergyTable E, specE;
  string etabFile = op.getString("o") + ".etab";
  string specEtabFile = op.getString("o") + ".spec.etab";

  // build the energy table
  if (!MstSys::fileExists(etabFile)) {
    dTERMen D(op.getString("c"));
    if (specContext.empty()) {
      E = D.buildEnergyTable(variable);
      E.writeToFile(etabFile);
    } else {
      E = D.buildEnergyTable(variable, vector<vector<string>>(), vector<vector<Residue*>>(), &specE, specContext);
      E.writeToFile(etabFile);
      specE.writeToFile(specEtabFile);
    }
  } else {
    cout << "reading previous energy table from " << etabFile << endl;
    E.readFromFile(etabFile);
    if (E.numSites() != variable.size()) MstUtils::error("pre-existing energy table has " + MstUtils::toString(E.numSites()) + " sites, while "  + MstUtils::toString(variable.size()) + " are selected for design");
  }

  vector<int> bestSol = E.mc(100, 1000000, 1.0, 0.01);
  mstreal lowE = E.scoreSolution(bestSol);
  cout << "lowest energy found is " << lowE << endl;
  Sequence bestSeq = E.solutionToSequence(bestSol);
  cout << "lowest-energy sequence: " << bestSeq.toString() << endl;
  cout << "mean energy is " << E.meanEnergy() << endl;
  cout << "estimated energy standard deviation is " << E.energyStdEst() << endl;

  // make a redesigned file, ready for repacking
  for (int i = 0; i < variable.size(); i++) {
    vector<Atom*> sidechain = MstUtils::setdiff(variable[i]->getAtoms(), RotamerLibrary::getBackbone(*(variable[i])));
    variable[i]->replaceAtoms(vector<Atom*>(), sidechain);
    variable[i]->setName(bestSeq.getResidue(i, true));
  }
  S.writePDB(op.getString("o") + ".red.pdb");
}
