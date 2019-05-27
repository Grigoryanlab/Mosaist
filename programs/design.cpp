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
  op.addOption("s", "selection of variable positions that should be designed. By default, all positions are designed. In symmetric design instances, this selection should be from the central unit cell.");
  op.addOption("ctx", "selection of fixed positions that should be treated as the fixed specificity context. If specified, will produce a specificity gap table.");
  op.addOption("sym", "design with crystallographic symmetry. Must specify a semicolon-separated list of unit cells, with the "
                      "first one being the central unit cell that design will focus on (it must have all surrounding unit cells "
                      "present). Each unit cell must be a space-separated list of chains. E.g., --sym 'A B; C D; E F' means that "
                      "the unit cell has two chains, with the central unit cell composed of chains A and B, and two other unit"
                      "cells (one composed of chains C and D and another one of chains E and F) also present in the structure.");
  op.addOption("seq", "skip design and simply put on this sequence, dumping the resulting PDB file.");
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

  // parse crystallographic symmetry information
  vector<vector<Residue*>> images;
  if (op.isGiven("sym")) {
    vector<string> unitCells = MstUtils::split(op.getString("sym"), ";");
    if (unitCells.size() < 2) MstUtils::error("at least two unit cells must be specified if --sym is used");
    images.resize(unitCells.size());
    for (int i = 0; i < unitCells.size(); i++) {
      vector<string> uc = MstUtils::split(MstUtils::trim(unitCells[i]), " ");
      if (uc.empty()) MstUtils::error("one of the unit cells is empty in symmetry string " + op.getString("sym"));
      for (int j = 0; j < uc.size(); j++) {
        Chain* C = S.getChainByID(uc[j]);
        if (C == NULL) MstUtils::error("chain '" + uc[i] + "' not found, specified in symmetry string " + op.getString("sym"));
        vector<Residue*> chainResis = C->getResidues();
        images[i].insert(images[i].end(), chainResis.begin(), chainResis.end());
      }
    }
  }

  // build the energy table
  Sequence bestSeq;
  if (op.isGiven("seq")) {
    bestSeq = Sequence(op.getString("seq"));
  } else {
    if (!MstSys::fileExists(etabFile)) {
      dTERMen D(op.getString("c"));
      if (specContext.empty()) {
        E = D.buildEnergyTable(variable, vector<vector<string>>(), images);
        E.writeToFile(etabFile);
      } else {
        E = D.buildEnergyTable(variable, vector<vector<string>>(), images, &specE, specContext);
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
    bestSeq = E.solutionToSequence(bestSol);
    cout << "lowest-energy sequence: " << bestSeq.toString() << endl;
    cout << "mean energy is " << E.meanEnergy() << endl;
    cout << "estimated energy standard deviation is " << E.energyStdEst() << endl;
  }

  // make a redesigned file, ready for repacking
  if (!op.isGiven("sym")) images.push_back(variable);
  for (int i = 0; i < variable.size(); i++) {
    // the variable selection must be from the first image
    int idx = find(images[0].begin(), images[0].end(), variable[i]) - images[0].begin();
    for (int j = 0; j < images.size(); j++) {
      Residue* varRes = images[j][idx];
      vector<Atom*> sidechain = MstUtils::setdiff(varRes->getAtoms(), RotamerLibrary::getBackbone(*varRes));
      varRes->replaceAtoms(vector<Atom*>(), sidechain);
      varRes->setName(bestSeq.getResidue(i, true));
    }
  }
  S.writePDB(op.getString("o") + ".red.pdb");
}
