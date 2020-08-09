#include "msttypes.h"
#include "dtermen.h"
#include "mstoptions.h"
#include "mstsystem.h"
#include "mstrotlib.h"
#include <chrono>

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Computes just pair energies using dTERMen design. Options:");
  op.addOption("p", "template PDB file.", true);
  op.addOption("c", "dTERMen configuration file.", true);
  op.addOption("sites", "a list of semicolon-separated MST selections specifying which site pairs to consider. If only one selection is given, then all residue pairs within it will be "
                        "computed. If two selections, then residue pairs between the two selections. If more than two, the number of selections should be even, in which case the pairs "
                        "considered will be between consecutive selection pairs (i.e., first with second, third with forth, etc.).", true);
  op.addOption("--new", "use the new style of pair-energy calculation, where homo-dimeric matches are accounted for in the statistics.");
  op.setOptions(argc, argv);

  Structure So(op.getString("p")), S;
  RotamerLibrary::extractProtein(S, So);
  vector<pair<Residue*, Residue*>> pairs;

  // identify pairs
  vector<string> selections = MstUtils::trim(MstUtils::split(op.getString("sites"), ";"));
  selector sel(S);
  if (selections.size() == 1) {
    vector<Residue*> residues = sel.selectRes(selections[0]);
    for (int i = 0; i < residues.size(); i++) {
      for (int j = i + 1; j < residues.size(); j++) {
        pairs.push_back(make_pair(residues[i], residues[j]));
      }
    }
  } else if (selections.size() % 2 == 0) {
    for (int k = 0; k < selections.size(); k += 2) {
      vector<Residue*> residuesI = sel.selectRes(selections[k]);
      vector<Residue*> residuesJ = sel.selectRes(selections[k + 1]);
      for (int i = 0; i < residuesI.size(); i++) {
        for (int j = 0; j < residuesJ.size(); j++) {
          pairs.push_back(make_pair(residuesI[i], residuesJ[j]));
        }
      }
    }
  } else {
    MstUtils::error("invalid number of selection specified in --sites");
  }

  dTERMen D(op.getString("c"));
  vector<res_t> alpha = D.getGlobalAlphabet();
  set<pair<Residue*, Residue*>> visited;
  vector<vector<mstreal>> pE;
  for (int i = 0; i < pairs.size(); i++) {
    if (visited.find(pairs[i]) != visited.end()) continue;
    if (op.isGiven("new")) pE = D.pairEnergiesNew(pairs[i].first, pairs[i].second);
    else pE = D.pairEnergies(pairs[i].first, pairs[i].second);
    for (int aai = 0; aai < pE.size(); aai++) {
      for (int aaj = 0; aaj < pE[aai].size(); aaj++) {
        cout << SeqTools::idxToTriple(alpha[aai]) << " " << SeqTools::idxToTriple(alpha[aaj]) << " " << pE[aai][aaj] << endl;
      }
    }
    visited.insert(pairs[i]);
  }
}
