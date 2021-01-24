#include "msttypes.h"
#include "dtermen.h"
#include "mstoptions.h"
#include "mstsystem.h"
#include "mstrotlib.h"
#include <chrono>

void getSelections(const string& selStr, vector<string>& selections) {
  vector<string> strSelections = MstUtils::trim(MstUtils::split(selStr, ";"));
  selections.insert(selections.end(), strSelections.begin(), strSelections.end());
}

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Computes just pair energies using dTERMen design. Options:");
  op.addOption("p", "template PDB file.", true);
  op.addOption("c", "dTERMen configuration file.", true);
  op.addOption("sites", "a list of semicolon-separated MST selections specifying which site pairs to consider. If only one selection is given, then all residue pairs within it will be "
                        "computed. If two selections, then residue pairs between the two selections. If more than two, the number of selections should be even, in which case the pairs "
                        "considered will be between consecutive selection pairs (i.e., first with second, third with forth, etc.).");
  op.addOption("sitesL", "a file containing a list of sites formatted like in --sites.");
  op.addOption("o", "output file (prints to STDOUT if not specified).");
  op.addOption("new", "use the new style of pair-energy calculation, where homo-dimeric matches are accounted for in the statistics.");
  op.addOption("new2", "use the new style of pair-energy calculation, where homo-dimeric matches are accounted for in the statistics.");
  op.setOptions(argc, argv);
  if (!(op.isGiven("sites") || op.isGiven("sitesL"))) MstUtils::error("either --sites or --sitesL must be specified!");

  Structure So(op.getString("p")), S;
  RotamerLibrary::extractProtein(S, So);
  vector<pair<Residue*, Residue*>> pairs;

  // identify pairs
  vector<string> selections;
  if (op.isGiven("sites")) getSelections(op.getString("sites"), selections);
  if (op.isGiven("sitesL")) getSelections(MstUtils::join("\n", MstUtils::fileToArray(op.getString("sitesL"))), selections);
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
  ofstream outFS;
  ostream* out;
  if (op.isGiven("o")) {
    string outFile = op.getString("o");
    string outDir = MstSys::splitPath(outFile, 0);
    MstSys::cmkdir(outDir, true);
    outFS.open(outFile);
    out = &outFS;
    cout << "Writing to " << outFile + " ..." << endl;
  } else {
    out = &cout;
    cout << "Writing to STDOUT ..." << endl;
  }
  for (int i = 0; i < pairs.size(); i++) {
    if (visited.find(pairs[i]) != visited.end()) continue;
    if (op.isGiven("new")) pE = D.pairEnergiesNew(pairs[i].first, pairs[i].second);
    else if (op.isGiven("new2")) pE = D.pairEnergiesNew2(pairs[i].first, pairs[i].second);
    else pE = D.pairEnergies(pairs[i].first, pairs[i].second);
    for (int aai = 0; aai < pE.size(); aai++) {
      for (int aaj = 0; aaj < pE[aai].size(); aaj++) {
        string res1 = pairs[i].first->getChainID() + "," + to_string(pairs[i].first->getNum()), res2 = pairs[i].second->getChainID() + "," + to_string(pairs[i].second->getNum());
        *out << res1 << " " << res2 << " " << SeqTools::idxToTriple(alpha[aai]) << " " << SeqTools::idxToTriple(alpha[aaj]) << " " << pE[aai][aaj] << endl;
      }
    }
    visited.insert(pairs[i]);
  }
  cout << "Done" << endl;
}
