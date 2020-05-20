#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <typeinfo>
#include <stdio.h>
#include <map>
#include <unistd.h>

#include "msttypes.h"
#include "mstoptions.h"
#include "mstrotlib.h"
#include "mstsequence.h"

using namespace std;
using namespace MST;

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Greedily clusters structures by backbone RMSD. Options:");
  op.addOption("l", "a file with a list of PDB files to read.", true);
  op.addOption("r", "RMSD cutoff to use for the clustering (greedy clustering is done).", true);
  op.addOption("oc", "optional: base name for outputing clusters as PDB files.");
  op.addOption("os", "optional: file name for outputing cluster sequences.");
  op.setOptions(argc, argv);
  double rmsdCut = op.getReal("r");
  string obasePDB = op.getString("oc", "");
  fstream of;

  // read PDB files and extract backbones as vectors of atom pointers
  vector<string> pdbFiles = MstUtils::fileToArray(op.getString("l"));
  vector<Structure> S(pdbFiles.size());
  vector<vector<Atom*>> backbones(S.size());
  for (int i = 0; i < pdbFiles.size(); i++) {
    S[i].readPDB(pdbFiles[i]);
    if (!RotamerLibrary::hasFullBackbone(S[i])) MstUtils::error("PDB file '" + pdbFiles[i] + "' does not have full backbone");
    backbones[i] = RotamerLibrary::getBackbone(S[i]);
  }
  pdbFiles.clear();

  // cluster
  Clusterer C;
  vector<vector<int>> clusters = C.greedyCluster(backbones, rmsdCut, 1000);

  // output clusters
  RMSDCalculator rc;
  for (int i = 0; i < clusters.size(); i++) {
    vector<int>& cluster = clusters[i];
    if (!obasePDB.empty()) MstUtils::openFile(of, obasePDB + MstUtils::toString(i+1) + ".pdb", fstream::out);
    for (int j = 0; j < cluster.size(); j++) {
      if (j == 0) cout << cluster[j] << ":";
      else cout << " " << cluster[j];
      if (!obasePDB.empty()) {
        Structure unit(S[cluster[j]]);
        if (j != 0) {
          rc.align(backbones[cluster[j]], backbones[cluster[0]], unit);
        }
        unit.writePDB(of);
      }
    }
    if (!obasePDB.empty()) of.close();
    cout << endl;
  }

  // output cluster sequences
  if (op.isGiven("os")) {
    MstUtils::openFile(of, op.getString("os"), fstream::out);
    for (int i = 0; i < clusters.size(); i++) {
      vector<int>& cluster = clusters[i];
      for (int j = 0; j < cluster.size(); j++) {
        of << Sequence(S[cluster[j]]).toString() << endl;
      }
      of << endl;
    }
    of.close();
  }
}
