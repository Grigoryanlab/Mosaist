#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <typeinfo>
#include <stdio.h>
#include <map>
#include <unistd.h>

#include "msttypes.h"
#include "mstoptions.h"
#include "mstfasst.h"

using namespace std;
using namespace MST;

void extractWindows(vector<vector<Atom*> >& windows, int w, const Structure& S) {
  const vector<string> bba = {"N", "CA", "C", "O"};
  for (int ci = 0; ci < S.chainSize(); ci++) {
    Chain& chain = S[ci];
    for (int ri = 0; ri < chain.residueSize() - w + 1; ri++) {
      vector<Atom*> win;
      for (int d = 0; d < w; d++) {
        Residue& res = chain[ri + d];
        for (int bi = 0; bi < bba.size(); bi++) {
          Atom* A = res.findAtom(bba[bi], false);
          if (A == NULL) break;
          win.push_back(A);
        }
      }
      if (win.size() != w * bba.size()) continue; // skip if any atoms missing
      windows.push_back(win);
    }
  }
}

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Clusters all local windows in given structures by backbone RMSD. Options:");
  op.addOption("w", "length of local window.", true);
  op.addOption("r", "RMSD cutoff to use for the clustering (greedy clustering is done).", true);
  op.addOption("i", "input PDB file.");
  op.addOption("db", "binary FASST database.");
  op.addOption("ob", "optional: output base name. If given, will write clusters as PDB files.");
  op.setOptions(argc, argv);

  int w = op.getInt("w");
  double rmsdCut = op.getReal("r");
  string obase = op.getString("ob", "");
  fstream of;
  vector<Structure*> S;

  vector<vector<Atom*> > windows;
  if (op.isGiven("i")) {
    S.push_back(new Structure(op.getString("i")));
    extractWindows(windows, w, *(S.back()));
  }

  if (op.isGiven("db")) {
    FASST F;
    F.readDatabase(op.getString("db"));
    for (int i = 0; i < F.numTargets(); i++) {
      S.push_back(new Structure(F.getTargetCopy(i)));
      extractWindows(windows, w, *(S.back()));
    }
  }
  cout << "extracted " << windows.size() << " windows. Clustering..." << endl;

  // cluster
  Clusterer C;
  vector<vector<int> > clusters = C.greedyCluster(windows, rmsdCut, 1000);

  // output clusters
  RMSDCalculator rc;
  for (int i = 0; i < clusters.size(); i++) {
    vector<int>& cluster = clusters[i];
    if (!obase.empty()) MstUtils::openFile(of, obase + MstUtils::toString(i+1) + ".pdb", fstream::out);
    for (int j = 0; j < cluster.size(); j++) {
      if (j == 0) cout << cluster[j] << ":";
      else cout << " " << cluster[j];
      if (!obase.empty()) {
        Structure unit(windows[cluster[j]]);
        if (j != 0) {
          rc.align(windows[cluster[j]], windows[cluster[0]], unit);
        }
        unit.writePDB(of);
      }
    }
    if (!obase.empty()) of.close();
    cout << endl;
  }
}
