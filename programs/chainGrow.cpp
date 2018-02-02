#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>
#include <typeinfo>
#include <stdio.h>

#include "msttypes.h"
#include "mstfasst.h"
#include "mstcondeg.h"
#include "mstfuser.h"
#include "mstcondeg.h"
#include "mstoptions.h"
#include "mstmagic.h"

using namespace std;
using namespace MST;

int main(int argc, char** argv) {
  MstOptions op;
  op.setTitle("Grows a chain from either the N- or the C-terminus, by the specified number of residues, aiming to make the overall structure as designable as possible. Options:");
  op.addOption("p", "starting PDB file.", true);
  op.addOption("L", "grow by this number of residues", true);
  op.addOption("n", "grow from the N-terminus.");
  op.addOption("c", "grow from the C-terminus.");
  op.addOption("rLib", "path to MST rotamer library.", true);
  op.addOption("o", "output base name.", true);
  op.addOption("b", "binary FASST database to use.", true);
  op.addOption("r", "grow randomly rather than deterministically. At each step, rather than picking the centroid of the top cluster, pick centroids randomly by cluster size.");
  op.setOptions(argc, argv);
  int pm = 2;
  mstreal cdCut = 0.01;
  int maxRetries = 10;
  int retry = 0;
  RMSDCalculator rc;
  Structure So(op.getString("p"));
  Structure S, Sp; RotamerLibrary::extractProtein(S, So);
  RotamerLibrary::standardizeBackboneNames(S);

  FASST F;
  F.setMemorySaveMode(true);
  F.readDatabase(op.getString("b"));
  F.pruneRedundancy(0.5);
  RotamerLibrary RL(op.getString("rLib"));
  Clusterer clusterGuru(false);
  fusionParams fuserOpts; fuserOpts.setNumIters(1000); fuserOpts.setVerbose(false);
  srand(time(NULL) + (int) getpid());

  // TERMify loop
  fstream out; MstUtils::openFile(out, op.getString("o") + ".grow.pdb", ios_base::out);
  out << "MODEL " << 0 << endl; S.writePDB(out); out << "ENDMDL" << endl;
  for (int nc = 0; nc < 2; nc++) {
    int c = 0;
    if ((nc == 0) && !op.isGiven("c")) continue;
    if ((nc == 1) && !op.isGiven("n")) continue;
    int N = op.getInt("L");
    for (int i = 0; i < N; i++) {
      cout << "growing in the " << (nc ? "N" : "C") << "-terminal direction, cycle " << i+1 << endl;
      ConFind C(&RL, S);

      // define TERM at the terminus
      cout << "selecting TERM..." << endl;
      vector<int> termResIndices;
      Residue& cres = nc ? S[0][0] : S.getChain(S.chainSize() - 1).getResidue(S[0].residueSize() - 1);
      Structure term = TERMUtils::selectTERM(cres, C, pm, cdCut, &termResIndices);
      int cresIdx = nc ? 0 : (term[0].residueSize() - 1); // where in the TERM the terminal residue is located
      cout << AtomPointerVector(term.getAtoms()) << endl;

      // find matches to it
      cout << "searching..." << endl;
      F.setQuery(term, false);
      F.setRMSDCutoff(RMSDCalculator::rmsdCutoff(term));
      // F.setMinNumMatches(1);
      F.setMaxNumMatches(1000);
      fasstSolutionSet matchList = F.search();
      cout << "\tfound " << matchList.size() << " matches" << endl;

      // cluster based on the location of the neighboring residue
      cout << "clustering..." << endl;
      vector<vector<Atom*> > nextRes; vector<int> nextResMatch;
      int del = nc ? -1 : 1; // relative location of "next" residue in the current direction
      for (int m = 0; m < matchList.size(); m++) {
        vector<int> resIndices = F.getMatchResidueIndices(matchList[m], FASST::matchType::REGION);
        Structure target = F.getMatchStructure(matchList[m], false, FASST::matchType::FULL, true);
        RotamerLibrary::standardizeBackboneNames(target);
        // does this match have a residue past the central residue?
        int nextIdx = resIndices[cresIdx] + del;
        if ((nextIdx >= target.residueSize()) || (nextIdx < 0) ||
            (target.getResidue(nextIdx).getChain() != target.getResidue(resIndices[cresIdx]).getChain())) continue;
        nextRes.push_back((vector<Atom*>) AtomPointerVector(target.getResidue(nextIdx).getAtoms()).clone());
        nextResMatch.push_back(m);
      }
      cout << "\tout of " << nextRes.size() << " cases" << endl;
      if (nextRes.size() == 0) {
        if ((c == 0) || (retry > maxRetries)) {
          cout << "\tno matches in first round, stopping current direction..." << endl;
          break;
        } else {
          cout << "\tbacking up one step to try a different random choice..." << endl;
          retry += 2; // we are about to take a step back, so will need two steps forward to advance
          cout << "\tbacking up one step to try a different random choice, try " << retry-1 << "..." << endl;
          i--; c--; S = Sp;
          continue;
        }
      }
      vector<vector<int>> clusters = clusterGuru.greedyCluster(nextRes, 0.3);
      cout << "\tresulted in " << clusters.size() << " clusters" << endl;
      for (int m = 0; m < nextRes.size(); m++) { AtomPointerVector(nextRes[m]).deletePointers(); } // clean up cloned fragments

      // take the centroid of top cluster and fuse it
      int pick = 0;
      if (op.isGiven("r")) {
        // pick a random cluster with probability proportional to size
        mstreal Ntot = 0, pc = 0.01;
        for (int i = 0; i < clusters.size(); i++) Ntot += (clusters[i].size() - 1 + pc);
        vector<mstreal> p(clusters.size()); // cumulative probabily
        for (int i = 0; i < clusters.size(); i++) {
          p[i] = (clusters[i].size() - 1 + pc)/Ntot;
          if (i > 0) p[i] += p[i-1];
        }
        mstreal r = MstUtils::randUnit();
        for (int i = 0; i < clusters.size(); i++) {
          if ((p[i] >= r) || (i == clusters.size() - 1)) { pick = i; break; }
        }
        cout << "\tpicking cluster " << pick << " of size " << clusters[pick].size() << endl;
      } else {
        cout << "\ttop cluster size " << clusters[0].size() << endl;
      }
      cout << "fusing..." << endl;
      int m = nextResMatch[clusters[pick][0]];
      vector<int> resIndices = F.getMatchResidueIndices(matchList[m], FASST::matchType::REGION);
      Structure match = F.getMatchStructure(matchList[m], false, FASST::matchType::FULL, true);
      vector<vector<Residue*> > resTopo(S.residueSize(), vector<Residue*>());
      for (int ri = 0; ri < S.residueSize(); ri++) resTopo[ri].push_back(&(S.getResidue(ri)));
      vector<bool> isFixed(resTopo.size(), true);
      for (int ri = 0; ri < resIndices.size(); ri++) {
        resTopo[termResIndices[ri]].push_back(&(match.getResidue(resIndices[ri])));
        for (int d = -pm; d <= pm; d++) {
          int idx = termResIndices[ri] + d;
          if ((idx >= 0) && (idx < S.residueSize())) isFixed[idx] = false;
        }
      }

      // extend the match by one residue
      int nextIdx = resIndices[cresIdx] + del;
      vector<Residue*> newPos(1, &(match.getResidue(nextIdx)));
      if (nc) {
        resTopo.insert(resTopo.begin(), newPos);
        isFixed.insert(isFixed.begin(), false);
      } else {
        resTopo.push_back(newPos);
        isFixed.push_back(false);
      }

      vector<int> fixed;
      for (int ri = 0; ri < isFixed.size(); ri++) {
        if (isFixed[ri]) fixed.push_back(ri);
      }
      Sp = S;
      S = Fuser::fuse(resTopo, fixed, fuserOpts);
      out << "MODEL " << c+1 << endl;
      S.writePDB(out); out << "ENDMDL" << endl;
      c++;
      if (retry > 0) retry--;
    }
  }
  out.close();
}
