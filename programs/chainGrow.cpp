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
  op.addOption("o", "output base name.", true);
  op.addOption("b", "binary FASST database to use.", true);
  op.setOptions(argc, argv);
  int pm = 2;
  mstreal cdCut = 0.01;
  RMSDCalculator rc;
  Structure S(op.getString("p"));
  Structure So = S; // RotamerLibrary::extractProtein(S, So);
  FASST F;
  F.setMemorySaveMode(true);
  F.readDatabase(op.getString("b"));
  F.pruneRedundancy(0.5);
  RotamerLibrary RL(op.getString("rLib"));
  Clusterer clusterGuru(false);
  fusionParams fuserOpts; fuserOpts.setNumIters(1000); fuserOpts.setVerbose(false);

  // TERMify loop
  for (int nc = 0; nc < 2; nc++) {
    if ((nc == 0) && !op.isGiven("n")) continue;
    if ((nc == 1) && !op.isGiven("c")) continue;
    int N = op.getInt("N");
    for (int i = 0; i < N; i++) {
      ConFind C(&RL, S);

      // define TERM at the terminus
      vector<int> termResIndices;
      Residue& cres = nc ? S[0][0] : S.getChain(S.chainSize() - 1).getResidue(S[0].residueSize() - 1);
      Structure term = TERMUtils::selectTERM(cres, C, pm, cdCut, &termResIndices);
      int cresIdx = nc ? 0 : (term[0].residueSize() - 1); // where in the TERM the terminal residue is located

      // find matches to it
      F.setQuery(term, false);
      F.setRMSDCutoff(RMSDCalculator::rmsdCutoff(term));
      F.setMaxNumMatches(1000);
      fasstSolutionSet matchList = F.search();
      if (matchList.size() == 0) {
        cout << "\tSTOP: could not find any more matches at the terminus, stopping current direction..." << endl;
        break;
      }

      // cluster based on the location of the neighboring residue
      vector<vector<Atom*> > nextRes; vector<int> nextResMatch;
      int del = nc ? -1 : 1; // relative location of "next" residue in the current direction
      for (int m = 0; m < matchList.size(); m++) {
        vector<int> resIndices = F.getMatchResidueIndices(matchList[m], FASST::matchType::REGION);
        Structure target = F.getMatchStructure(matchList[m], false, FASST::matchType::FULL, true);
        // does this match have a residue past the central residue?
        int nextIdx = resIndices[cresIdx] + del;
        if ((nextIdx >= target.residueSize()) || (nextIdx < 0)) continue;
        nextRes.push_back(target.getResidue(nextIdx).getAtoms());
        nextResMatch.push_back(m);
      }
      vector<vector<int>> clusters = clusterGuru.greedyCluster(nextRes, 0.1);

      // take the centroid of top cluster and fuse it
      // TODO: enable a Markov-chain walk by randomizing this choice!
      int m = nextResMatch[clusters[0][0]];
      vector<int> resIndices = F.getMatchResidueIndices(matchList[m], FASST::matchType::REGION);
      Structure match = F.getMatchStructure(matchList[m], false, FASST::matchType::REGION, true);
      vector<vector<Residue*> > resTopo(S.residueSize(), vector<Residue*>());
      for (int ri = 0; ri < S.residueSize(); ri++) resTopo[ri].push_back(&(S.getResidue(ri)));
      vector<bool> isFixed(S.residueSize(), true);
      for (int ri = 0; ri < match.residueSize(); ri++) {
        resTopo[termResIndices[ri]].push_back(&(match.getResidue(ri)));
        for (int d = -pm; d <= pm; d++) {
          int idx = termResIndices[ri] + d;
          if ((idx >= 0) && (idx < S.residueSize())) isFixed[idx] = false;
        }
      }
      vector<int> fixed;
      for (int ri = 0; ri < isFixed.size(); ri++) {
        if (isFixed[ri]) fixed.push_back(ri);
      }
      S = Fuser::fuse(resTopo, fixed, fuserOpts);
    }
  }
}
