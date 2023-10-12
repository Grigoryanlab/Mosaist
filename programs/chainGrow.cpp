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
#include "mstlinalg.h"

using namespace std;
using namespace MST;

/* Calculates the cosine between local chain tangent vectors, as a function of
 * residue separation. If the chain is too short, returns an empty vector. The
 * given vector of Residue pointers is assumed to be in N-to-C order. All atoms
 * within the residues are considered (i.e., no filtering for backbone is done). */
vector<mstreal> chainAutocorr(const vector<Residue*> chain) {
  int winL = 5; // window length for getting the local tange
  if (chain.size() < winL + 1) return vector<mstreal>();

  vector<CartesianPoint> pAxes(chain.size() - winL + 1);
  for (int i = 0; i < pAxes.size(); i++) {
    AtomPointerVector window;
    for (int k = i; k < i + winL; k++) {
      window.push_back(chain[k]);
    }
    Matrix axes = MstLinAlg::getPrincipalAxes(window);
    pAxes[i] = CartesianPoint((vector<mstreal>) axes.column(0)).getUnit();
    CartesianPoint ncVec = CartesianPoint(window.back()) - CartesianPoint(window.front());
    if (pAxes[i].dot(ncVec) < 0) pAxes[i] = -pAxes[i]; // make sure points in the N-to-C direction
  }
  vector<mstreal> cosines(chain.size() - winL, 0);
  vector<int> counts(chain.size() - winL, 0);
  for (int i = 0; i < pAxes.size(); i++) {
    for (int j = i+1; j < pAxes.size(); j++) {
      cosines[j-i-1] += pAxes[i].dot(pAxes[j]);
      counts[j-i-1]++;
    }
  }
  for (int i = 0; i < cosines.size(); i++) cosines[i] /= counts[i];
  return cosines;
}

/* Picks a suitable match for extending the structure S from within the set of
 * FASST solutions matchList. F is the FASST object used to locate the matches.
 * endIdx is the index of the residue (in the searched TERM) from which the
 * chain is supposed to grow. del is the amount (and direction) of the growth,
 * in residues, and op contains the options passed to the program. */
int pickMatch(fasstSolutionSet& matchList, Structure& S, FASST& F, int endIdx, int del, MstOptions& op) {
  // find suitable matches (those with enough context)
  vector<vector<Atom*> > nextRes; vector<int> matchInds;
  for (int m = 0; m < matchList.size(); m++) {
    vector<int> resIndices = F.getMatchResidueIndices(matchList[m], FASST::matchType::REGION);
    Structure target = F.getMatchStructure(matchList[m], false, FASST::matchType::FULL, true);
    RotamerLibrary::standardizeBackboneNames(target);
    // does this match have a residue past the central residue?
    int nextIdx = resIndices[endIdx] + del;
    if ((nextIdx >= target.residueSize()) || (nextIdx < 0) ||
        (target.getResidue(nextIdx).getChain() != target.getResidue(resIndices[endIdx]).getChain())) continue;
    nextRes.push_back((vector<Atom*>) AtomPointerVector(target.getResidue(nextIdx).getAtoms()).clone());
    matchInds.push_back(m);
  }
  cout << "\tout of these, found " << matchInds.size() << " suitable matches" << endl;
  if (matchInds.size() == 0) return -1;
  CartesianPoint weights(matchInds.size(), 1.0);
  // vector<vector<int>> clusters = clusterGuru.greedyCluster(nextRes, 0.3);
  // cout << "\tresulted in " << clusters.size() << " clusters" << endl;

  // unless we were asked to flatten matches, weight by RMSD
  if (!op.isGiven("f")) {
    mstreal rmsdFactor = F.getRMSDCutoff() / 10.0;
    for (int i = 0; i < matchInds.size(); i++) {
      weights[i] *= exp(-matchList[matchInds[i]].getRMSD() / rmsdFactor);
    }
  }

  // if a radius of gyration factor was given, weight by compactness
  if (op.isGiven("R")) {
    mstreal Rcoeff = op.getReal("R");
    AtomPointerVector atoms = S.getAtoms();
    for (int i = 0; i < matchInds.size(); i++) {
      AtomPointerVector tempAtoms = atoms;
      tempAtoms.insert(tempAtoms.end(), nextRes[i].begin(), nextRes[i].end());
      weights[i] *= exp(-tempAtoms.radiusOfGyration()/Rcoeff);
    }
  }

  // impose natural persistance length behavior of the chain
  if (op.isGiven("P")) {
    // TODO: just do the part that is growing! (or maybe at least clean out the backbone before calling this)
    Residue res;
    vector<Residue*> residues = S.getResidues();
    residues.push_back(&res);
    mstreal L = 10; // ideal native persistance length
    mstreal D = 33; // ideal native period of oscillations
    CartesianPoint errors(matchInds.size(), 0.0); // deviation from expectation for this solution
    vector<mstreal> cosines;
    for (int i = 0; i < matchInds.size(); i++) {
      res.copyAtoms(nextRes[i]);
      cosines = chainAutocorr(residues);
      for (int n = 0; n < cosines.size(); n++) {
        mstreal meanCosine = exp(-(n+1)/L)*cos((n+1)*2*M_PI/D); // valie expected for native proteins
        errors[i] += (cosines[n] - meanCosine)*(cosines[n] - meanCosine);
      }
    }
    if (!cosines.empty()) { // assuming there is enough of a chain to compute some inter-positional cosines
      for (int i = 0; i < matchInds.size(); i++) weights[i] *= exp(-errors[i]/MstUtils::max(errors));
      mstreal tot = weights.sum();
      for (int i = 0; i < matchInds.size(); i++) weights[i] /= tot;
cout << "weights: " << MstUtils::vecToString(weights) << endl;
    }
  }

  // pick either the best scoring option or probabilistically, depending on
  // whether randomization was requested
  int pick = 0;
  if (op.isGiven("r")) {
    weights = weights / weights.sum();
    for (int i = 1; i < weights.size(); i++) weights[i] += weights[i-1];
    mstreal r = MstUtils::randUnit();
    for (int i = 0; i < weights.size(); i++) {
      if ((weights[i] >= r) || (i == weights.size() - 1)) { pick = i; break; }
    }
  } else {
    MstUtils::max((vector<mstreal>) weights, 0, weights.size()-1, &pick);
  }
  for (int m = 0; m < nextRes.size(); m++) { AtomPointerVector(nextRes[m]).deletePointers(); } // clean up cloned fragments
  cout << "\tpicked match with index " << matchInds[pick] << endl;

  return matchInds[pick];
}

int main(int argc, char** argv) {
  MstOptions op;
  op.setTitle("Grows a chain from either the N- or the C-terminus, by the specified number of residues, aiming to make the overall structure as designable as possible. Options:");
  op.addOption("p", "starting PDB file.", true);
  op.addOption("L", "grow by this number of residues", true);
  op.addOption("n", "grow from the N-terminus.");
  op.addOption("c", "grow from the C-terminus.");
  op.addOption("cid", "chain ID of the chain to grow (if more than one chain in structure).");
  op.addOption("sid", "segment ID of the chain to grow (if more than one chain in structure). At most one of --cid or --sid should be specified.");
  op.addOption("rLib", "path to MST rotamer library.", true);
  op.addOption("o", "output base name.", true);
  op.addOption("b", "binary FASST database to use.", true);
  op.addOption("r", "grow randomly rather than deterministically. At each step, rather than picking the centroid of the top cluster, pick centroids randomly by cluster size.");
  op.addOption("m", "enforce a minimum of this many number of matches at each branch. By default, uses only an automatic RMSD cutoff, so some branches will terminate.");
  op.addOption("R", "radius of gyration factor. If specified, will preferrentially choose more compact solutions, with this exponential pre-factor (i.e., exp[-Rg/factor]).");
  op.addOption("P", "impose native-like persistence length of the protein chain. Probably best to use either --R or --P.");
  op.setOptions(argc, argv);
  if (op.isGiven("R") && (!op.isReal("R") || (op.getReal("R") <= 0))) MstUtils::error("-R must be a positive number!");
  int pm = 2;
  mstreal cdCut = 0.01;
  int maxRetries = 10;
  int retry = 0;
  int dN = 1; // number of residues to grow each time
  RMSDCalculator rc;
  Structure So(op.getString("p"));
  Structure S, Sp; RotamerLibrary::extractProtein(S, So);
  RotamerLibrary::standardizeBackboneNames(S);

  FASST F;
  F.readDatabase(op.getString("b"), 2);
  F.setRedundancyCut(0.5);
  RotamerLibrary RL(op.getString("rLib"));
  Clusterer clusterGuru(false);
  fusionParams fuserOpts; fuserOpts.setNumIters(1000); fuserOpts.setVerbose(false);
  srand(time(NULL) + (int) getpid());

  // TERMify loop
  fstream out; MstUtils::openFile(out, op.getString("o") + ".grow.pdb", ios_base::out);
  out << "MODEL " << 0 << endl; S.writePDB(out); out << "ENDMDL" << endl;
  for (int nc = 0; nc < 2; nc++) {
    int del = nc ? -1 : 1;
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
      Chain* cChain = &(nc ? S[0] : S.getChain(S.chainSize() - 1)); // default chain of residues to grow
      if (op.isGiven("cid")) {
        Chain* cChain = S.getChainByID(op.getString("cid"));
        MstUtils::assertCond(cChain != NULL, "did not find chain with ID " + op.getString("cid"));
      } else if (op.isGiven("sid")) {
        Chain* cChain = S.getChainBySegID(op.getString("sid"));
        MstUtils::assertCond(cChain != NULL, "did not find segment with ID " + op.getString("sid"));
      }
      Residue* cres = &(nc ? cChain->getResidue(0) : cChain->getResidue(cChain->residueSize() - 1));
      Structure term = TERMUtils::selectTERM({cres}, C, pm, cdCut, &termResIndices);
      int cresIdx = MstUtils::indexMap(termResIndices).at(cres->getResidueIndex());
      cout << AtomPointerVector(term.getAtoms()) << endl;

      // find matches to it
      cout << "searching..." << endl;
      F.setQuery(term, false);
      F.setRMSDCutoff(RMSDCalculator::rmsdCutoff(term));
      if (op.isGiven("m")) F.setMinNumMatches(op.getInt("m"));
      F.setMaxNumMatches(1000);
      fasstSolutionSet matchList = F.search();
      cout << "\tfound " << matchList.size() << " matches" << endl;

      int m = pickMatch(matchList, S, F, cresIdx, dN*del, op); // pick a suitable match for growing next residue
      if (m < 0) {
        if ((c == 0) || (retry > maxRetries)) {
          cout << "\tno suitable matches in first round, stopping current direction..." << endl;
          break;
        } else {
          cout << "\tbacking up one step to try a different random choice..." << endl;
          retry += 2; // we are about to take a step back, so will need two steps forward to advance
          cout << "\tbacking up one step to try a different random choice, try " << retry-1 << "..." << endl;
          i--; c--; S = Sp;
          continue;
        }
      }
      cout << "fusing..." << endl;
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

      // extend the match by dN residues
      for (int ii = 1; ii <= dN; ii++) {
        int nextIdx = resIndices[cresIdx] + ii*del;
        vector<Residue*> newPos(1, &(match.getResidue(nextIdx)));
        if (nc) {
          resTopo.insert(resTopo.begin(), newPos);
          isFixed.insert(isFixed.begin(), false);
        } else {
          resTopo.push_back(newPos);
          isFixed.push_back(false);
        }
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
// vector<mstreal> v = chainAutocorr(S.getResidues());
// cout << S.residueSize() << "\t" << v.size() << endl << MstUtils::vecToString(v, "\n") << endl; exit(-1);
    }
  }
  out.close();
}
