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

using namespace std;
using namespace MST;

AtomPointerVector getBackbone(const Structure& S, vector<int> residues) {
  AtomPointerVector atoms;
  vector<string> bba = {"N", "CA", "C", "O"};
  for (int i = 0; i < residues.size(); i++) {
    Residue& res = S.getResidue(residues[i]);
    for (int j = 0; j < bba.size(); j++) {
      atoms.push_back(res.findAtom(bba[j]));
    }
  }
  return atoms;
}

mstreal getRadius(const Structure& S) {
  selector sel(S);
  AtomPointerVector atoms = sel.select("name N or name CA or name C or name O");
  mstreal rad = 0;
  for (int i = 0; i < atoms.size(); i++) {
    for (int j = i+1; j < atoms.size(); j++) {
      mstreal d = atoms[i]->distance(atoms[j]);
      if (d > rad) rad = d;
    }
  }
  return rad;
}

void selectAround(const vector<Residue*>& cenRes, int pm, Structure& frag, vector<int>& fragResIdx) {
  Structure* S = cenRes[0]->getChain()->getParent();
  vector<bool> selected(S->residueSize(), false);
  for (int i = 0; i < cenRes.size(); i++) {
    Residue& res = *(cenRes[i]);
    Chain* C = res.getChain();
    int ri = res.getResidueIndex();
    int Li = C->getResidue(C->residueSize() - 1).getResidueIndex(); // last residue index in the chain
    for (int k = ri - pm; k <= ri + pm; k++) {
      if ((k < 0) || (k > Li)) continue;
      selected[k] = true;
    }
  }
  Chain* newChain = frag.appendChain("A", true);
  for (int k = 0; k < selected.size(); k++) {
    if (selected[k]) {
      // where there is a break in the selection, start a new chain
      if ((newChain->residueSize() > 0) && ((!selected[k-1]) || (S->getResidue(k-1).getChain() != S->getResidue(k).getChain()))) {
        newChain = frag.appendChain("A", true);
      }
      newChain->appendResidue(new Residue(S->getResidue(k)));
      fragResIdx.push_back(k);
    }
  }
}

void addMatches(FASST& F, Structure& frag, const vector<int>& fragResIdx, vector<Structure*>& allMatches, vector<vector<Residue*> >& resTopo) {
  F.setQuery(frag, false);
  F.setRMSDCutoff(RMSDCalculator::rmsdCutoff(frag));
  F.setMaxNumMatches(10);
  F.setMinNumMatches(2);
  F.search();
  fasstSolutionSet matches = F.getMatches();
  cout << "found " << matches.size() << " matches" << endl;
  for (auto it = matches.begin(); it != matches.end(); ++it) {
    allMatches.push_back(new Structure(F.getMatchStructure(*it, false, FASST::matchType::REGION)));
    Structure& match = *(allMatches.back());
    for (int k = 0; k < fragResIdx.size(); k++) {
      resTopo[fragResIdx[k]].push_back(&(match.getResidue(k)));
    }
  }
}

int main(int argc, char** argv) {
  // TODO: debug Neilder-Meid optimization by comparing a simple case with Matlab
  // TODO: enable a setting in Fuser, whereby fully overlapping segments are scored
  //       (in terms of RMSD) via a weighted average, such that the lowest-RMSD
  //       segment makes a dominant contribution to the score
  MstOptions op;
  op.setTitle("Starting from some arbitrary conformation of a chain, iteratively build a TERM-based compact structure. Options:");
  op.addOption("p", "starting conformation PDB file.", true);
  op.addOption("rLib", "a path to an MST-formatter rotamer library.", true);
  op.addOption("d", "a database file with a list of PDB files.");
  op.addOption("b", "a binary database file. If both --d and --b are given, will overwrite this file with a corresponding binary database.");
  op.addOption("o", "output base name.", true);
  op.addOption("rad", "compactness radius. Default will be based on protein length.");
  op.addOption("cyc", "number of iteration cycles (10 by default).");
  op.addOption("f", "a quoted, space-separated list of 0-initiated residue integers to fix.");
  op.setOptions(argc, argv);
  RMSDCalculator rc;
  Structure I(op.getString("p"));
  FASST F;
  vector<int> fixed;
  F.setMemorySaveMode(true);
  if (op.isGiven("d")) {
    F.addTargets(MstUtils::fileToArray(op.getString("d")));
    if (op.isGiven("b")) {
      F.writeDatabase(op.getString("b"));
    }
  } else if (op.isGiven("b")) {
    F.readDatabase(op.getString("b"));
  } else {
    MstUtils::error("either --b or --d must be given!");
  }
  if (op.isGiven("f")) {
    fixed = MstUtils::splitToInt(op.getString("f"));
  }
  F.pruneRedundancy(0.5);
  RotamerLibrary RL(op.getString("rLib"));
  int pmSelf = 2, pmPair = 1;
  int Ni = 1000;
  fusionScores scores;

  // TERMify loop
  Structure S = I;
  mstreal R0 = getRadius(I);
  mstreal Rf = op.getReal("rad", pow(I.residueSize() * 1.0, 1.0/3)*5.0);
  int Ncyc = op.getInt("cyc", 10);
  fstream out; MstUtils::openFile(out, op.getString("o") + ".traj.pdb", ios_base::out);
  out << "MODEL " << 0 << endl; S.writePDB(out); out << "ENDMDL" << endl;
  for (int c = 0; c < Ncyc; c++) {
    out << "MODEL " << c+1 << endl;
    cout << "Cycle " << c+1 << "..." << endl;
    /* --- Decorate the current conformation with TERMs --- */
    // first self TERMs
    cout << "Searching for self TERMs..." << endl;
    vector<vector<Residue*> > resTopo(I.residueSize());
    // by inserting the starting structure into the topology first, any fixed
    // regions will be fixed to these starting coordinates
    for (int ri = 0; ri < S.residueSize(); ri++) resTopo[ri].push_back(&(S.getResidue(ri)));
    vector<Structure*> allMatches;
    for (int ci = 0; ci < S.chainSize(); ci++) {
      Chain& C = S[ci];
      for (int ri = 0; ri < C.residueSize(); ri++) {
        Structure frag; vector<int> fragResIdx;
        selectAround({&C[ri]}, pmSelf, frag, fragResIdx);
// frag.writePDB("/tmp/self." + MstUtils::toString(ri) + ".pdb");
        cout << "TERM around " << C[ri] << " ";
        addMatches(F, frag, fragResIdx, allMatches, resTopo);
      }
    }

    // then pair TERMs
    cout << "Searching for pair TERMs..." << endl;
    ConFind cfd(&RL, S);
    contactList L = cfd.getContacts(S, 0.01);
    vector<pair<Residue*, Residue*> > list = L.getOrderedContacts();
    for (int k = 0; k < list.size(); k++) {
      Residue* resA = list[k].first;
      Residue* resB = list[k].second;
      Structure frag; vector<int> fragResIdx;
      selectAround({resA, resB}, pmPair, frag, fragResIdx);
// frag.writePDB("/tmp/pair." + MstUtils::toString(k) + ".pdb");
      cout << "TERM around " << *resA << " x " << *resB << " ";
      addMatches(F, frag, fragResIdx, allMatches, resTopo);
    }

// fstream ofs;
// MstUtils::openFile(ofs, "/tmp/allmatches.pdb", ios_base::out);
// S.writePDB(ofs);
// for (int i = 0; i < allMatches.size(); i++) allMatches[i]->writePDB(ofs);
// ofs.close();

    /* --- Fuse --- */
    fusionParams opts; opts.setNumIters(Ni); opts.setVerbose(false);
    opts.setGradientDescentFlag(true);
    opts.setRepFC(1);
    opts.setCompFC(0.1);
    mstreal compactnessRadius = Rf;
    // (R0 - Rf)*exp(-c/10.0) + Rf; // exponential scaling
    // (Rf*(c + 1) + R0*(Ncyc - c - 1))/Ncyc; // linear scaling

    opts.setCompRad(compactnessRadius);
    cout << "will be trying to combine structure to a radius of " << compactnessRadius << endl;
    Structure fused = Fuser::fuse(resTopo, scores, fixed, opts);
    cout << "\t" << scores << endl;
    // opts.setStartingStructure(fused);
    // opts.setGradientDescentFlag(false);
    // fused = Fuser::fuse(resTopo, scores, vector<int>(), opts);
    // cout << "\t" << scores << endl;

    // align based on the fixed part, if anything was fixed (for ease of visualization)
    if (fixed.size() > 0) {
      AtomPointerVector before = getBackbone(S, fixed);
      AtomPointerVector after = getBackbone(fused, fixed);
      rc.align(after, before, fused);
    }

    /* --- write intermediate result and clean up--- */
    fused.writePDB(out); out << "ENDMDL" << endl;
    S = fused;
    for (int mi = 0; mi < allMatches.size(); mi++) delete(allMatches[mi]);
  }
  out.close();
  S.writePDB(op.getString("o") + ".fin.pdb");
}
