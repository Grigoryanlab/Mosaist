#include "msttypes.h"
#include "mstrotlib.h"
#include "mstcondeg.h"
#include "mstsequence.h"
#include "mstoptions.h"
#include "mstfasst.h"
#include <unistd.h>

using namespace MST;

void writeMatrix(const string& fileName, const vector<vector<int> >& subMatrix, bool echo = true) {
  fstream ofs;
  MstUtils::openFile(ofs, fileName, ios::out);
  for (int i = 0; i < subMatrix.size(); i++) {
    ofs << SeqTools::idxToSingle(i) << "\t";
    if (echo) cout << SeqTools::idxToSingle(i) << "\t";
  }
  ofs << endl;
  if (echo) cout << endl;
  for (int i = 0; i < subMatrix.size(); i++) {
    for (int j = 0; j < subMatrix.size(); j++) {
      ofs << subMatrix[i][j] << "\t";
      if (echo) cout << subMatrix[i][j] << "\t";
    }
    ofs << endl;
    if (echo) cout << endl;
  }
  ofs.close();
}

void selectAround(const vector<Residue*>& cenRes, int pm, Structure& frag, bool ignoreGaps = false, vector<int>* fragResIdx = NULL, vector<int>* fragCenResIdx = NULL) {
  Structure* S = cenRes[0]->getChain()->getParent();
  vector<int> selected(S->residueSize(), 0);
  for (int i = 0; i < cenRes.size(); i++) {
    Residue& res = *(cenRes[i]);
    Chain* C = res.getChain();
    int ri = res.getResidueIndex();
    int Li = C->getResidue(C->residueSize() - 1).getResidueIndex(); // last residue index in the chain
    for (int k = ri - pm; k <= ri + pm; k++) {
      if ((k < 0) || (k > Li)) continue;
      if (!selected[k]) selected[k] = 1;
      if (k == ri) selected[k] = 2;
    }
  }
  Chain* newChain = frag.appendChain("A", true);
  vector<Residue*> fragCenRes;
  for (int k = 0; k < selected.size(); k++) {
    if (selected[k] > 0) {
      // where there is a break in the selection, start a new chain
      if ((newChain->residueSize() > 0) && ((!ignoreGaps && (selected[k-1] == 0)) || (S->getResidue(k-1).getChain() != S->getResidue(k).getChain()))) {
        newChain = frag.appendChain("A", true);
      }
      Residue* newRes = new Residue(S->getResidue(k));
      newChain->appendResidue(newRes);
      if (fragResIdx != NULL) fragResIdx->push_back(k);
      if (selected[k] == 2) fragCenRes.push_back(newRes);
      if (ignoreGaps) newRes->setNum(k+1);
    }
  }
  if (fragCenResIdx != NULL) {
    for (int i = 0; i < fragCenRes.size(); i++) fragCenResIdx->push_back(fragCenRes[i]->getResidueIndex());
  }
}

Structure getTERM(Residue& res, RotamerLibrary& RL, vector<int>& keyResidues, vector<vector<int> >& resSpacing, mstreal cdcut = 0.0, bool verbose = false) {
  Structure& S = *(res.getStructure());
  ConFind C(&RL, S);
  contactList list = C.getContacts(&res, cdcut);
  // cout << "found " << list.size() << " contacts for residue " << res << " with CD cutoff " << cdcut << endl;
  vector<Residue*> residues(1, &res);
  for (int i = 0; i < list.size(); i++) {
    residues.push_back(list.residueB(i));
  }
  Structure frag, fragNoGaps; // the latter will not break chains due to selection gaps
  selectAround(residues, 1, frag, false, NULL, &keyResidues);
  if (verbose) frag.writePDB(cout);
  selectAround(residues, 1, fragNoGaps, true);

  // for each chain in the TERM, identify
  resSpacing.clear(); resSpacing.resize(fragNoGaps.chainSize());
  for (int i = 0; i < fragNoGaps.chainSize(); i++) {
    Chain& C = fragNoGaps[i];
    for (int ri = 0; ri < C.residueSize(); ri++) resSpacing[i].push_back(C[ri].getNum());
    if (verbose) cout << MstUtils::vecToString(resSpacing[i]) << endl;
  }
  return frag;
}

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Compute the TERM-induced amino-acid substitution matrix. Options:");
  op.addOption("d", "a FASST database.", true);
  op.addOption("N", "number of structures to randomly sample from the database.", true);
  op.addOption("n", "number of residues to randomly visit for each structure.", true);
  op.addOption("r", "path to rotamer library file.", true);
  op.addOption("o", "output base for the substitution table and the log file.", true);
  op.addOption("c", "contact degree cutoff for defining contacts.");
  op.addOption("v", "verbose output.");
  op.setOptions(argc, argv);
  RotamerLibrary RL(op.getString("r"));
  RMSDCalculator rc;
  FASST F;
  int N = 1000; // number of structures to randomly sample
  int n = 1;    // number
  int centResIdx;
  mstreal cdcut = op.getReal("c", 0.05);
  string matFile = op.getString("o") + ".mat";
  string logFile = op.getString("o") + ".log";
  bool verbose = op.isGiven("v");
  // srand(123456789);
  srand(time(NULL) + (int) getpid());
  fstream logfs; MstUtils::openFile(logfs, logFile, ios::out);

  // prepare substitution counts matrix
  vector<vector<int> > subMatrix(20, vector<int>(20, 0));

  // prepare FASST object
  F.setMemorySaveMode(true);
  F.readDatabase(op.getString("d"));
  F.options().setRedundancyCut(0.5);
  // F.setMaxNumMatches(1000); // these many matches are enough to compute reasonable frequencies
  // F.setMinNumMatches(20);   // need at least this many to compute anything reasonable

  for (int c = 0; c < op.getInt("N"); c++) {
    // pick a random structure from the database
    int ti = MstUtils::randInt(0, F.numTargets() - 1);
    Structure S(F.getTargetCopy(ti));

    for (int r = 0; r < op.getInt("n"); c++) {
      // pick a random residue
      int ri = MstUtils::randInt(0, S.residueSize() - 1);
      Residue& res = S.getResidue(ri);
      int aai = SeqTools::aaToIdx(SeqTools::toSingle(res.getName()));
      if (aai >= subMatrix.size()) continue;
      /* Will contain indices of the residues defining the TERM (first the
       * central residue, and then residues it is in contact with, if any).
       * Indices are into the TERM structure, not the original structure. */
      vector<int> keyResidues;
      vector<vector<int> > resSpacing; // will store the spacing between residues, in the original structure
      Structure term = getTERM(res, RL, keyResidues, resSpacing, cdcut, verbose);
      F.setQuery(term);
      F.setRMSDCutoff(rc.rmsdCutoff(resSpacing));
      if (verbose) cout << "term.ressidueSize() = " << term.residueSize() << ", term.chainSize() = " << term.chainSize() << ", RMSD cutoff = " << rc.rmsdCutoff(resSpacing) << endl;
      fasstSolutionSet matches = F.search();
      if (verbose) cout << "found " << matches.size() << " matches" << endl;
      logfs << S.getName() << " " << ti << " " << ri << " " << res << " " << keyResidues.size() << " "
            << term.residueSize() << " " << term.chainSize() << " " << matches.size();
      int Na = 0; // number of accepted matches
      // split matches by target
      map<int, vector<fasstSolution> > matchesByTarget;
      for (int i = 0; i < matches.size(); i++) {
        int idx = matches[i].getTargetIndex();
        if (idx == ti) continue; // don't count matches from the same structure as the query
        matchesByTarget[idx].push_back(matches[i]);
      }

      for (auto it = matchesByTarget.begin(); it != matchesByTarget.end(); ++it) {
        int idx = it->first;
        vector<fasstSolution>& sols = it->second;
        Structure target = F.getTargetCopy(idx);
        ConFind C(&RL, target);
        for (int i = 0; i < sols.size(); i++) {
          vector<int> residues = F.getMatchResidueIndices(sols[i], FASST::matchType::REGION);

          // make sure that the central residue within this match does not have
          // any additional contacts within its structure
          Residue& subRes = target.getResidue(residues[keyResidues[0]]);
          vector<Residue*> conts = C.getContactingResidues(&subRes, cdcut);
          set<int> keyResidueSet;
          for (int ri = 0; ri < keyResidues.size(); ri++) keyResidueSet.insert(residues[keyResidues[ri]]);
          bool reject = false;
          for (int ci = 0; ci < conts.size(); ci++) {
            int ind = conts[ci]->getResidueIndex();
            if (keyResidueSet.find(ind) == keyResidueSet.end()) {
              reject = true; // reject this match, because it has an extra contact with the central residue
            }
          }
          if (reject) {
            if (verbose) cout << "\tmatch " << i+1 << " from target " << idx << " rejected due to extra contacts..." << endl;
            continue;
          }
          if (verbose) cout << "\tmatch " << i+1 << " accepted!" << endl;

          // match is accepted, so count the substitution
          int aaj = SeqTools::aaToIdx(SeqTools::toSingle(subRes.getName()));
          if (aaj >= subMatrix.size()) continue;
          if (verbose) cout << "\t\tcounting substitution between " << res.getName() << " and " << subRes.getName() << endl;
          subMatrix[aai][aaj]++;
          Na++;
          writeMatrix(matFile, subMatrix, verbose);
        }
      }
      logfs << " " << Na << endl;
    }
  }
  logfs.close();
}
