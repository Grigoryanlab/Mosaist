#include "msttypes.h"
#include "mstrotlib.h"
#include "mstcondeg.h"
#include "mstsequence.h"
#include "mstoptions.h"
#include "mstfasst.h"

using namespace MST;

void selectAround(const vector<Residue*>& cenRes, int pm, Structure& frag, vector<int>* fragResIdx = NULL, vector<int>* fragCenResIdx = NULL) {
  Structure* S = cenRes[0]->getChain()->getParent();
  vector<int> selected(S->residueSize(), 0);
  for (int i = 0; i < cenRes.size(); i++) {
    Residue& res = *(cenRes[i]);
    Chain* C = res.getChain();
    int ri = res.getResidueIndex();
    int Li = C->getResidue(C->residueSize() - 1).getResidueIndex(); // last residue index in the chain
    for (int k = ri - pm; k <= ri + pm; k++) {
      if ((k < 0) || (k > Li)) continue;
      selected[k] = 1;
      if (k == ri) selected[k] = 2;
    }
  }
  Chain* newChain = frag.appendChain("A", true);
  vector<Residue*> fragCenRes;
  for (int k = 0; k < selected.size(); k++) {
    if (selected[k] > 0) {
      // where there is a break in the selection, start a new chain
      if ((newChain->residueSize() > 0) && ((selected[k-1] == 0) || (S->getResidue(k-1).getChain() != S->getResidue(k).getChain()))) {
        newChain = frag.appendChain("A", true);
      }
      Residue* newRes = new Residue(S->getResidue(k));
      newChain->appendResidue(newRes);
      if (fragResIdx != NULL) fragResIdx->push_back(k);
      if (selected[k] == 2) fragCenRes.push_back(newRes);
    }
  }
  if (fragCenResIdx != NULL) {
    for (int i = 0; i < fragCenRes.size(); i++) fragCenResIdx->push_back(fragCenRes[i]->getResidueIndex());
  }
}

Structure getTERM(Residue& res, RotamerLibrary& RL, vector<int>& keyResidues, mstreal cdcut = 0.0) {
  Structure& S = *(res.getStructure());
  ConFind C(&RL, S);
  contactList list = C.getContacts(&res, cdcut);
cout << "found " << list.size() << " contacts for residue " << res << " with CD cutoff " << cdcut << endl;
  vector<Residue*> residues(1, &res);
  for (int i = 0; i < list.size(); i++) {
    residues.push_back(list.residueB(i));
  }
  Structure frag;
  selectAround(residues, 2, frag, NULL, &keyResidues);
cout << "frag.ressidueSize() = " << frag.residueSize() << endl; cout << "frag.chainSize() = " << frag.chainSize() << endl;
frag.writePDB("/tmp/frag.pdb");
cout << MstUtils::vecToString(keyResidues) << endl;
  return frag;
}

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Compute the TERM-induced amino-acid substitution matrix. Options:");
  op.addOption("d", "a FASST database.", true);
  op.addOption("N", "number of structures to randomly sample from the database.", true);
  op.addOption("n", "number of residues to randomly visit for each structure.", true);
  op.addOption("r", "path to rotamer library file.", true);
  op.addOption("c", "contact degree cutoff for defining contacts.", false);
  op.setOptions(argc, argv);
  RotamerLibrary RL(op.getString("r"));
  RMSDCalculator rc;
  FASST F;
  int N = 100; // number of structures to randomly sample
  int n = 10;  // number
  int centResIdx;
  mstreal cdcut = op.getReal("c", 0.05);
  srand(123456789);

  // prepare substitution counts matrix
  vector<vector<int> > subMatrix(20, vector<int>(20, 0));

  // prepare FASST object
  F.setMemorySaveMode(true);
  F.readDatabase(op.getString("d"));
  F.pruneRedundancy(0.5);

  for (int c = 0; c < op.getInt("N"); c++) {
    // pick a random structure from the database
    int ti = MstUtils::randInt(0, F.numTargets() - 1);
    Structure S(F.getTarget(ti));

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
      Structure term = getTERM(res, RL, keyResidues, cdcut);
      set<int> keyResidueSet = MstUtils::contents(keyResidues);
      F.setQuery(term);
      F.setRMSDCutoff(rc.rmsdCutoff(term));
      fasstSolutionSet matches = F.search();
      cout << "found " << matches.size() << " matches" << endl;
      for (int i = 0; i < matches.size(); i++) {
        if (matches[i].getTargetIndex() == ti) continue; // don't count matches from the same structure as the query
        Structure match = F.getMatchStructure(matches[i], false, FASST::matchType::FULL);
        vector<int> residues = F.getMatchResidueIndices(matches[i], FASST::matchType::REGION);

        // make sure that the central residue within this match does not have
        // any additional contacts within its structure
        ConFind C(&RL, match);
        Residue& subRes = match.getResidue(residues[keyResidues[0]]);
        contactList list = C.getContacts(&subRes, cdcut);
        bool reject = false;
        for (int ci = 0; ci < list.size(); ci++) {
          int ind = list.residueB(ci)->getResidueIndex();
          if (keyResidueSet.find(ind) == keyResidueSet.end()) {
            reject = true; // reject this match, because it has an extra contact with the central residue
          }
        }
        if (reject) {
          cout << "\tmatch " << i+1 << " rejected due to extra contacts..." << endl;
          continue;
        }
        cout << "\tmatch " << i+1 << " accepted!" << endl;

        // match is accepted, so count the substitution
        int aaj = SeqTools::aaToIdx(SeqTools::toSingle(subRes.getName()));
        if (aaj >= subMatrix.size()) continue;
        cout << "\t\tcounting substitution between " << res.getName() << " and " << subRes.getName() << endl;
        subMatrix[aai][aaj]++;
      }
    }
  }
}
