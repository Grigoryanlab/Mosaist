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
#include "mstoptions.h"
#include "mstmagic.h"
#include "mstrotlib.h"

using namespace std;
using namespace MST;

class bindOptions {
  public:
    bindOptions(FASST* _F = NULL, RotamerLibrary* _RL = NULL) {
      F = _F;
      RL = _RL;
      pm = 1;
      pcut = 0.05;
      setDefaultRMSDCutoffs();
    }
    FASST* getFASST() { return F; }
    RotamerLibrary* getRotamerLibrary() { return RL; }
    int contextLen() { return pm; }
    mstreal getRMSDCut() { return rmsdCut; }
    mstreal getRMSDCut2() { return rmsdCut2; }
    mstreal minProb() { return pcut; }
    string getContSectName() const { return contSec; }

    void setDefaultRMSDCutoffs() {
      rmsdCut = rc.rmsdCutoff((vector<int>) {2*pm + 1});
      rmsdCut2 = rc.rmsdCutoff({2*pm + 1, 2*pm + 1});
    }
    void setContextLen(int _pm) { pm = _pm; }
    void setFASST(FASST* _F) { F = _F; }
    void setRotamerLibrary(RotamerLibrary* _RL) { RL = _RL; }
    void setMinProb(mstreal _pcut) { pcut = _pcut; }
    void setContSectName(const string& _contSec) { contSec = _contSec; }

  private:
    FASST* F;
    RotamerLibrary* RL;
    int pm;
    mstreal pcut, rmsdCut, rmsdCut2;
    RMSDCalculator rc;
    string contSec;
};

class attachment {
  public:
    attachment() { score = 0; }
    attachment(const attachment& A) { S = A.S; bbSource = A.bbSource; bbAtt = A.bbAtt; score = A.score; }
    attachment(const Structure& _S, mstreal _score);
    bool isAnExampleOf(const Structure& ex);
    mstreal getScore() const { return score; }
    Structure& getStructure() { return S; }

  private:
    Structure S;
    AtomPointerVector bbSource, bbAtt; // backbone segments of the source and the attachment portions
    mstreal score;
};

attachment::attachment(const Structure& _S, mstreal _score) {
  score = _score;
  S = _S;
  MstUtils::assert(S.chainSize() >=2, "an attachment must have at least two chains", "attachment::attachment");
  RotamerLibrary::standardizeBackboneNames(S);
  MstUtils::assert(RotamerLibrary::hasFullBackbone(S), "not all backbone atoms are defined!", "attachment::attachment");
  AtomPointerVector bb = RotamerLibrary::getBackbone(S);
  int nbba = bb.size()/S.residueSize();
  bbSource = bb.subvector(0, S[0].residueSize() * nbba - 1);
  bbAtt = bb.subvector(S[0].residueSize() * nbba, (S.residueSize() - S[0].residueSize()) * nbba - 1);
}

bool attachment::isAnExampleOf(const Structure& ex) {
  AtomPointerVector bbEx = ex.getAtoms();
  MstUtils::assert(bbEx.size() == bbAtt.size(), "given structure not of the right size", "attachment::isAnExampleOf");
  RMSDCalculator rc;
  return rc.rmsd(bbAtt, bbEx) < 0.5;
}

vector<attachment> getAttachments(Residue* sR, bindOptions& opts);
mstreal scoreAttachment(Structure& pose, bindOptions& opts);
vector<attachment> getAttachmentsInContext(Residue* sR, bindOptions& opts);
mstreal scoreAttachment(attachment& A, bindOptions& opts);

int main(int argc, char** argv) {
  MstOptions op;
  op.setTitle("Given some specific surface site(s) on a structure, binds the best binding poses. Options:");
  op.addOption("p", "target PDB structure.", true);
  op.addOption("s", "surface site selection (procedure will be repeated for each residue in the selection).", true);
  op.addOption("db", "a binary FASST database file. This database needs to have the \"cont\" residue property section populated with contacts.", true);
  op.addOption("rLib", "rotamer library file path.", true);
  op.addOption("o", "output base name.", true);
  op.setOptions(argc, argv);
  string contSec = "conts";
  RMSDCalculator rc;
  bindOptions opts;
  RotamerLibrary RL;
  RL.readRotamerLibrary(op.getString("rLib"));

  // read all inputs
  Structure T(op.getString("p"));
  selector sel(T);
  vector<Residue*> surf = sel.selectRes(op.getString("s"));
  FASST F;
  F.setMemorySaveMode(true); // backbone only
  F.readDatabase(op.getString("db"));
  if (!F.isResiduePairPropertyPopulated(contSec)) MstUtils::error("the FASST database does not appear to have a contact section");
  F.pruneRedundancy(0.5);
  F.setMaxNumMatches(10000);
  opts.setFASST(&F);
  opts.setContSectName(contSec);

  for (int i = 0; i < surf.size(); i++) {
    Residue* sR = surf[i];
    cout << "visiting residue " << *sR << endl;

    vector<attachment> A = getAttachments(sR, opts);
    cout << "found " << A.size() << " attachments" << endl;
    vector<mstreal> scores(A.size(), 0);
    for (int i = 0; i < A.size(); i++) {
      scores[i] = scoreAttachment(A[i], opts);
      cout << "\tscore of attachment " << i << " --> " << scores[i] << endl;
    }
    if (scores.size() > 0) {
      int bi = 0;
      MstUtils::max(scores, 0, scores.size() - 1, &bi);
      cout << "best attachment has score " << scores[bi] << ", writing..." << endl;
      A[bi].getStructure().writePDB("/tmp/att.pdb");
    }

    // Structure anchor;
    // TERMUtils::selectTERM(vector<Residue*>(1, sR), anchor, opts.contextLen());
    // F.setRMSDCutoff(opts.getRMSDCut());
    // F.setQuery(anchor);
    // fasstSolutionSet sols = F.search();
    // vector<Structure> matches; F.getMatchStructures(sols, matches);
    // vector<vector<Atom*> > contactTERMs;
    // int Ne = 0, Nc = 0;
    // for (int k = 0; k < matches.size(); k++) {
    //   Residue cR = matches[k].getResidue(opts.contextLen());
    //   if (!cR.isNamed(sR->getName())) continue;
    //
    //   // iterate over all contacts
    //   int ti = sols[k].getTargetIndex();
    //   int ri = (sols[k].getAlignment())[0] + opts.contextLen();
    //   Structure* mT = F.getTarget(ti);
    //   if (F.hasResiduePairProperties(ti, contSec, ri)) {
    //     map<int, mstreal> C = F.getResiduePairProperties(ti, contSec, ri);
    //     for (auto rj = C.begin(); rj != C.end(); ++rj) {
    //       vector<Atom*> contactTERM;
    //       TERMUtils::exciseTERM({&(mT->getResidue(ri)), &(mT->getResidue(rj->first))}, contactTERM, opts.contextLen());
    //       contactTERMs.push_back(contactTERM);
    //       Nc++;
    //     }
    //   } else {
    //     Ne++; // counts as non-contacting "exposed" residue
    //   }
    // }
    //
    // // cluster contact TERMs
    // cout << "\tclustering " << contactTERMs.size() << " terms" << endl;
    // vector<vector<int> > cIs = clust.greedyCluster(contactTERMs, opts.getRMSDCut2(), 10000);
    // for (int ci = 0; ci < cIs.size(); ci++) {
    //   printf("cluster %02d: %d out of %d + %d = %f\n", ci, (int) cIs[ci].size(), Nc, Ne, (cIs[ci].size()*1.0)/(Nc + Ne));
    //   mstreal p = (cIs[ci].size()*1.0)/(Nc + Ne);
    //   if (p > opts.minProb()) {
    //     AtomPointerVector cent(contactTERMs[cIs[ci][0]]);
    //     Structure centS(cent);
    //     rc.align(cent.subvector(0, cent.size()/2), F.getQuerySearchedAtoms(), centS);
    //     centS.writePDB(op.getString("o") + "." + MstUtils::toString(i) + "-" + MstUtils::toString(ci) + ".pdb");
    //   } else {
    //     exit(-1);
    //   }
    // }
  }
}

mstreal scoreAttachment(attachment& A, bindOptions& opts) {
  return scoreAttachment(A.getStructure(), opts);
}

mstreal scoreAttachment(Structure& pose, bindOptions& opts) {
  if (pose.chainSize() != 2) MstUtils::error("expected two chains", "scoreBindingPose");
  mstreal score = 0;
  for (int ord = 0; ord < 2; ord++) {
    int ri = ord ? 1 : 0;
    int rj = ord ? 0 : 1;
    Residue& sR = pose[ri].getResidue(pose[ri].residueSize()/2);
    Residue& iR = pose[rj].getResidue(pose[rj].residueSize()/2);
    vector<attachment> A = getAttachments(&(sR), opts);
    bool found = false;
    for (int i = 0; i < A.size(); i++) {
      if (A[i].isAnExampleOf(pose[rj].getAtoms())) {
        score += A[i].getScore();
        found = true;
        break;
      }
    }
    if (!found) return 10E10;
  }
  return score;
}

vector<attachment> getAttachments(Residue* sR, bindOptions& opts) {
  cout << "getting attachments for residue " << *sR << "..." << endl;
  Clusterer clust(true);
  RMSDCalculator rc;

  vector<attachment> A;
  Structure anchor;
  TERMUtils::selectTERM(vector<Residue*>(1, sR), anchor, opts.contextLen());
  FASST& F = *(opts.getFASST());
  F.setRMSDCutoff(opts.getRMSDCut());
  F.setQuery(anchor);
  fasstSolutionSet sols = F.search();
  vector<Structure> matches; F.getMatchStructures(sols, matches);
  vector<vector<Atom*> > contactTERMs;
  int Ne = 0, Nc = 0;
  for (int k = 0; k < matches.size(); k++) {
    Residue cR = matches[k].getResidue(opts.contextLen());
    if (!cR.isNamed(sR->getName())) continue;

    // iterate over all contacts
    int ti = sols[k].getTargetIndex();
    int ri = (sols[k].getAlignment())[0] + opts.contextLen();
    Structure* mT = F.getTarget(ti);
    if (F.hasResiduePairProperties(ti, opts.getContSectName(), ri)) {
      map<int, mstreal> C = F.getResiduePairProperties(ti, opts.getContSectName(), ri);
      for (auto rj = C.begin(); rj != C.end(); ++rj) {
        vector<Atom*> contactTERM;
        TERMUtils::exciseTERM({&(mT->getResidue(ri)), &(mT->getResidue(rj->first))}, contactTERM, opts.contextLen());
        contactTERMs.push_back(contactTERM);
        Nc++;
      }
    } else {
      Ne++; // counts as non-contacting "exposed" residue
    }
  }

  // cluster contact TERMs
  cout << "\tclustering " << contactTERMs.size() << " terms" << endl;
  vector<vector<int> > cIs = clust.greedyCluster(contactTERMs, opts.getRMSDCut2(), 10000);
  for (int ci = 0; ci < cIs.size(); ci++) {
    printf("\t\tcluster %02d: %d out of %d + %d = %f\n", ci, (int) cIs[ci].size(), Nc, Ne, (cIs[ci].size()*1.0)/(Nc + Ne));
    mstreal p = (cIs[ci].size()*1.0)/(Nc + Ne);
    if (p > opts.minProb()) {
      AtomPointerVector cent(contactTERMs[cIs[ci][0]]);
      Structure centS(cent);
      rc.align(cent.subvector(0, cent.size()/2), F.getQuerySearchedAtoms(), centS);
      A.push_back(attachment(centS, log(p)));
    }
  }
  return A;
}

// vector<attachment> getAttachmentsInContext(Residue* sR, bindOptions& opts) {
//   vector<attachment> A = getAttachments(sR, opts);
//   for (int i = 0; i < A.size(); i++) {
//     Structure* parentStructure = sR->getStructure();
//     MstUtils::assert(parentStructure != NULL, "source residue does not belong to a structure!", "getAttachmentsInContext");
//     Structure S(parentStructure);
//     Chain* nch = S.addChain("A");
//     nch->appendResidue(new Residue(sR));
//     ConFind C(S, opts.getRotamerLibrary());
//     vector<Residue*> contResidues = C.getContactingResidues(sR, 0.001???);
//     // TODO: also do based on backbone-backbone distance!
//     // for each contacting residue
//     // isolate attachment with sR
//     // score that attachment
//     // if individual score falls below some threshold, give up on this attachment
//     // otherwise update its score to increment by the score of this contact
//   }
//   return A;
// }
