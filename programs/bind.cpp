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
      pcut = 0.0;
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
    attachment(const attachment& A) : attachment(A.S, A.score) { }
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
  MstUtils::assertCond(S.chainSize() >=2, "an attachment must have at least two chains", "attachment::attachment");
  RotamerLibrary::standardizeBackboneNames(S);
  MstUtils::assertCond(RotamerLibrary::hasFullBackbone(S), "not all backbone atoms are defined!", "attachment::attachment");
  AtomPointerVector bb = RotamerLibrary::getBackbone(S);
  int nbba = bb.size()/S.residueSize();
  bbSource = bb.subvector(0, S[0].residueSize() * nbba);
  bbAtt = bb.subvector(S[0].residueSize() * nbba, S.residueSize() * nbba);
}

bool attachment::isAnExampleOf(const Structure& ex) {
  AtomPointerVector bbEx = ex.getAtoms();
  MstUtils::assertCond(bbEx.size() == bbAtt.size(), "given structure not of the right size", "attachment::isAnExampleOf");
  RMSDCalculator rc;
cout << "\t\tinside attachment::isAnExampleOf: " << rc.rmsd(bbAtt, bbEx) << endl;
if ((rc.rmsd(bbAtt, bbEx) > 0.5) && (rc.rmsd(bbAtt, bbEx) < 2.0)) {
  Structure(bbAtt).writePDB("/tmp/bb-att.pdb");
  S.writePDB("/tmp/bb-all.pdb");
  Structure(bbEx).writePDB("/tmp/bb-ex.pdb");
  exit(-1);
}
  return rc.rmsd(bbAtt, bbEx) < 0.5;
}

vector<attachment> getAttachments(Residue* sR, bindOptions& opts);
vector<attachment> getAttachmentsInContext(Residue* sR, bindOptions& opts);
mstreal scoreAttachment(Structure& pose, bindOptions& opts);
mstreal scoreAttachment(attachment& A, bindOptions& opts);
mstreal scoreAttachment(Residue* sR, Structure& contactTERM, bindOptions& opts);
vector<vector<int> > clusterContactTERMs(const vector<Structure>& contactTERMs, bindOptions& opts);
int mineAttachments(Residue* sR, bindOptions& opts, vector<Structure>& contactTERMs);

int main(int argc, char** argv) {
  MstOptions op;
  op.setTitle("Given some specific surface site(s) on a structure, binds the best binding poses. Options:");
  op.addOption("p", "target PDB structure.", true);
  op.addOption("s", "surface site selection (procedure will be repeated for each residue in the selection).", true);
  op.addOption("db", "a binary FASST database file. This database needs to have the \"conts\" residue property section populated with contacts.", true);
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
  F.readDatabase(op.getString("db"), 2);
  if (!F.isResiduePairPropertyPopulated(contSec)) MstUtils::error("the FASST database does not appear to have a contact section");
  F.setRedundancyCut(0.5);
  F.setMaxNumMatches(1000);
  opts.setFASST(&F);
  opts.setContSectName(contSec);

  for (int i = 0; i < surf.size(); i++) {
    Residue* sR = surf[i];
    cout << "visiting residue " << *sR << endl;

    vector<attachment> A = getAttachments(sR, opts);
    cout << "found " << A.size() << " attachments" << endl;
    vector<mstreal> scores;
    for (int i = 0; i < MstUtils::min((int) A.size(), 8); i++) {
      scores.push_back(scoreAttachment(A[i], opts));
      cout << "\tscore of attachment " << i << " --> " << scores[i] << endl;
    }
    if (scores.size() > 0) {
      int bi = 0;
      MstUtils::min(scores, 0, scores.size() - 1, &bi);
      cout << "best attachment has score " << scores[bi] << ", writing..." << endl;
      A[bi].getStructure().writePDB("/tmp/att.pdb");
    }
  }
}

mstreal scoreAttachment(attachment& A, bindOptions& opts) {
  return scoreAttachment(A.getStructure(), opts);
}

mstreal scoreAttachment(Structure& pose, bindOptions& opts) {
  if (pose.chainSize() != 2) MstUtils::error("expected two chains", "scoreBindingPose");
  Residue& rA = pose[0].getResidue(pose[0].residueSize()/2);
  Residue& rB = pose[1].getResidue(pose[1].residueSize()/2);
  Structure& sAB = pose;
  Structure sBA;
  sBA.appendChain(new Chain(pose[1]));
  sBA.appendChain(new Chain(pose[0]));
  return scoreAttachment(&rA, sAB, opts) + scoreAttachment(&rB, sBA, opts);
}

int mineAttachments(Residue* sR, bindOptions& opts, vector<Structure>& contactTERMs) {
  cout << "mining attachments for residue " << *sR << "..." << endl;
  Structure anchor;
  TERMUtils::selectTERM(vector<Residue*>(1, sR), anchor, opts.contextLen());
  FASST& F = *(opts.getFASST());
  F.setRMSDCutoff(opts.getRMSDCut());
  F.setQuery(anchor);
  cout << "\tsearching for a local " << anchor.chainSize() << "-segment TERM..." << endl;
  fasstSolutionSet sols = F.search();
  vector<Structure> matches; F.getMatchStructures(sols, matches);
  int Ne = 0, Nc = 0;
  cout << "\tfound " << matches.size() << " matches, excising local context..." << endl;
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
        Structure contactTERM;
        if (!TERMUtils::exciseTERM({&(mT->getResidue(ri)), &(mT->getResidue(rj->first))}, contactTERM, opts.contextLen())) continue;
        if (contactTERM.residueSize() != 2*(2*opts.contextLen() + 1)) continue; // if ran into end of chain, some segments are not complete
        contactTERMs.push_back(contactTERM);
        Nc++;
      }
    } else {
      Ne++; // counts as non-contacting "exposed" residue
    }
  }

  return Ne;
}

vector<vector<int> > clusterContactTERMs(const vector<Structure>& contactTERMs, bindOptions& opts) {
  Clusterer clust(true);
  vector<vector<Atom*> > contactTERMatoms(contactTERMs.size());
  for (int i = 0; i < contactTERMs.size(); i++) contactTERMatoms[i] = contactTERMs[i].getAtoms();
  cout << "\tclustering " << contactTERMs.size() << " TERMs..." << endl;
  return clust.greedyCluster(contactTERMatoms, opts.getRMSDCut2(), 10000);
}

vector<attachment> getAttachments(Residue* sR, bindOptions& opts) {
  // given the soure residue, get the most typical attachment geometries
  vector<Structure> contactTERMs;
  int Ne = mineAttachments(sR, opts, contactTERMs);
  int Nc = contactTERMs.size();

  // cluster these
  vector<vector<int> > cIs = clusterContactTERMs(contactTERMs, opts);

  // create attachment objects with scores
  RMSDCalculator rc;
  vector<attachment> A;
  for (int ci = 0; ci < cIs.size(); ci++) {
    printf("\t\tcluster %02d: %d out of %d + %d = %f\n", ci, (int) cIs[ci].size(), Nc, Ne, (cIs[ci].size()*1.0)/(Nc + Ne));
    mstreal p = (cIs[ci].size()*1.0)/(Nc + Ne);
    if (p > opts.minProb()) {
      Structure cent(contactTERMs[cIs[ci][0]]);
      rc.align(cent[0].getAtoms(), opts.getFASST()->getQuerySearchedAtoms(), cent);
      A.push_back(attachment(cent, -log(p)));
    }
  }
  return A;
}

mstreal scoreAttachment(Residue* sR, Structure& contactTERM, bindOptions& opts) {
  // given the soure residue, get the most typical attachment geometries
  vector<Structure> contactTERMs;
  int Ne = mineAttachments(sR, opts, contactTERMs);
  int Nc = contactTERMs.size();

  // cluster these
  vector<vector<int> > cIs = clusterContactTERMs(contactTERMs, opts);

  // see where the geometry in question clusters
  RMSDCalculator rc;
  for (int ci = 0; ci < cIs.size(); ci++) {
    Structure cent = contactTERMs[cIs[ci][0]];
    mstreal p = (cIs[ci].size()*1.0)/(Nc + Ne);
    if (rc.bestRMSD(cent.getAtoms(), contactTERM.getAtoms()) < opts.getRMSDCut2()) {
      return -log(p);
    }
  }
  return 10e10;
}

// mstreal scoreContext(Chain& _nC, const Structure& _S, bindOptions& opts) {
//   Structure S = _S;
//   Chain* nC = new Chain(_nC);
//   S.appendChain(nC);
//   ConFind C(S, opts.getRotamerLibrary());
//
//   // iterate over all residues of the new chain and check all of its contacts
//   for (int i = opts.contextLen(); i < nC->residueSize() - opts.contextLen(); i++) {
//     Residue& sR = (*nC)[i];
//     vector<Residue*> contResidues = C.getContactingResidues(&sR, 0.01);
//   }
// }

// vector<attachment> getAttachmentsInContext(Residue* sR, bindOptions& opts) {
//   vector<attachment> A = getAttachments(sR, opts);
//   for (int i = 0; i < A.size(); i++) {
//     Structure* parentStructure = sR->getStructure();
//     MstUtils::assertCond(parentStructure != NULL, "source residue does not belong to a structure!", "getAttachmentsInContext");
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
