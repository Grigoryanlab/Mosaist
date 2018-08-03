#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>
#include <typeinfo>
#include <stdio.h>

#include "msttypes.h"
#include "mstfasst.h"
#include "mstfasstcache.h"
#include "mstcondeg.h"
#include "mstfuser.h"
#include "mstcondeg.h"
#include "mstoptions.h"
#include "mstmagic.h"
#include "mstsystem.h"

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

vector<Structure*> getMatches(fasstCache& C, Structure& frag, vector<int>& fragResIdx, int need = 5, const vector<int>& centIdx = vector<int>()) {
  // want at least "need" matches in the end; estimate how many to ask for at first
  int seqConsts = 0;
  for (int i = 0; i < centIdx.size(); i++) {
    if (!SeqTools::isUnknown(frag.getResidue(centIdx[i]).getName())) seqConsts++;
  }
  int Nmin = need*(pow(20, seqConsts)*2 + 3);
  FASST* F = C.getFASST();
  F->setQuery(frag, false);
  F->setRMSDCutoff(RMSDCalculator::rmsdCutoff(frag));
  F->options().setMaxNumMatches(Nmin*2);
  // F->options().setMinNumMatches(Nmin);

  vector<Structure*> matchStructures;
  // limit iterations, because it is technically possible that a match meeting
  // the sequence constraints does not exist, at which point we will just give up
  for (int c = 0; c < 10; c++) { // TODO: instead, iterate until consecutive searches do not increase number of hits
    fasstSolutionSet matches = C.search(true);
    for (auto it = matches.begin(); it != matches.end(); ++it) {
      if (centIdx.size() > 0) {
        // check for sequence compatibility
        Sequence mseq = F->getMatchSequence(*it);
        bool comp = true;
        for (int i = 0; i < centIdx.size(); i++) {
          Residue& res = frag.getResidue(centIdx[i]);
          if (SeqTools::isUnknown(res.getName())) continue;
          if (SeqTools::aaToIdx(res.getName()) != mseq[centIdx[i]]) { comp = false; break; }
        }
        if (!comp) continue; // skip non-compatible solutions
      }
      matchStructures.push_back(new Structure(F->getMatchStructure(*it, false, FASST::matchType::REGION)));
      Structure& match = *(matchStructures.back());
      MstUtils::assert(match.residueSize() == fragResIdx.size(), "unexpected match size");
      for (int k = 0; k < match.residueSize(); k++) {
        match.getResidue(k).setNum(fragResIdx[k]); // make residue numbers store indices into the original structure
      }
    }
    cout << "\tfound " << matchStructures.size() << " matches" << endl;
    if (matchStructures.size() >= need) break;
    for (int i = 0; i < matchStructures.size(); i++) delete(matchStructures[i]);
    matchStructures.clear();
    if (matches.size() == F->options().getMaxNumMatches()) {
      F->options().setMaxNumMatches(F->options().getMaxNumMatches()*2);
    } else {
      F->options().setRMSDCutoff(F->options().getRMSDCutoff()*1.1);
    }
    cout << "\t\tinsufficient, increasing params to " << F->options().getRMSDCutoff() << " / " << F->options().getMaxNumMatches() << endl;
  }
  return matchStructures;
}

void addMatches(vector<Structure*>& matchStructures, vector<vector<Residue*> >& resTopo, fstream& matchOut, int ri = -1) {
  for (int i = 0; i < matchStructures.size(); i++) {
    if ((ri >= 0) && (i != ri)) continue;
    Structure& match = *(matchStructures[i]);
    if (matchOut.is_open()) {
      match.writePDB(matchOut);
      matchOut << "END" << endl;
    }
    for (int k = 0; k < match.residueSize(); k++) {
      Residue& res = match.getResidue(k);
      resTopo[res.getNum()].push_back(&res); // residue numbers store indices into the originating structure
    }
  }
}

bool mc(mstreal oldScore, mstreal newScore, mstreal kT) {
  return ((newScore < oldScore) || (MstUtils::randUnit() < exp((oldScore - newScore)/kT)));
}

fusionTopology getTopo(int L, vector<vector<Structure*> >& allMatches, vector<int>& picks, fstream& matchOut, Structure* global = NULL) {
  fusionTopology resTopo(L);
  for (int si = 0; si < allMatches.size(); si++) {
    Structure& match = *(allMatches[si][picks[si]]);
    resTopo.addFragment(match);
    if (matchOut.is_open()) { match.writePDB(matchOut); matchOut << "END" << endl; }
  }
  if (global != NULL) {
    MstUtils::assert(global->residueSize() == L, "the global target structure specified has an unexpected number of residues for the topology");
    resTopo.addFragment(*global, MstUtils::range(0, L), -10.0);
    if (matchOut.is_open()) { global->writePDB(matchOut); matchOut << "END" << endl; }
  }
  return resTopo;
}

fusionTopology getTopo(int L, vector<vector<Structure*> >& allMatches, vector<vector<int> >& picks, fstream& matchOut, Structure* global = NULL) {
  fusionTopology resTopo(L);
  for (int si = 0; si < allMatches.size(); si++) {
    for (int j = 0; j < picks[si].size(); j++) {
      Structure& match = *(allMatches[si][picks[si][j]]);
      resTopo.addFragment(match);
      if (matchOut.is_open()) { match.writePDB(matchOut); matchOut << "END" << endl; }
    }
  }
  if (global != NULL) {
    MstUtils::assert(global->residueSize() == L, "the global target structure specified has an unexpected number of residues for the topology");
    resTopo.addFragment(*global, MstUtils::range(0, L), -10.0);
    if (matchOut.is_open()) { global->writePDB(matchOut); matchOut << "END" << endl; }
  }
  return resTopo;
}

AtomPointerVector getCorrespondingAtoms(Structure& from, Structure& like) {
  MstUtils::assert(from.residueSize() == like.residueSize(), "the two structures must have the same number of residues", "getCorrespondingAtoms()");
  AtomPointerVector atoms;
  for (int ri = 0; ri < like.residueSize(); ri++) {
    Residue& fromRes = from.getResidue(ri);
    Residue& likeRes = like.getResidue(ri);
    for (int ai = 0; ai < likeRes.atomSize(); ai++) {
      atoms.push_back(fromRes.findAtom(likeRes[ai].getName()));
    }
  }
  return atoms;
}

mstreal totalScore(fusionScores& scoreObj, Structure& fused, AtomPointerVector& init, bool report = false) {
  RMSDCalculator rc;
  // mstreal a = scoreObj.getScore()/10;
  // mstreal a = scoreObj.getTotRMSDScore();
  // mstreal a = scoreObj.getScore();
  // mstreal b = rc.bestRMSD(init, fused.getAtoms());
  // if (report) cout << "fuser: " << scoreObj << "; total = " << a << " - " << b;
  if (report) cout << "fuser: " << scoreObj << "; RMSD from start = " << rc.bestRMSD(init, fused.getAtoms());
  // return a - b;
  return scoreObj.getScore();
}

void copySequence(const Structure& from, Structure& to) {
  if (from.residueSize() != to.residueSize()) MstUtils::error("structures of different sizes", "copySequence");
  for (int i = 0; i < from.residueSize(); i++) {
    to.getResidue(i).setName(from.getResidue(i).getName());
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
  op.addOption("s", "account for sequence by working only with TERMs that contain relevant amino acids from the target. Amino acids named UNK will be ignored, whichh is convenient for disregarding parts of the sequence (e.g., those parts that will be designed later).");
  op.addOption("n", "pick this many matches for each TERM for fusion (default is 1). At each iteration, one match for a randomly-selected TERM will be randomly substituted.");
  op.addOption("r", "if specified, will randomly pick the matches for the first iteration. By default, will take the top match for each TERM or the top --n matches, if --n is specified.");
  op.addOption("m", "if specified, will try to move away from the original structure by including a term in the objective function that rewards large RMSDs to it. Default is no.");
  op.addOption("cyc", "number of cycles--i.e., number of times fresh TERMs are searched for (10 by default).");
  op.addOption("iter", "number of iterations per cycle (1 by default). At the start of each iteration, the overall structure is reinitialized to the current structure");
  op.addOption("a", "alternate between optimizing without and with internal coordinate constraints. Can be useful for getting out of local minima.");
  op.addOption("f", "a quoted, space-separated list of 0-initiated residue integers to fix.");
  op.addOption("fs", "a selection string for residues to fix.");
  op.addOption("us", "a selection string for residues to mark as having unknown identity (i.e., their identity will not matter if accounting for sequence).");
  op.addOption("rad", "compactness radius. Default will be based on protein length.");
  op.addOption("c", "path to a FASST cache file to use for initializing the cache.");
  op.addOption("w", "flag; if specified, the FASST cache will be periodically updated.");
  op.addOption("app", "flag; if specified, will append to the output PDB file (e.g., for the purpose of accumulating a trajectory from multiple runs).");
  if (op.isGiven("f") && op.isGiven("fs")) MstUtils::error("only one of --f or --fs can be given!");
  MstUtils::setSignalHandlers();
  op.setOptions(argc, argv);
  RMSDCalculator rc;
  Structure I(op.getString("p"));
  FASST F;
  vector<int> fixed;
  F.setMemorySaveMode(true);
  if (op.isGiven("c") && op.getString("c").empty()) MstUtils::error("--c must be a valid file path");
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
    cout << "fix specification gave " << fixed.size() << " residues, fixing..." << endl;
  }
  if (op.isGiven("fs")) {
    selector sel(I);
    vector<Residue*> fixedResidues = sel.selectRes(op.getString("fs"));
    cout << "fix selection gave " << fixedResidues.size() << " residues, fixing..." << endl;
    map<Residue*, int> indices = MstUtils::indexMap(I.getResidues());
    fixed.resize(fixedResidues.size());
    for (int i = 0; i < fixedResidues.size(); i++) fixed[i] = indices[fixedResidues[i]];
  }
  if (op.isGiven("us")) {
    selector sel(I);
    vector<Residue*> unkResidues = sel.selectRes(op.getString("us"));
    cout << "unknown sequence selection gave " << unkResidues.size() << " residues, marking unknown..." << endl;
    for (int i = 0; i < unkResidues.size(); i++) unkResidues[i]->setName("UNK");
  }
  int numPerTERM = op.getInt("n", 1);

  F.setRedundancyCut(0.5);
  RotamerLibrary RL(op.getString("rLib"));
  int pmSelf = 2, pmPair = 1;
  int Ni = 1000, lastWriteTime;
  fusionScores bestScore, currScore;
  fstream out, shellOut, dummy;
  contactList L;
  string tag = "TERMify-" + MstSys::getUserName();

  // If fixed residues were defined, then not all contacts will be needed; only
  // contacts that can ultimately implicate a non-fixed residue are of interest.
  // These are either contacts with non-fixed residues OR contacts with fixed
  // residues that are within +/- pmPair of non-fixed residues.
  vector<int> contResis;
  if (!fixed.empty()) {
    vector<bool> getConts(I.residueSize(), false);
    set<int> nonFixed = MstUtils::contents(MstUtils::setdiff(MstUtils::range(0, I.residueSize()), fixed));
    int k = 0;
    for (int ci = 0; ci < I.chainSize(); ci++) {
      Chain& C = I[ci];
      for (int ri = 0; ri < C.residueSize(); ri++) {
        if (nonFixed.find(k) != nonFixed.end()) {
          getConts[k] = true;
        } else {
          for (int rri = MstUtils::max(0, ri - pmSelf); rri < MstUtils::min(C.residueSize(), ri + pmSelf + 1); rri++) {
            if (nonFixed.find(k + rri - ri) != nonFixed.end()) {
              getConts[k] = true;
              break;
            }
          }
        }
        k++;
      }
    }
    for (int i = 0; i < getConts.size(); i++) {
      if (getConts[i]) contResis.push_back(i);
    }
  }

  // create a cache with size proprtional to the number of flexible residues
  fasstCache cache(&F, (I.residueSize() - fixed.size())*10);
  if (op.isGiven("c") && MstSys::fileExists(op.getString("c"))) {
    MstSys::getNetLock(tag, true);
    cout << "reading cache from " << op.getString("c") << "... " << endl;
    cache.read(op.getString("c"));
    MstSys::releaseNetLock(tag);
  }

  // TERMify loop
  Structure S = I.reassignChainsByConnectivity(); // fixed residues are already selected, but this does not change residue order
  RotamerLibrary::standardizeBackboneNames(S);
  mstreal R0 = getRadius(I);
  mstreal Rf = op.getReal("rad", pow(I.residueSize() * 1.0, 1.0/3)*5.0);
  int Ncyc = op.getInt("cyc", 10);
  MstUtils::openFile(out, op.getString("o") + ".traj.pdb", op.isGiven("app") ? ios::app : ios::out);
  out << "MODEL " << 0 << endl; S.writePDB(out); out << "ENDMDL" << endl;
  for (int c = 0; c < Ncyc; c++) {
    cout << "Cycle " << c+1 << "..." << endl;
    if (c == 0) {
      MstUtils::openFile(shellOut, op.getString("o") + ".shell.init.pdb", ios::out);
    } else if (c == Ncyc - 1) {
      MstUtils::openFile(shellOut, op.getString("o") + ".shell.pdb", ios::out);
    }
    if (op.isGiven("c") && op.isGiven("w")) {
      if ((c == 0) || (time(NULL) - lastWriteTime > 5*60)) { // write every five minutes or so
        cout << "dumping cache to " << op.getString("c") << "... " << endl;
        MstSys::getNetLock(tag);
        cache.write(op.getString("c"));
        MstSys::releaseNetLock(tag);
        lastWriteTime = time(NULL);
      }
    }

    /* --- Decorate the current conformation with TERMs --- */
    // first self TERMs
    cout << "Searching for self TERMs..." << endl;
    vector<vector<Structure*>> allMatches;
    for (int ci = 0; ci < S.chainSize(); ci++) {
      Chain& C = S[ci];
      for (int ri = 0; ri < C.residueSize(); ri++) {
        Structure frag; vector<int> fragResIdx;
        vector<int> centIdx = TERMUtils::selectTERM({&C[ri]}, frag, pmSelf, &fragResIdx);
        if (MstUtils::setdiff(fragResIdx, fixed).empty()) continue; // TERMs composed entirely of fixed residues have no impact
        cout << "TERM around " << C[ri] << endl;
        F.setRMSDCutoff(RMSDCalculator::rmsdCutoff(fragResIdx, S)); // account for spacing between residues from the same chain
        vector<Structure*> matches = getMatches(cache, frag, fragResIdx, numPerTERM, op.isGiven("s") ? centIdx : vector<int>());
        if (!matches.empty()) allMatches.push_back(matches);
      }
    }

    // then pair TERMs
    cout << "Searching for pair TERMs..." << endl;
    ConFind cfd(&RL, S);
    if (fixed.empty()) {
      L = cfd.getContacts(S, 0.01);
    } else {
      vector<Residue*> residues;
      for (int i = 0; i < contResis.size(); i++) residues.push_back(&(S.getResidue(contResis[i])));
      L = cfd.getContacts(residues, 0.01);
    }
    vector<pair<Residue*, Residue*> > contactList = L.getOrderedContacts();
    for (int k = 0; k < contactList.size(); k++) {
      Residue* resA = contactList[k].first;
      Residue* resB = contactList[k].second;
      Structure frag; vector<int> fragResIdx;
      vector<int> centIdx = TERMUtils::selectTERM({resA, resB}, frag, pmPair, &fragResIdx);
      if (MstUtils::setdiff(fragResIdx, fixed).empty()) continue; // TERMs composed entirely of fixed residues have no impact
      cout << "TERM around " << *resA << " x " << *resB << endl;
      F.setRMSDCutoff(RMSDCalculator::rmsdCutoff(fragResIdx, S)); // account for spacing between residues from the same chain
      vector<Structure*> matches = getMatches(cache, frag, fragResIdx, numPerTERM, op.isGiven("s") ? centIdx : vector<int>());
      if (!matches.empty()) allMatches.push_back(matches);
    }
    // fuser options
    fusionParams opts; opts.setNumIters(Ni); opts.setVerbose(false);
    opts.setMinimizerType(fusionParams::gradDescent);
    opts.setRepFC(1);
    opts.setCompFC(0.1);
    mstreal compactnessRadius = Rf;
    // (R0 - Rf)*exp(-c/10.0) + Rf; // exponential scaling
    // (Rf*(c + 1) + R0*(Ncyc - c - 1))/Ncyc; // linear scaling
    opts.setCompRad(compactnessRadius);
    cout << "will be trying to combine structure to a radius of " << compactnessRadius << endl;

    /* --- pick a random combination of TERMs and fuse --- */
    vector<vector<int> > currPicks(allMatches.size()), bestPicks;
    vector<vector<Residue*> > resTopo(I.residueSize());
    for (int si = 0; si < allMatches.size(); si++) {
      for (int j = 0; j < MstUtils::min(numPerTERM, (int) allMatches[si].size()); j++) {
        if (op.isGiven("r")) {
          currPicks[si].push_back(MstUtils::randInt(allMatches[si].size()));
        } else {
          currPicks[si].push_back(j);
        }
      }
    }

    // if there are fixed residues, then their starting conformation (and thus,
    // the final one also) is chosen from the original structure
    if (!fixed.empty()) opts.setStartingStructure(S);

    /* --- do an MC simulation to find a good combo of TERMs --- */
    fusionTopology bestTopo, currTopo;
    Structure bestFused, currFused;
    fusionScores bestScore, currScore, propScore;
    mstreal kT = 0.001;
    for (int it = 0; it < op.getInt("iter", 1); it++) {
      vector<vector<int> > propPicks = currPicks;
      // make a "mutation"
      if (it != 0) {
        int si = MstUtils::randInt(allMatches.size());
        int mi = MstUtils::randInt(propPicks[si].size());
        propPicks[si][mi] = MstUtils::randInt(allMatches[si].size());
      }
      fusionTopology propTopo = getTopo(I.residueSize(), allMatches, propPicks, (it == 0) ? shellOut : dummy, op.isGiven("m") ? &S : NULL);
      propTopo.addFixedPositions(fixed);
      Structure propFused;
      if (op.isGiven("a")) {
        vector<mstreal> ics = opts.getIntCoorFCs();
        opts.setIntCoorFCs({0, 0, 0});
        propFused = Fuser::fuse(propTopo, propScore, opts);
        opts.setIntCoorFCs(ics);
        opts.setStartingStructure(propFused);
        propFused = Fuser::fuse(propTopo, propScore, opts);
      } else {
        propFused = Fuser::fuse(propTopo, propScore, opts);
      }
      copySequence(S, propFused);
      AtomPointerVector init = getCorrespondingAtoms(S, propFused);
      cout << "\titeration " << it << " => "; totalScore(propScore, propFused, init, true); cout << endl;
      if ((it == 0) || mc(totalScore(currScore, currFused, init), totalScore(propScore, propFused, init), kT)) {
        cout << "\t\taccepted" << endl;
        currFused = propFused;
        currTopo = propTopo;
        currScore = propScore;
        currPicks = propPicks;
      }
      if ((it == 0) || (totalScore(propScore, propFused, init) < totalScore(bestScore, bestFused, init))) {
        cout << "\t\t\tnew best" << endl;
        bestFused = propFused;
        bestTopo = propTopo;
        bestScore = propScore;
        bestPicks = propPicks;
      }
    }

    // align based on the fixed part, if anything was fixed (for ease of visualization)
    if (fixed.size() > 0) {
      AtomPointerVector before = getBackbone(S, fixed);
      AtomPointerVector after = getBackbone(bestFused, fixed);
      rc.align(after, before, bestFused);
    }

    /* --- write intermediate result and clean up--- */
    S = bestFused.reassignChainsByConnectivity();
    for (int si = 0; si < allMatches.size(); si++) {
      for (int mi = 0; mi < allMatches[si].size(); mi++) {
        delete(allMatches[si][mi]);
      }
    }
    if (shellOut.is_open()) shellOut.close();

    out << "MODEL " << c+1 << endl;
    S.writePDB(out);
    out << "ENDMDL" << endl;
  }
  out.close();
  S.writePDB(op.getString("o") + ".fin.pdb");
}
