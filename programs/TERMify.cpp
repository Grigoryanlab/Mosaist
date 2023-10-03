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

void numberResidues(Structure& S, const vector<int>& resIdx) {
  for (int k = 0; k < S.residueSize(); k++) S.getResidue(k).setNum(resIdx[k]);
}

vector<Structure*> getMatches(FASST& C, Structure& frag, const vector<int>& fragResIdx, int need = 5, const vector<int>& centIdx = vector<int>()) {
  vector<Structure*> matchStructures;
  if (need == 0) return matchStructures;
  C.setQuery(frag, false);
  C.setRMSDCutoff(RMSDCalculator::rmsdCutoff(frag));
  C.options().setMaxNumMatches(need);
  C.options().setMinNumMatches(need);
  C.options().unsetSequenceConstraints();

  // add sequence constraints, as needed
  fasstSeqConstSimple seqConst(centIdx.size());
  Structure splitQuery = C.getQuery();
  for (int i = 0; i < centIdx.size(); i++) {
    if (!SeqTools::isUnknown(frag.getResidue(centIdx[i]).getName())) {
      const Residue& res = splitQuery.getResidue(centIdx[i]);
      seqConst.addConstraint(res.getChain()->getIndex(), res.getResidueIndexInChain(), {res.getName()});
    }
  }
  if (seqConst.hasConstraints()) C.options().setSequenceConstraints(seqConst);

  // limit iterations, because it is technically possible that a match meeting
  // the sequence constraints does not exist, at which point we will just give up
  while (true) {
    fasstSolutionSet matches = C.search();
    for (auto it = matches.begin(); (it != matches.end()) && (matchStructures.size() != need); ++it) {
      matchStructures.push_back(new Structure(C.getMatchStructure(*it, false, FASST::matchType::REGION)));
      Structure& match = *(matchStructures.back());
      if (!RotamerLibrary::hasFullBackbone(match)) {
        delete(matchStructures.back()); matchStructures.pop_back();
        continue;
      }
      RotamerLibrary::standardizeBackboneNames(match);
      MstUtils::assertCond(match.residueSize() == fragResIdx.size(), "unexpected match size");
      numberResidues(match, fragResIdx); // make residue numbers store indices into the original structure
    }
    if (C.isVerbose()) cout << "\tfound " << matchStructures.size() << " matches" << endl;
    if (matchStructures.size() == need) break;

    cout << "\t\tneed to search again..." << endl;
    for (Structure* m : matchStructures) delete(m);
    matchStructures.clear();
    int newNeeded = (int) ceil(1.5 * C.options().getMaxNumMatches());
    C.options().setMaxNumMatches(newNeeded);
    C.options().setMinNumMatches(newNeeded);
  }

  return matchStructures;
}

bool mc(mstreal oldScore, mstreal newScore, mstreal kT) {
  return ((newScore < oldScore) || (MstUtils::randUnit() < exp((oldScore - newScore)/kT)));
}

fusionTopology getTopo(int L, vector<vector<Structure*> >& allMatches, vector<int>& picks, fstream& matchOut, Structure* global = NULL) {
  fusionTopology resTopo(L);
  for (int si = 0; si < allMatches.size(); si++) {
    Structure& match = *(allMatches[si][picks[si]]);
    resTopo.addFragment(match);
    if (matchOut.is_open()) match.writePDB(matchOut);
  }
  if (global != NULL) {
    MstUtils::assertCond(global->residueSize() == L, "the global target structure specified has an unexpected number of residues for the topology");
    resTopo.addFragment(*global, MstUtils::range(0, L), -10.0);
    if (matchOut.is_open()) global->writePDB(matchOut);
  }
  return resTopo;
}

fusionTopology getTopo(int L, vector<vector<Structure*> >& allMatches, vector<vector<int> >& picks, fstream& matchOut, Structure* global = NULL) {
  fusionTopology resTopo(L);
  for (int si = 0; si < allMatches.size(); si++) {
    for (int j = 0; j < picks[si].size(); j++) {
      Structure& match = *(allMatches[si][picks[si][j]]);
      resTopo.addFragment(match);
      if (matchOut.is_open()) match.writePDB(matchOut);
    }
  }
  if (global != NULL) {
    MstUtils::assertCond(global->residueSize() == L, "the global target structure specified has an unexpected number of residues for the topology");
    resTopo.addFragment(*global, MstUtils::range(0, L), -10.0);
    if (matchOut.is_open()) global->writePDB(matchOut);
  }
  return resTopo;
}

AtomPointerVector getCorrespondingAtoms(Structure& from, Structure& like) {
  MstUtils::assertCond(from.residueSize() == like.residueSize(), "the two structures must have the same number of residues", "getCorrespondingAtoms()");
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

mstreal totalScore(fusionOutput& scoreObj, Structure& fused, AtomPointerVector& init, bool report = false) {
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

bool sameSize(Structure& A, Structure& B) {
  if (A.chainSize() != B.chainSize()) return false;
  for (int i = 0; i < A.chainSize(); i++) {
    if (A[i].residueSize() != B[i].residueSize()) return false;
  }
  return true;
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
  op.addOption("frag", "a selection that will be used to define a single \"TERM\" to include in the fuser topology. This can be used to restrain some portion of the structure durring minimization/dynamics. This option can be given multiple times and will result in multiple suchh TERMs being defined.");
  op.addOption("us", "a selection string for residues to mark as having unknown identity (i.e., their identity will not matter if accounting for sequence).");
  op.addOption("rad", "compactness radius. Default will be based on protein length.");
  op.addOption("c", "path to a FASST cache file to use for initializing the cache.");
  op.addOption("w", "flag; if specified, the FASST cache will be periodically updated.");
  op.addOption("app", "flag; if specified, will append to the output PDB file (e.g., for the purpose of accumulating a trajectory from multiple runs).");
  op.addOption("dyn", "use dynamics rather than optimization to search for a solution. If a number is specified, it is interpreted as the length of the dynamics simulation (relative to the length of a typical minimization run); default is 100.");
  op.addOption("alt", "alternative conformation PDB file (must have the same length/topology as the starting conformation file). If given, for each TERM the corresponding segment of structure will be added as a \"match\".");
  op.addOption("orig", "If given, each TERM's original conformation (from the starting structure) will be explicitly added as a \"match\".");
  op.addOption("v", "set verbose output flag.");
  op.addOption("cycCheck","flag; if given, will check whether the fused structure has converged and will potentially quick early. Convergence is established by comparing the RMSD resultant from the current cycle to the average RMSD from the first 10 cycles. If the latter is less than a third of the former, the cycling is said to have converged.");
  if (op.isGiven("f") && op.isGiven("fs")) MstUtils::error("only one of --f or --fs can be given!");

  op.setOptions(argc, argv);
  RMSDCalculator rc;
  Structure I(op.getString("p")), A;
  vector<int> fixed;
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

  // define any custom TERMs
  vector<Structure> customTERMs(op.timesGiven("frag"));
  for (int fi = 0; fi < op.timesGiven("frag"); fi++) {
    selector sel(I);
    vector<Residue*> customResidues = sel.selectRes(op.getString("frag", "", fi));
    cout << "defining custom TERM " << fi << " with " << customResidues.size() << " residues..." << endl;
    customTERMs[fi] = Structure(customResidues);
    vector<Residue*> newResis = customTERMs[fi].getResidues();
    for (int ri = 0; ri < newResis.size(); ri++) newResis[ri]->setNum(customResidues[ri]->getResidueIndex());
  }

  // figure out which FASST object will be used
  if (op.isGiven("s") && (op.isGiven("c") || op.isGiven("w"))) MstUtils::error("cannot specify --s with caching");
  fasstCache withCache((I.residueSize() - fixed.size())*10); FASST woCache;
  string tag = "TERMify-" + MstSys::getUserName();
  bool useCache = op.isGiven("c") || op.isGiven("w") || !op.isGiven("s");
  FASST& search = useCache ? withCache : woCache;
  if (useCache) {
    withCache.setStrictEquivalence(true);
    withCache.setVerbose(op.isGiven("v"));
    if (op.isGiven("c") && op.getString("c").empty()) MstUtils::error("--c must be a valid file path");
    // read initial cache
    if (op.isGiven("c")) {
      if (op.getString("c").empty()) MstUtils::error("--c must be a valid file path");
      if (MstSys::fileExists(op.getString("c"))) MstUtils::error("--c is not an existing file");
      MstSys::getNetLock(tag, true);
      cout << "reading cache from " << op.getString("c") << "... " << endl;
      withCache.read(op.getString("c"));
      MstSys::releaseNetLock(tag);
    }
  }
  if (op.isGiven("d")) {
    search.addTargets(MstUtils::fileToArray(op.getString("d")));
    if (op.isGiven("b")) {
      search.writeDatabase(op.getString("b"));
    }
  } else if (op.isGiven("b")) {
    search.readDatabase(op.getString("b"), 2);
  } else {
    MstUtils::error("either --b or --d must be given!");
  }

  if (op.isGiven("us")) {
    selector sel(I);
    vector<Residue*> unkResidues = sel.selectRes(op.getString("us"));
    cout << "unknown sequence selection gave " << unkResidues.size() << " residues, marking unknown..." << endl;
    for (int i = 0; i < unkResidues.size(); i++) unkResidues[i]->setName("UNK");
  }
  if (op.isGiven("alt")) {
    A.readPDB(op.getString("alt"));
    if (!sameSize(A, I)) MstUtils::error(op.getString("alt") + " and " + op.getString("p") + " do not have the same topology");
  }
  int numPerTERM = op.getInt("n", 1);

  if (search.isResidueRelationshipPopulated("sim")) search.setRedundancyProperty("sim");
  else search.setRedundancyCut(0.5);
  RotamerLibrary RL(op.getString("rLib"));
  int pmSelf = 2, pmPair = 1;
  int Ni = 1000, lastWriteTime;
  fusionOutput bestScore, currScore;
  fstream out, shellOut, dummy;
  contactList L;

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

  // TERMify loop
  // NOTE: we reassign by connectivity at the start so that we do not have to
  // worry about connectivity during the cycling (things may get temporarily
  // broken as TERMs fight each other, but we will always interpret topology as
  // connectivity past this step).
  Structure S = I.reassignChainsByConnectivity(); // fixed residues are already selected, but this does not change residue order
  RotamerLibrary::standardizeBackboneNames(S);
  vector<Structure*> O = {op.isGiven("alt") ? &A : NULL, op.isGiven("orig") ? &S : NULL}; // any "other" structures from which we should get matches
  mstreal R0 = getRadius(I);
  mstreal Rf = op.getReal("rad", pow(I.residueSize() * 1.0, 1.0/3)*5.0);
  int Ncyc = op.getInt("cyc", 10);
  MstUtils::openFile(out, op.getString("o") + ".traj.pdb", op.isGiven("app") ? ios::app : ios::out);
  out << "MODEL " << 0 << endl; S.writePDB(out); out << "ENDMDL" << endl;
  //create a vector of RMSD_final of each cycle
  vector<mstreal> cyc_rmsd;
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
        withCache.write(op.getString("c"));
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
        vector<int> centIdx = TERMUtils::selectTERM({&C[ri]}, frag, pmSelf, &fragResIdx, false);
        if (MstUtils::setdiff(fragResIdx, fixed).empty()) continue; // TERMs composed entirely of fixed residues have no impact
        cout << "TERM around " << C[ri] << endl;
        search.setRMSDCutoff(RMSDCalculator::rmsdCutoff(fragResIdx, S)); // account for spacing between residues from the same chain
        vector<Structure*> matches = getMatches(search, frag, fragResIdx, numPerTERM, op.isGiven("s") ? centIdx : vector<int>());
        for (int ii = 0; ii < O.size(); ii++) {
          if (O[ii] == NULL) continue;
          Structure* altFrag = new Structure();
          TERMUtils::selectTERM({&(O[ii]->getResidue(C[ri].getResidueIndex()))}, *altFrag, pmSelf, NULL, false);
          numberResidues(*altFrag, fragResIdx);
          matches.push_back(altFrag);
        }
        if (!matches.empty()) {
          allMatches.push_back(matches);
          int lastIdx = MstUtils::min(MstUtils::max(numPerTERM, 1), (int) matches.size());
          Structure* last = matches[lastIdx - 1];
          if (op.isGiven("v")) cout << "\tRMSD of match " << lastIdx << " is " << rc.bestRMSD(RotamerLibrary::getBackbone(*last), RotamerLibrary::getBackbone(frag)) << endl;
        }
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
      vector<int> centIdx = TERMUtils::selectTERM({resA, resB}, frag, pmPair, &fragResIdx, false);
      if (MstUtils::setdiff(fragResIdx, fixed).empty()) continue; // TERMs composed entirely of fixed residues have no impact
      cout << "TERM around " << *resA << " x " << *resB << endl;
      search.setRMSDCutoff(RMSDCalculator::rmsdCutoff(fragResIdx, S)); // account for spacing between residues from the same chain
      vector<Structure*> matches = getMatches(search, frag, fragResIdx, numPerTERM, op.isGiven("s") ? centIdx : vector<int>());
      for (int ii = 0; ii < O.size(); ii++) {
        if (O[ii] == NULL) continue;
        Structure* altFrag = new Structure();
        TERMUtils::selectTERM({&(O[ii]->getResidue(resA->getResidueIndex())), &(O[ii]->getResidue(resB->getResidueIndex()))}, *altFrag, pmPair, NULL, false);
        numberResidues(*altFrag, fragResIdx);
        matches.push_back(altFrag);
      }
      if (!matches.empty()) {
        allMatches.push_back(matches);
        int lastIdx = MstUtils::min(MstUtils::max(numPerTERM, 1), (int) matches.size());
        Structure* last = matches[lastIdx - 1];
        if (op.isGiven("v")) cout << "\tRMSD of match " << lastIdx << " is " << rc.bestRMSD(RotamerLibrary::getBackbone(*last), RotamerLibrary::getBackbone(frag)) << endl;
      }
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
      for (int j = 0; j < allMatches[si].size(); j++) {
        if (op.isGiven("r")) {
          currPicks[si].push_back(MstUtils::randInt(allMatches[si].size()));
          break;
        } else {
          currPicks[si].push_back(j);
        }
      }
    }


    // if (!fixed.empty()) opts.setStartingStructure(S); // if there are fixed residues, then their conformation is taken from the original structure
    opts.setStartingStructure(S);

    /* --- do an MC simulation to find a good combo of TERMs --- */
    fusionTopology bestTopo, currTopo;
    Structure bestFused, currFused;
    fusionOutput bestScore, currScore, propScore;
    mstreal kT = 0.001;
    mstreal rmsd_final;
    for (int it = 0; it < op.getInt("iter", 1); it++) {
      vector<vector<int> > propPicks = currPicks;
      // make a "mutation"
      if (it != 0) {
        int si = MstUtils::randInt(allMatches.size());
        int mi = MstUtils::randInt(propPicks[si].size());
        propPicks[si][mi] = MstUtils::randInt(allMatches[si].size());
      }
      fusionTopology propTopo = getTopo(I.residueSize(), allMatches, propPicks, (it == 0) ? shellOut : dummy, op.isGiven("m") ? &S : NULL);
      for (int fi = 0; fi < customTERMs.size(); fi++) propTopo.addFragment(customTERMs[fi]);
      propTopo.addFixedPositions(fixed);
      Structure propFused;
      if (op.isGiven("dyn")) {
        opts.setMinimizerType(fusionParams::langevinDyna);
        opts.setNumIters((op.isInt("dyn") ? op.getInt("dyn") : 100)*Ni);
        opts.setLogBase(op.getString("o"));
        opts.setThermalEnergy(1.0);
        opts.setSaveInterval(MstUtils::max(1, opts.numIters()/1000));
        opts.setAdaptiveWeighting(true);
        propFused = Fuser::fuse(propTopo, propScore, opts);
        opts.setMinimizerType(fusionParams::gradDescent);
        opts.setNumIters(Ni);
        opts.setStartingStructure(propFused);
      }
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
        if (op.isGiven("v")) cout << "\t\taccepted" << endl;
        currFused = propFused;
        currTopo = propTopo;
        currScore = propScore;
        currPicks = propPicks;
      }
      if ((it == 0) || (totalScore(propScore, propFused, init) < totalScore(bestScore, bestFused, init))) {
        if (op.isGiven("v")) cout << "\t\t\tnew best" << endl;
        bestFused = propFused;
        bestTopo = propTopo;
        bestScore = propScore;
        bestPicks = propPicks;
      }

      // calculate the RMSD of this cycle for cycCheck before S gets updated below
      RMSDCalculator checkPoint;
      cyc_rmsd.push_back(checkPoint.bestRMSD(init, bestFused.getAtoms()));
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
    Structure::combine(S, I).writePDB(out);
    out << "ENDMDL" << endl;
    // if rmsd is less than one third of the rmsd average of the first 10 cycles, break off the termify cycle--
    int initCycs = 10;
    if (op.isGiven("cycCheck") && (c >= initCycs)) {
      //cutoff = (the average of the RMSDs for the first initCycs cycles)/3
      mstreal r_tot = 0;
      for (int r = 0; r < initCycs; r++) r_tot += cyc_rmsd[r];
      cout << "Cycle check cutoff: " << (r_tot/initCycs)/3 << endl;
      if (cyc_rmsd.back() < (r_tot/initCycs)/3) break;
    }
  }
  out.close();
  Structure::combine(S, I).writePDB(op.getString("o") + ".fin.pdb");
}
