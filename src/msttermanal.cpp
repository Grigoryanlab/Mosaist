#include "mstmagic.h"
#include "mstrotlib.h"
#include "msttermanal.h"

mstreal TERMANAL::structureScore(const Structure& term, const vector<Residue*>& central, bool verbose) {
  pair<mstreal, mstreal> parts = structureScoreParts(term, central, verbose);
  return combineScoreParts(parts);
}

pair<mstreal, mstreal> TERMANAL::structureScoreParts(const Structure& term, const vector<Residue*>& central, bool verbose) {
  if (F == NULL) MstUtils::error("FASST object not set", "TERMANAL::structureScoreParts");
  if (verbose) {
    cout << "\tvisiting TERM (" << MstUtils::vecPtrToString(central) << "): " << MstUtils::vecPtrToString(term.getResidues()) << endl;
    // term.writePDB(cout);
  }

  // get the top hits
  fasstSearchOptions origOpts = setupSearch(term);
  fasstSolutionSet matches = F->search();
  vector<fasstSolution*> topMatches;
  vector<mstreal> rmsds;
  if (compatMode) {
    vector<Atom*> termBB = RotamerLibrary::getBackbone(term);
    topMatches = getTopMatches(F, matches, termBB, &rmsds);
  } else {
    int numMatches = matches.size();
    topMatches.resize(numMatches);
    rmsds.resize(numMatches);
    for (int i = 0; i < numMatches; i++) {
      topMatches[i] = &matches[i];
      rmsds[i] = matches[i].getRMSD();
    }
  }

  // sequence likelihood = fraction of top hits that share the right amino acid
  mstreal seqLikeComb = calcSeqLikelihood(topMatches, central, verbose);

  // structure frequency = number of matches with RMSD below that of the Nth
  // match to the top native representative of the query, divided by N=matchCount
  mstreal structFreq = calcStructFreq(topMatches, rmsds, verbose);

  // recover original search options
  F->setOptions(origOpts);

  return make_pair(seqLikeComb, structFreq);
}

mstreal TERMANAL::combineScoreParts(const pair<mstreal, mstreal>& parts) {
  return log(parts.first * parts.second + pseudoCount);
}

vector<mstreal> TERMANAL::scoreStructure(const Structure& S, const vector<Residue*>& subregion, vector<pair<mstreal, mstreal>>* scoreParts, bool verbose) {
  if (!RL.isLoaded()) MstUtils::error("Rotamer library not loaded", "TERMANAL::scoreStructure");

  // This map will store the score components for various visited TERMs that. The
  // key is the residue index of the central residue of the corresponding TERM.
  map<int, pair<mstreal, mstreal>> structScoreParts;

  // NOTE: we need to identify all TERMs, even if we are interested in scores at
  // a subset of residues. This is because the final smoothed score of a residue
  // depends on all TERMs that it is a part of, whether the central residues of
  // those TERMs happen to be in our subset or not.
  // Variable resOverlaps keeps track of which TERMs each residue belongs to,
  // index by full residue index into S.
  ConFind C(&RL, S);
  vector<Residue*> centrals;
  vector<vector<int>> resOverlaps;
  vector<Structure> terms = collectTERMs(C, S.getResidues(), centrals, resOverlaps);

  // Smooth and combine each pair of score parts by incorporating the scores from each TERM that a residue belongs to
  vector<mstreal> structScores(subregion.size());
  if (scoreParts != NULL) scoreParts->resize(subregion.size());
  for (int i = 0; i < subregion.size(); i++) {
    int idx = subregion[i]->getResidueIndex();
    if (verbose) cout << "Scoring " << *(subregion[i]) << " ..." << endl;

    // visit any TERMs necessary for computing the scores at residue index idx
    for (int ti : resOverlaps[idx]) {
      if (structScoreParts.find(ti) == structScoreParts.end()) {
        structScoreParts[ti] = structureScoreParts(terms[ti], {centrals[ti]}, verbose);
      }
    }

    // apply smoothing to score components
    pair<mstreal, mstreal> parts = smoothScores(structScoreParts, resOverlaps[idx]);

    // compute and output total score
    structScores[i] = combineScoreParts(parts);
    mstreal designScore = log(parts.first + pseudoCount), abundanceScore = log(parts.second + pseudoCount);
    if (scoreParts != NULL) (*scoreParts)[i] = make_pair(designScore, abundanceScore);
    if (verbose) cout << "scores for " << *(subregion[i]) << ": " << *(subregion[i]) << "\t" << designScore << "\t" << abundanceScore << "\t" << structScores[i] << endl;
  }
  return structScores;
}

fasstSearchOptions TERMANAL::setupSearch(const Structure& S) {
  fasstSearchOptions origOpts = F->options(); // in case the same FASST object is being shared by others
  F->setOptions(fasstSearchOptions());
  F->options().setRedundancyCut(0.7);
  F->options().setRMSDCutoff(rmsdCut);
  F->options().setMinNumMatches(0);
  int maxNumMatches = compatMode ? compatSearchLimit : matchCount;
  F->options().setMaxNumMatches(maxNumMatches);
  F->setQuery(S);
  return origOpts;
}

vector<fasstSolution*> TERMANAL::getTopMatches(FASST* F, fasstSolutionSet& matches, vector<Atom*>& queryA, vector<mstreal>* topRmsds) {
  RMSDCalculator RC;
  int numMatches = matches.size();
  vector<pair<mstreal, int>> rmsds(numMatches);
  for (int i = 0; i < numMatches; i++) {
    Structure matchS = F->getMatchStructure(matches[i]);
    AtomPointerVector matchA = matchS.getAtoms();
    rmsds[i] = make_pair(RC.bestRMSD(queryA, matchA), i);
  }
  sort(rmsds.begin(), rmsds.end());
  int numToKeep = min(matchCount, numMatches);
  vector<fasstSolution*> topMatches(numToKeep);
  if (topRmsds != NULL) topRmsds->resize(numToKeep);
  for (int i = 0; i < numToKeep; i++) {
    topMatches[i] = &matches[rmsds[i].second];
    if (topRmsds != NULL) (*topRmsds)[i] = rmsds[i].first;
  }
  return topMatches;
}

mstreal TERMANAL::calcSeqLikelihood(vector<fasstSolution*>& matches, const vector<Residue*>& central, bool verbose) {
  int numMatches = matches.size();
  vector<Sequence> matchSeqs(numMatches);
  for (int i = 0; i < numMatches; i++) matchSeqs[i] = F->getMatchSequence(*matches[i]);
  vector<mstreal> seqLike(central.size(), 0);
  for (int i = 0; i < central.size(); i++) seqLike[i] = calcSeqFreq(central[i], matchSeqs);
  mstreal seqLikeComb = !central.empty() ? CartesianPoint(seqLike).mean() : 1.0;
  if (verbose) cout << "\tsequence likelihood = [" << MstUtils::vecToString(seqLike) << "], overall = " << seqLikeComb << endl;
  return seqLikeComb;
}

mstreal TERMANAL::calcSeqFreq(Residue* res, vector<Sequence>& matchSeqs) {
  int idx = res->getResidueIndex();
  res_t aa = SeqTools::aaToIdx(res->getName());
  if (aa == SeqTools::unknownIdx()) MstUtils::error("unknown residue name in central residue " + MstUtils::toString(*res) + ", which is a problem for computing sequence likelihood", "TERMANAL::calcSeqFreq");
  int seqCount = 0;
  for (int i = 0; i < matchSeqs.size(); i++) if (matchSeqs[i][idx] == aa) ++seqCount;
  mstreal seqFreq = 1.0*seqCount / matchCount;
  return seqFreq;
}

mstreal TERMANAL::calcStructFreq(vector<fasstSolution*>& matches, vector<mstreal>& rmsds, bool verbose) {
  int n;
  if (matches.size() >= 1) {
    Structure firstMatch = F->getMatchStructure(*matches[0]);
    F->setQuery(firstMatch);
    fasstSolutionSet matchesToClosestNative = F->search();
    mstreal r;
    if (compatMode) {
      vector<Atom*> firstMatchBB = RotamerLibrary::getBackbone(firstMatch);
      vector<mstreal> nativeRmsds;
      getTopMatches(F, matchesToClosestNative, firstMatchBB, &nativeRmsds);
      r = nativeRmsds.back();
    } else r = matchesToClosestNative.worstRMSD();
    for (n = 0; n < matches.size(); n++) {
      if (rmsds[n] > r) break;
    }
  } else n = 0;
  mstreal structFreq = min(1.0, (1.0*n)/matchCount);
  if (verbose) cout << "\tstructure frequency = " << structFreq << endl;
  return structFreq;
}

vector<Structure> TERMANAL::collectTERMs(ConFind& C, const vector<Residue*>& subregion, vector<Residue*>& centrals, vector<vector<int>>& resOverlaps) {
  int numTerms = subregion.size();
  vector<Structure> terms(numTerms);
  Structure* S = subregion[0]->getStructure();
  centrals.resize(S->residueSize());
  resOverlaps.resize(S->residueSize());
  for (int i = 0; i < numTerms; i++) {
    vector<int> fragResIdx;
    int resInd = TERMUtils::selectTERM(*(subregion[i]), C, terms[i], pad, cdCut, &fragResIdx);
    centrals[i] = &terms[i].getResidue(resInd);
    for (int j = 0; j < fragResIdx.size(); j++) resOverlaps[fragResIdx[j]].push_back(i);
  }
  return terms;
}

pair<mstreal, mstreal> TERMANAL::smoothScores(map<int, pair<mstreal, mstreal>>& structScoreParts, vector<int>& overlapResInds) {
  mstreal seqLike = 0.0, structFreq = 0.0;
  for (int i : overlapResInds) {
    if (structScoreParts.find(i) == structScoreParts.end()) MstUtils::error("not all required structure score parts are precalculated (ex: residue index " + MstUtils::toString(i) + ")", "TERMANAL::smoothScores");
    seqLike += structScoreParts[i].first;
    structFreq += structScoreParts[i].second;
  }
  int numOverlaps = overlapResInds.size();
  seqLike /= numOverlaps;
  structFreq /= numOverlaps;
  return make_pair(seqLike, structFreq);
}
