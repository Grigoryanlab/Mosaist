#include "mstmagic.h"
#include "msttermanal.h"

mstreal TERMANAL::structureScore(const Structure& S, const vector<Residue*>& central, bool verbose) {
  pair<mstreal, mstreal> parts = structureScoreParts(S, central, verbose);
  return combineScoreParts(parts);
}

pair<mstreal, mstreal> TERMANAL::structureScoreParts(const Structure& S, const vector<Residue*>& central, bool verbose) {
  if (F == NULL) MstUtils::error("FASST object not set", "TERMANAL::structureScoreParts");

  // get the top hits for S
  fasstSearchOptions origOpts = setupSearch(S);
  fasstSolutionSet matches = F->search();
  vector<Sequence> matchSeqs = F->getMatchSequences(matches);
  if (matchSeqs.size() != matchCount) MstUtils::error("Not enough matches returned for TERM '" + S.getName() + "', which is a problem for computing sequence likelihood", "TERMANAL::structureScore");

  // sequence likelihood = fraction of top hits that share the right amino acid
  mstreal seqLikeComb = calcSeqLikelihood(central, matchSeqs, verbose);

  // structure frequency = number of matches with RMSD below that of the Nth
  // match to the top native representative of the query, divided by N=matchCount
  mstreal structFreq = calcStructFreq(matches, verbose);

  // recover original search options
  F->setOptions(origOpts);

  return make_pair(seqLikeComb, structFreq);
}

mstreal TERMANAL::combineScoreParts(const pair<mstreal, mstreal>& parts) {
  return -log(parts.first * parts.second + pseudoCount);
}

vector<mstreal> TERMANAL::scoreStructure(const Structure& S, const vector<Residue*>& subregion, bool verbose) {
  if (!RL.isLoaded()) MstUtils::error("Rotamer library not loaded", "TERMANAL::scoreStructure");

  // Create a TERM for each residsue and keep track of which TERMs each residue belongs to via 'resOverlaps'
  ConFind C(&RL, S);
  vector<vector<Residue*>> termCentrals;
  vector<vector<int>> resOverlaps;
  vector<Structure> terms = collectTERMs(C, subregion, termCentrals, resOverlaps);

  // Calculate the score parts for each TERM
  int numTerms = terms.size();
  vector<pair<mstreal, mstreal>> structScoreParts(numTerms);
  for (int i = 0; i < numTerms; i++) structScoreParts[i] = structureScoreParts(terms[i], termCentrals[i], verbose);

  // Smooth and combine each pair of score parts by incorporating the scores from each TERM that a residue belongs to
  vector<mstreal> structScores(numTerms);
  for (int i = 0; i < numTerms; i++) structScores[i] = smoothScore(structScoreParts, i, resOverlaps[i]);
  return structScores;
}

fasstSearchOptions TERMANAL::setupSearch(const Structure& S) {
  fasstSearchOptions origOpts = F->options(); // in case the same FASST object is being shared by others
  F->setOptions(fasstSearchOptions());
  F->options().setRedundancyCut(0.6);
  F->options().setRMSDCutoff(RMSDCalculator::rmsdCutoff(S, 1.1, 15));
  F->options().setMinNumMatches(matchCount);
  F->options().setMaxNumMatches(matchCount);
  F->setQuery(S);
  return origOpts;
}

mstreal TERMANAL::calcSeqLikelihood(const vector<Residue*>& central, vector<Sequence>& matchSeqs, bool verbose) {
  vector<mstreal> seqLike(central.size(), 0);
  for (int i = 0; i < central.size(); i++) seqLike[i] = calcSeqFreq(central[i], matchSeqs);
  mstreal seqLikeComb = !central.empty() ? CartesianPoint(seqLike).mean() : 1.0;
  if (verbose) cout << "sequence likelihood = [" << MstUtils::vecToString(seqLike) << "], overall = " << seqLikeComb << endl;
  return seqLikeComb;
}

mstreal TERMANAL::calcSeqFreq(Residue* res, vector<Sequence>& matchSeqs) {
  int idx = res->getResidueIndex();
  res_t aa = SeqTools::aaToIdx(res->getName());
  if (aa == SeqTools::unknownIdx()) MstUtils::error("unknown residue name in central residue " + MstUtils::toString(*res) + ", which is a problem for computing sequence likelihood", "TERMANAL::calcSeqFreq");
  mstreal seqFreq = 0.0;
  for (int i = 0; i < matchSeqs.size(); i++) if (matchSeqs[i][idx] == aa) seqFreq += 1.0;
  seqFreq /= matchSeqs.size();
  return seqFreq;
}

mstreal TERMANAL::calcStructFreq(fasstSolutionSet& matches, bool verbose) {
  F->setQuery(F->getMatchStructure(matches[0]));
  fasstSolutionSet matchesToClosestNative = F->search();
  mstreal r = matchesToClosestNative[matchCount - 1].getRMSD();
  int n;
  for (n = 0; n < matches.size(); n++) {
    if (matches[n].getRMSD() > r) break;
  }
  mstreal structFreq = min(1.0, (1.0*n)/matchCount);
  if (verbose) cout << "structure frequency = " << structFreq << endl;
  return structFreq;
}

vector<Structure> TERMANAL::collectTERMs(ConFind& C, const vector<Residue*>& subregion, vector<vector<Residue*>>& termCentrals, vector<vector<int>>& resOverlaps) {
  int numTerms = subregion.size();
  vector<Structure> terms(numTerms);
  termCentrals.resize(numTerms);
  Structure* S = subregion[0]->getStructure();
  resOverlaps.resize(S->residueSize());
  for (int i = 0; i < numTerms; i++) {
    vector<int> fragResIdx;
    int resInd = TERMUtils::selectTERM(*subregion[i], C, terms[i], pad, cdCut, &fragResIdx);
    Residue* termCentRes = &terms[i].getResidue(resInd);
    termCentrals[i] = {termCentRes};
    for (int j = 0; j < fragResIdx.size(); j++) resOverlaps[i].push_back(fragResIdx[j]);
  }
  return terms;
}

mstreal TERMANAL::smoothScore(vector<pair<mstreal, mstreal>>& structScoreParts, int resInd, vector<int>& overlapResInds) {
  mstreal seqLike = structScoreParts[resInd].first, structFreq = structScoreParts[resInd].second;
  for (int i : overlapResInds) {
    seqLike += structScoreParts[i].first;
    structFreq += structScoreParts[i].second;
  }
  int numOverlaps = overlapResInds.size();
  seqLike /= numOverlaps;
  structFreq /= numOverlaps;
  return combineScoreParts(make_pair(seqLike, structFreq));
}
