#include "msttermanal.h"

// Returns a map from each residue in S to its set of corresponding residues in each Structure in otherS
// Residues with no corresponding residues in other structures do not show up as keys in the map
// This is needed to perform the smoothing done in the paper's methods
map<Residue*, vector<Residue*>> TERMANAL::getOverlappingTerms(const Structure& S, const vector<Structure*>& otherS) {
  return getOverlappingTerms(S.getResidues(), otherS);
}

map<Residue*, vector<Residue*>> TERMANAL::getOverlappingTerms(const vector<Residue*>& R, const vector<Structure*>& otherS) {
  map<Residue*, vector<Residue*>> overlapSets;
  for (const Structure* Si: otherS) {
    map<Residue*, Residue*> overlaps = getOverlaps(R, *Si);
    for (auto j = overlaps.begin(); j != overlaps.end(); ++j) overlapSets[j->first].push_back(j->second);
  }
  return overlapSets;
}

// Returns a map from each residue in R to its corresponding residue in otherS (if such a residue exists)
map<Residue*, Residue*> TERMANAL::getOverlaps(const vector<Residue*>& R, const Structure& otherS) {
  map<pair<char, int>, Residue*> resMap1 = mapResInds(R), resMap2 = mapResInds(otherS.getResidues());
  map<Residue*, Residue*> overlaps;
  for (auto i = resMap1.begin(); i != resMap1.end(); ++i) {
    const pair<char, int>& id = i->first;
    Residue* res = i->second;
    if (resMap2.find(id) != resMap2.end()) overlaps[res] = resMap2[id];
  }
  return overlaps;
}

// Creates a map from each residue's unique ID (chain, residue index) to its pointer in the Structure
map<pair<char, int>, Residue*> TERMANAL::mapResInds(const vector<Residue*>& R) {
  map<pair<char, int>, Residue*> resMap;
  for (Residue* res: R) {
    char chainID = res->getChainID().empty() ? ' ' : res->getChainID()[0];
    resMap[make_pair(chainID, res->getNum())] = res;
  }
  return resMap;
}

mstreal TERMANAL::structureScore(const Structure& S, const vector<Residue*>& central, const map<Residue*, vector<Residue*>>& overlapSets, bool verbose) {
  if (F == NULL) MstUtils::error("FASST object not set", "TERMANAL::structureScore");
  fasstSearchOptions origOpts = F->options(); // in case the same FASST object is being shared by others
  F->setOptions(fasstSearchOptions());
  F->options().setRedundancyCut(0.6);

  // find top hits for each overlapping central residue
  int N = 50;
  F->options().setRMSDCutoff(RMSDCalculator::rmsdCutoff(S, 1.1, 15));
  F->options().setMinNumMatches(N);
  F->options().setMaxNumMatches(N);

  // sequence likelihood = fraction of top hits that share the right amino acid
  vector<mstreal> seqLike(central.size(), 0);
  map<Structure*, fasstSolutionSet> closestMatches;
  for (int i = 0; i < central.size(); i++) {
    Residue* res = central[i];
    vector<Residue*> overlaps({res});
    const auto& overlapSet = overlapSets.find(res);
    if (overlapSet != overlapSets.end()) overlaps.insert(overlaps.end(), overlapSet->second.begin(), overlapSet->second.end());
    for (int j = 0; j < overlaps.size(); j++) seqLike[i] += calcSeqFreq(overlaps[j], closestMatches);
    seqLike[i] /= overlaps.size();
  }
  mstreal seqLikeComb = (!central.empty()) ? CartesianPoint(seqLike).mean() : 1.0;
  if (verbose) cout << "sequence likelihood = [" << MstUtils::vecToString(seqLike) << "], overall = " << seqLikeComb << endl;

  // structure frequency = number of matches with RMSD below that of the Nth
  // match to the top native representative of the query, divided by N
  int numS = closestMatches.size();
  mstreal structFreq = 0.0;
  for (auto i = closestMatches.begin(); i != closestMatches.end(); ++i) {
    fasstSolutionSet& matches = i->second;
    F->setQuery(F->getMatchStructure(matches[0]));
    fasstSolutionSet matchesToClosestNative = F->search();
    mstreal r = matchesToClosestNative[N - 1].getRMSD();
    int n;
    for (n = 0; n < matches.size(); n++) {
      if (matches[n].getRMSD() > r) break;
    }
    structFreq += min(1.0, (1.0*n)/N);
  }
  structFreq /= numS;
  if (verbose) cout << "structure frequency = " << structFreq << endl;

  // structure score combines the two
  mstreal eps = 0.01;
  mstreal ss = -log(structFreq * seqLikeComb + eps);

  // recover original search options
  F->setOptions(origOpts);

  return ss;
}

mstreal TERMANAL::calcSeqFreq(Residue* res, map<Structure*, fasstSolutionSet>& closestMatches) {
  int idx = res->getResidueIndex();
  res_t aa = SeqTools::aaToIdx(res->getName());
  if (aa == SeqTools::unknownIdx()) MstUtils::error("unknown residue name in central residue " + MstUtils::toString(*res) + ", which is a problem for computing sequence likelihood", "TERMANAL::calcSeqFreq");
  Structure* S = res->getStructure();
  F->setQuery(*S);
  fasstSolutionSet matches = F->search();
  vector<Sequence> matchSeqs = F->getMatchSequences(matches);
  if (matchSeqs.empty()) MstUtils::error("No matches returned for residue " + MstUtils::toString(*res) + ", which is a problem for computing sequence likelihood", "TERMANAL::calcSeqFreq");
  mstreal seqFreq = 0.0;
  for (int i = 0; i < matchSeqs.size(); i++) if (matchSeqs[i][idx] == aa) seqFreq += 1.0;
  seqFreq /= matchSeqs.size();
  closestMatches[S] = matches;
  return seqFreq;
}

fasstSolutionSet TERMANAL::findTopN(const Structure& S, int N) {
  F->options().setRMSDCutoff(RMSDCalculator::rmsdCutoff(S, 1.1, 15));
  F->options().setMaxNumMatches(N);
  F->setQuery(S);
  fasstSolutionSet matches;
  while (true) {
    matches = F->search();
    if (matches.size() >= N) break;
    F->options().setRMSDCutoff(1.1*F->options().getRMSDCutoff());
  }
  return matches;
}
