#include "msttermanal.h"

mstreal TERMANAL::structureScore(const Structure& S, const vector<Residue*>& central, bool verbose) {
  if (F == NULL) MstUtils::error("FASST object not set", "TERMANAL::structureScore");
  fasstSearchOptions origOpts = F->options(); // in case the same FASST object is being shared by others
  F->setOptions(fasstSearchOptions());
  F->options().setRedundancyCut(0.6);

  // find top hits
  int N = 50;
  F->options().setRMSDCutoff(RMSDCalculator::rmsdCutoff(S, 1.1, 15));
  F->setQuery(S);
  F->options().setMinNumMatches(N);
  F->options().setMaxNumMatches(N);
  fasstSolutionSet matches = F->search();
  vector<Sequence> matchSeqs = F->getMatchSequences(matches);

  // sequence likelihood = fraction of top hits that share the right amino acid
  vector<mstreal> seqLike(central.size(), 0);
  for (int i = 0; i < central.size(); i++) {
    Residue* res = central[i];
    int idx = res->getResidueIndex();
    if (idx >= S.residueSize()) MstUtils::error("index of central residue beyond capacity of TERM passed", "TERMANAL::structureScore");
    res_t aa = SeqTools::aaToIdx(res->getName());
    if (aa == SeqTools::unknownIdx()) MstUtils::error("uknown residue name in central residue " + MstUtils::toString(*res) + ", which is a problem for computing sequence likelihood", "TERMANAL::structureScore");
    for (int k = 0; k < matchSeqs.size(); k++) {
      if (matchSeqs[k][idx] == aa) seqLike[i] += 1.0;
    }
    seqLike[i] /= matchSeqs.size();
  }
  mstreal seqLikeComb = (!central.empty()) ? CartesianPoint(seqLike).mean() : 1.0;
  if (verbose) cout << "sequence likelihood = [" << MstUtils::vecToString(seqLike) << "], overall = " << seqLikeComb << endl;

  // structure frequency = number of matches with RMSD below that of the Nth
  // match to the top native representative of the query, divided by N
  F->setQuery(F->getMatchStructure(matches[0]));
  fasstSolutionSet matchesToClosestNative = F->search();
  mstreal r = matchesToClosestNative[N - 1].getRMSD();
  int n = 0;
  for (; n < matches.size(); n++) {
    if (matches[n].getRMSD() > r) break;
  }
  mstreal structFreq = (1.0*n)/N;
  if (verbose) cout << "structure frequency = " << structFreq << endl;

  // structure score combines the two
  mstreal eps = 0.01;
  mstreal ss = -log(structFreq * seqLikeComb + eps);

  // recover original search options
  F->setOptions(origOpts);

  return ss;
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
