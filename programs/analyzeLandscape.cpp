#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>
#include <typeinfo>
#include <stdio.h>

#include "msttypes.h"
#include "mstoptions.h"
#include "dtermen.h"

using namespace std;
using namespace MST;

class Landscape {
  typedef vector<int> sol_t;
  public:
    Landscape(mstreal _lowE = -20, mstreal _dE = 1, int _N = 40, int _maxPerBin = -1) {
      lowE = _lowE; dE = _dE; N = _N;
      seqsByEnergy.resize(N);
      visitedByEnergy.resize(N);
      maxPerBin = _maxPerBin;
      numHits.resize(N, 0);
      resetMeans();
    }
    void resetMeans() { runN = 0; meanEner = 0; }
    int getNumLevels() const { return N; }
    int getNumSeqsInLevel(int i) const { return seqsByEnergy[i].size(); }
    vector<pair<sol_t, mstreal> > getSeqsInLevel(int i) const { return seqsByEnergy[i]; }
    vector<vector<pair<sol_t, mstreal> > > getSeqsInAllLevels() const { return seqsByEnergy; }
    vector<int> getNumHits() const { return numHits; }
    int getTotalNeed(int cap) const {
      int need = 0;
      for (int i = 0; i < numHits.size(); i++) need += MstUtils::max(0, cap - numHits[i]);
      return need;
    }
    int getNumLevelsInNeed(int cap) const {
      int need = 0;
      for (int i = 0; i < numHits.size(); i++) {
        if ((cap - numHits[i]) > 0) need++;
      }
      return need;
    }
    mstreal getMeanDefitiency(int cap) const { // on average, where are we missing hits?
      if (MstUtils::max(numHits) < cap) {
        mstreal num = 0, den = 0;
        for (int i = 0; i < numHits.size(); i++) {
          num += MstUtils::max(0, cap - numHits[i])*(lowE + (i+ 0.5)*dE);
          den += MstUtils::max(0, cap - numHits[i]);
        }
        if (den == 0) return lowE + (N/2.0)*dE;
        return num/den;
      }

      /* if we are dealing with a multi-modal distribution of deficiencies, then
       * switch to finding the first peak of deficiency. That way, we won't move
       * on until we clear that peak. */
      int maxNeed = 0; int maxNeedLevel = numHits.size() - 1;
      for (int i = numHits.size() - 1; i >= 0; i--) {
        int need = MstUtils::max(0, cap - numHits[i]);
        if (need > maxNeed) { maxNeed = need; maxNeedLevel = i; }
        if (maxNeed > 0) break;
      }
      return lowE + (maxNeedLevel + 0.5)*dE;
    }
    mstreal getMeanEnergy() const { return meanEner; }
    int energyLevelIndex(mstreal ener) const { return (ener - lowE)/dE; }

    void addSequence(const sol_t& seq, mstreal ener) {
      meanEner = (meanEner*runN + ener)/(runN + 1); runN++;
      int levelIdx = energyLevelIndex(ener);
      if (levelIdx < 0) {
        cerr << "\tin Landscape::addSequence -> lowest energy is " << lowE << ", but trying to add a sequence with energy " << ener << ". Will ignore..." << endl;
        return;
      }
      if (levelIdx >= N) return;

      set<sol_t >& visitedInLevel = visitedByEnergy[levelIdx];
      if (visitedInLevel.find(seq) == visitedInLevel.end()) {
        numHits[levelIdx]++;
        vector<pair<sol_t, mstreal> >& seqsInLevel = seqsByEnergy[levelIdx];
        if ((maxPerBin > 0) && (seqsInLevel.size() >= maxPerBin)) {
          // erase a random sequence
          int ri = MstUtils::randInt(seqsInLevel.size());
          visitedInLevel.erase(seqsInLevel[ri].first);
          // replace it with the new sequence
          visitedInLevel.insert(seq);
          seqsInLevel[ri].first = seq;
          seqsInLevel[ri].second = ener;
        } else {
          // insert a new sequence
          visitedInLevel.insert(seq);
          seqsInLevel.push_back(pair<sol_t, mstreal>(seq, ener));
        }
      }
    }

    vector<int> numberOfSeqsByEnergy() const {
      vector<int> nums(seqsByEnergy.size());
      for (int i = 0; i < seqsByEnergy.size(); i++) nums[i] = seqsByEnergy[i].size();
      return nums;
    }

    static void recordSolution(void* land, const sol_t& sol, mstreal ener) {
      ((Landscape*) land)->addSequence(sol, ener);
    }

  private:
    mstreal lowE, dE;    // low energy limit and bin energy width
    int N;               // number of bins
    int maxPerBin;       // max number of counts to keep per bin. Will randomly
                         // replace an existing solution in the bin if this limit
                         // is reached and a new unique solution is added
    vector<vector<pair<sol_t, mstreal> > > seqsByEnergy;
    vector<set<sol_t> > visitedByEnergy;
    vector<int> numHits; // total number of unique solutions that was added to the bin

    // counters for computing mean properties (over all solutions, no matter whether unique)
    int runN;            // running total number of solutions added
    mstreal meanEner;    // the running mean energy
};

int distance(const vector<int>& seqi, const vector<int>& seqj) {
  int d = 0;
  for (int i = 0; i < seqi.size(); i++) {
    d += (seqi[i] != seqj[i]);
  }
  return d;
}

int main(int argc, char** argv) {
  MstOptions op;
  op.setTitle("Analyzes the sequence landscape encoded by a given sequence-level pseudo-energy table. Options:");
  op.addOption("e", "energy table file.", true);
  op.addOption("o", "output base for Matlab analysis.", true);
  op.addOption("s", "step in energy units for building \"contour lines\" in the energy landscape. Default is 1.0.");
  op.addOption("n", "number of intervals of this size to track. Default is 40.");
  op.addOption("k", "number of sequences to sub-sample in each interval at the end to compute an embedding. Default is 1000.");
  op.addOption("nat", "native sequence, single-letter.");
  op.setOptions(argc, argv);
  MstUtils::setSignalHandlers();

  mstreal dE = op.getReal("s", 1.0);
  int N = op.getInt("n", 40);
  MstUtils::assertCond(dE > 0, "energy step must be positive!");
  MstUtils::assertCond(N > 0, "number of energy intervals must be positive integer!");
  srand(time(NULL) + (int) getpid());

  EnergyTable Etab(op.getString("e"));

  // first, run a long-ish MC to try to get the best energy
  Sequence natSeq;
  if (op.isGiven("nat")) {
    natSeq = Sequence(op.getString("nat"));
    cout << "native sequence energy is " << Etab.scoreSequence(natSeq) << endl;
  }
  vector<int> bestSol = Etab.mc(100, 1000000, 1.0, 0.01);
  mstreal lowE = Etab.scoreSolution(bestSol);
  cout << "lowest energy found is " << lowE << endl;
  cout << "lowest-energy sequence: " << (Etab.solutionToSequence(bestSol)).toString() << endl;
  cout << "mean energy is " << Etab.meanEnergy() << endl;
  cout << "estimated energy standard deviation is " << Etab.energyStdEst() << endl;
  int Nsub = op.getInt("k", 1000); // number of sub-sampled sequences we are looking for in each energy bin
  Landscape L(lowE, dE, N, Nsub);

  // then run a series of MCs to discover sequences at various energy levels
  int Ni = 1000000, Ne = 1000000;
  mstreal kT = 1.0;
  int cap = 1000*Nsub; // we'd like to stop when the number of hits in each energy
                       // level reaches many times the number of sub-sampled sequences per bin
  int stuck = 0;       // for how many cycles have we been stuck? (no improvement)
  int prevNeed = L.getTotalNeed(cap); int need = prevNeed;
  while ((MstUtils::min(L.getNumHits()) < cap) || (stuck > 100)) {
    Etab.mc(1, Ni, kT, kT, 1, &L, &Landscape::recordSolution, Ne);
    mstreal Ed = L.getMeanDefitiency(cap);
    mstreal Es = L.getMeanEnergy();
    int need = L.getTotalNeed(cap);
    cout << "need = " << need << " over " << L.getNumLevelsInNeed(cap) << " level(s); deficiency at " << Ed << ", while sampling with kT = " << kT << " was at " << Es << endl;
    if (Ed < Es) kT /= 1.02;
    else kT *= 1.05;
    kT = MstUtils::max(kT, 1E-2);
    L.resetMeans();

    // are we improving?
    if (need == prevNeed) stuck++;
    else stuck = 0;
    prevNeed = need;
    if (stuck > 100) break;
  }

  // get all sequences sampled
  vector<vector<int> > allSeqs;
  vector<mstreal> allEnergies;
  for (int i = 0; i < L.getNumLevels(); i++) {
    vector<pair<vector<int>, mstreal> > seqsInLevel = L.getSeqsInLevel(i);
    for (int j = 0; j < seqsInLevel.size(); j++) {
      allSeqs.push_back(seqsInLevel[j].first);
      allEnergies.push_back(seqsInLevel[j].second);
    }
  }

  // add the native
  if (op.isGiven("nat")) {
    allSeqs.insert(allSeqs.begin(), Etab.sequenceToSolution(natSeq));
    allEnergies.insert(allEnergies.begin(), Etab.scoreSequence(natSeq));
  }

  // output distances for Matlab analysis
  cout << "calculating all-by-all for a subset of " << allSeqs.size() << " sequences..." << endl;
  fstream of; MstUtils::openFile(of, op.getString("o") + "_mat.dat", ios::out);
  for (int i = 0; i < allSeqs.size(); i++) {
    for (int j = 0; j < allSeqs.size(); j++) {
      of << distance(allSeqs[i], allSeqs[j]) << " ";
    }
    of << endl;
  }
  of.close();

  // output energies
  MstUtils::openFile(of, op.getString("o") + "_ener.dat", ios::out);
  for (int i = 0; i < allEnergies.size(); i++) of << allEnergies[i] << endl;
  of.close();

  // output sequences
  MstUtils::openFile(of, op.getString("o") + "_seq.dat", ios::out);
  for (int i = 0; i < allSeqs.size(); i++) of << (Etab.solutionToSequence(allSeqs[i])).toString() << endl;
  of.close();
}
