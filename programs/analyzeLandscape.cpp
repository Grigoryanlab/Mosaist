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
  public:
    Landscape(mstreal _lowE = -20, mstreal _dE = 1, int _N = 40) {
      lowE = _lowE; dE = _dE; N = _N;
      seqsByEnergy.resize(N);
    }
    int getNumLevels() const { return N; }
    int getNumSeqsInLevel(int i) const { return seqsByEnergy[i].size(); }
    map<vector<int>, mstreal> getSeqsInLevel(int i) const { return seqsByEnergy[i]; }
    vector<map<vector<int>, mstreal> > getSeqsInAllLevels() const { return seqsByEnergy; }

    void addSequence(const vector<int>& seq, mstreal ener) {
      int levelIdx = (ener - lowE)/dE;
      if (levelIdx < 0) {
        cerr << "\tin Landscape::addSequence -> lowest energy is " << lowE << ", but trying to add a sequence with energy " << ener << ". Will ignore..." << endl;
        return;
      }
      if (levelIdx >= N) return;
      seqsByEnergy[levelIdx][seq] = ener;
    }

    vector<int> numberOfSeqsByEnergy() const {
      vector<int> nums(seqsByEnergy.size());
      for (int i = 0; i < seqsByEnergy.size(); i++) nums[i] = seqsByEnergy[i].size();
      return nums;
    }

    static void recordSolution(void* land, const vector<int>& sol, mstreal ener) {
      ((Landscape*) land)->addSequence(sol, ener);
    }

  private:
    mstreal lowE, dE;
    int N;
    vector<map<vector<int>, mstreal > > seqsByEnergy;
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

  mstreal dE = op.getReal("s", 1.0);
  int N = op.getInt("n", 40);
  MstUtils::assert(dE > 0, "energy step must be positive!");
  MstUtils::assert(N > 0, "number of energy intervals must be positive integer!");
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
  Landscape L(lowE, dE, N);

  // then run a series of MCs to discover sequences at various energy levels
  Etab.mc(10, 1000000, 1.0, 0.01, 1, &L, &Landscape::recordSolution);
  cout << "number of seqs by energy:\n" << MstUtils::vecToString(L.numberOfSeqsByEnergy(), "\n") << endl;
  // vector<vector<int> > seqs = MstUtils::keys(L.getSeqsInLevel(0));
  // for (int i = 0; i < seqs.size(); i++) {
  //   cout << Etab.solutionToSequence(seqs[i]).toString() << endl;
  // }

  // sub-sample to get a smaller number of sequences per energy region
  map<vector<int>, mstreal> subLand;
  if (op.isGiven("nat")) {
    // add the native
    subLand[Etab.sequenceToSolution(natSeq)] = Etab.scoreSequence(natSeq);
  }
  int n = op.getInt("k", 1000);
  for (int i = 0; i < L.getNumLevels(); i++) {
    map<vector<int>, mstreal> seqsMap = L.getSeqsInLevel(i);
    vector<vector<int> > seqs = MstUtils::keys(seqsMap);
    vector<int> inds = MstUtils::range(0, (int) seqsMap.size());
    MstUtils::shuffle(inds);
    for (int j = 0; j < MstUtils::min((int) seqs.size(), n); j++) {
      subLand[seqs[j]] = seqsMap[seqs[j]];
    }
  }

  // output distances for Matlab analysis
  vector<vector<int> > allSeqs = MstUtils::keys(subLand);
  cout << "calculating all-by-all for a subset of " << allSeqs.size() << " sequences..." << endl;
  fstream of; MstUtils::openFile(of, op.getString("o") + "_mat.dat", ios::out);
  for (int i = 0; i < allSeqs.size(); i++) {
    for (int j = 0; j < allSeqs.size(); j++) {
      of << distance(allSeqs[i], allSeqs[j]) << " ";
    }
    of << endl;
  }
  of.close();
  MstUtils::openFile(of, op.getString("o") + "_ener.dat", ios::out);
  for (int i = 0; i < allSeqs.size(); i++) of << subLand[allSeqs[i]] << endl;
  of.close();
}
