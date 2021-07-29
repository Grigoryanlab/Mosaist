#include "msttypes.h"
#include "dtermen.h"
#include "mstoptions.h"
#include "mstsequence.h"
#include <chrono>

mstreal sequenceComplexityPenalty(void* extra, const vector<int>& seq, EnergyTable& tab, int mutSite, int mutAA) {
  mstreal s = -(*((mstreal*) extra));
  return s * SeqTools::complexity(seq, mutSite, mutAA);
}

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Compute various things using an per-built energy table. Options:");
  op.addOption("e", "energy table file.", true);
  op.addOption("p", "PDB file to score the sequence of.");
  op.addOption("s", "single-letter amino-acid sequence to score.");
  op.addOption("opt", "optimize energy using MC sampling with default parameters. If an integer is specified, will use this many iterations per cycle (1E6 by default).");
  op.addOption("kTi", "if --opt is given, this will set the initial sampling temperature (default is 1.0).");
  op.addOption("kTf", "if --opt is given, this will set the final annealed temperature (default is 0.1).");
  op.addOption("lc", "add a low-complexity penalty to the energy scaled by this factor (should be positive).");
  op.addOption("cyc", "if --opt is given, this will set the number of MC cycles to run (default is 100).");
  op.addOption("o", "output file name of the energy table in case it needs to be written.");
  op.setOptions(argc, argv);
  if (op.isGiven("lc") && !op.isReal("lc")) {
    cout << op.usage() << endl;
    MstUtils::error("if given, --lc must be real");
  }

  EnergyTable E;
  E.readFromFile(op.getString("e"));

  // read sequences from various sources
  vector<Sequence> seqs;
  if (op.isGiven("p")) seqs.push_back(Sequence(Structure(op.getString("p"))));
  if (op.isGiven("s")) seqs.push_back(Sequence(op.getString("s"), "seq", op.getString("sd", "")));

  // score sequences
  for (int i = 0; i < seqs.size(); i++) {
    mstreal score = E.scoreSequence(seqs[i]);
    Sequence& seq = seqs[i];
    cout << score << " " << seq.toString() << endl;
  }

  // print mean and standard deviation
  if (op.isGiven("m")) cout << "mean " << E.meanEnergy() << endl;
  if (op.isGiven("std")) cout << "stdev " << E.energyStdEst(op.isInt("std") ? op.getInt("std") : 1000) << endl;

  if (op.isGiven("opt")) {
    int Ni = 1000000;
    if (op.isInt("opt")) Ni = op.getInt("opt");
    vector<int> lastSol;
    auto recordLast = [](void* cont, const vector<int>& sol, mstreal ener) { *((vector<int>*) cont) = sol; };
    vector<int> bestSol;
    if (op.isGiven("lc")) {
      mstreal lcsf = op.getReal("lc");
      bestSol = E.mc(op.getInt("cyc", 100), Ni, op.getReal("kTi", 10.0), op.getReal("kTf", 0.1), 1, &lastSol, recordLast, -1, &lcsf, &sequenceComplexityPenalty);
    } else {
      bestSol = E.mc(op.getInt("cyc", 100), Ni, op.getReal("kTi", 10.0), op.getReal("kTf", 0.1), 1, &lastSol, recordLast);
    }
    if (!lastSol.empty()) cout << "last sequence visited: " << (E.solutionToSequence(lastSol)).toString() << endl;
    cout << "lowest-energy sequence found: " << (E.solutionToSequence(bestSol)).toString() << endl;
    cout << "lowest-energy score found: " << E.scoreSolution(bestSol) << endl;
  }
}
