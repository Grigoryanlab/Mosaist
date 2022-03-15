#include "msttypes.h"
#include "dtermen.h"
#include "mstoptions.h"
#include "mstsequence.h"
#include <chrono>

mstreal sequenceComplexityPenalty(void* extra, const vector<int>& seq, EnergyTable& tab, int mutSite, int mutAA) {
  mstreal s = -(*((mstreal*) extra));
  return s * SeqTools::complexity(seq, mutSite, mutAA);
}

mstreal sequenceComplexityPenaltySimple(void* extra, const vector<int>& seq, EnergyTable& tab, int mutSite, int mutAA) {
  mstreal s = -((mstreal*) extra)[0];
  mstreal fcut = ((mstreal*) extra)[1];

  // get counts of all letters
  vector<int> counts(100, 0);
  int maxCount = 0;
  for (int aa : seq) {
    if (aa >= counts.size()) counts.resize(aa + 1, 0);
    counts[aa]++;
    maxCount = max(maxCount, counts[aa]);
  }

  mstreal f = (maxCount*1.0)/seq.size();
  if (f < fcut) return 0;

  return s * (f - fcut)*(f - fcut)/((1 - fcut)*(1 - fcut));
}

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Loads a pre-built energy table and scores sequences or performs MCMC optimization to design a new sequence. Options:");
  op.addOption("e", "Energy table file.", true);
  op.addOption("p", "PDB file. If provided, will score the sequence of the structure. Note: must have the same number of residues as the energy table.");
  op.addOption("s", "Single-letter amino-acid sequence. If provided, will score. Must have the same number of residues as the energy table");
  op.addOption("opt", "If provided, will perform MCMC simulated annealing to find the optimal sequence with default parameters. If an integer is specified, will use this many iterations per cycle (otherwise 1E6 by default).");
  op.addOption("kTi", "if --opt is given, this will set the initial sampling temperature (default is 1.0).");
  op.addOption("kTf", "if --opt is given, this will set the final annealed temperature (default is 0.1).");
  op.addOption("lc", "if --opt is givem, will add a low-complexity penalty to the energy scaled by this factor (should be positive).");
  op.addOption("fcut", "if --opt is given, will set a limit on the fraction of positions allowed to be occupied by a single amino acid type. If specified, will use the simpler complexity penalty rather than the one based on number of arrangements of the letter distribution.");
  op.addOption("cyc", "if --opt is given, this will set the number of MC cycles to run (default is 100).");
  op.addOption("randomSeed","If --randomSeed is given, will set a new random seed each time the program is run. Otherwise will use the same random seed and provide consistent results");
  op.setOptions(argc, argv);
  if (op.isGiven("lc") && !op.isReal("lc")) {
    cout << op.usage() << endl;
    MstUtils::error("if given, --lc must be real");
  }
  if (op.isGiven("randomSeed")) srand(time(NULL) + int(getpid()));
  else srand(1);

  srand(time(NULL) + int(getpid()));

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

      if (op.isGiven("fcut")) {
        mstreal params[] = {op.getReal("lc"), op.getReal("fcut")};
        bestSol = E.mc(op.getInt("cyc", 100), Ni, op.getReal("kTi", 10.0), op.getReal("kTf", 0.1), 1, &lastSol, recordLast, -1, &params, &sequenceComplexityPenaltySimple);
      } else {
        mstreal lcsf = op.getReal("lc");
        bestSol = E.mc(op.getInt("cyc", 100), Ni, op.getReal("kTi", 10.0), op.getReal("kTf", 0.1), 1, &lastSol, recordLast, -1, &lcsf, &sequenceComplexityPenalty);
      }
    } else {
      bestSol = E.mc(op.getInt("cyc", 100), Ni, op.getReal("kTi", 10.0), op.getReal("kTf", 0.1), 1, &lastSol, recordLast);
    }
    if (!lastSol.empty()) cout << "last sequence visited: " << (E.solutionToSequence(lastSol)).toString() << endl;
    cout << "lowest-energy sequence found: " << (E.solutionToSequence(bestSol)).toString() << endl;
    cout << "lowest-energy score found: " << E.scoreSolution(bestSol) << endl;
  }
}
