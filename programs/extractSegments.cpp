#include "msttypes.h"
#include "mstsequence.h"
#include "mstoptions.h"
#include "dtermen.h"
#include "mstlinalg.h"

using namespace MST;

mstreal scoreAlignmentSegment(const vector<Sequence>& M, const vector<int>& positions, bool verbose = false) {
  // TODO: if we sort M based on positions, it will be already almost sorted from before, so will be linear
  // (so need to keep one copy of M for each segment and keep using it only for the sake of storing the order)
  // compute the frequency of each variant in the segment
  int N = M.size();
  map<Sequence, mstreal> f;
  for (int i = 0; i < M.size(); i++) f[M[i].subSequence(positions)] += 1.0/N;

  mstreal score = 0;
  int Neff;
  if (false) {
    // score is the sum of squares of frequencies (maximized when only one variant exists)
    for (const auto& it : f) score += (it.second) * (it.second);
    Neff = (int) round(1/score);
  } else {
    // score is Shannon information, with small sample correction (maximized when only one variant exists)
    int P = f.size(); // number of categories
    mstreal H = 0;    // entropy
    for (const auto& it : f) {
      H -= it.second * log(it.second);
    }
    H += P/(2*N);        // small sample correction (doi:10.1088/1751-8113/41/20/202001)
    // score = log(P) - H;  // information (log(P) is maximal entropy)
    score = -H;             // negative entropy
    Neff = (int) round(exp(H));
  }

  if (verbose) {
    cout << "segment score information:\n";
    cout << "\tscore = " << score << ", " << Neff << " effective sequences, " << f.size() << " unique sequences, " << M.size() << " sequences. Top few are:" << endl;
    vector<Sequence> uniqSeqs = MstUtils::keys(f);
    sort(uniqSeqs.begin(), uniqSeqs.end(), [f](const Sequence& i, const Sequence& j) { return f.at(i) > f.at(j); });
    for (int i = 0; i < MstUtils::min(10, (int) uniqSeqs.size()); i++) {
      cout << "\t\t" << uniqSeqs[i].toString() << " -> " << f[uniqSeqs[i]] << " (" << f[uniqSeqs[i]] * M.size() << ")" << endl;
    }
  }
  return score;
}

mstreal scoreAlignmentSegments(const vector<Sequence>& M, const vector<vector<int>>& segments) {
  mstreal score = 1;
  for (int i = 0; i < segments.size(); i++) score = MstUtils::min(score, scoreAlignmentSegment(M, segments[i]));
  return score;
}

vector<vector<int>> optimalSplit(const vector<Sequence>& M, const vector<int>& positions, bool verbose = false) {
  if (positions.size() < 2) MstUtils::error("have to have at least two positions to split!");
  vector<vector<int>> bestSegments; mstreal bestScore;
  for (int c = 0; c < 10; c++) {
    // start with a random split
    vector<vector<int>> segments(2);
    if (MstUtils::randUnit() < 0.5) { segments[0].push_back(positions[0]); segments[1].push_back(positions[1]); }
    else { segments[0].push_back(positions[1]); segments[1].push_back(positions[0]); }
    for (int i = 2; i < positions.size(); i++) {
      if (MstUtils::randUnit() < 0.5) segments[0].push_back(positions[i]);
      else segments[1].push_back(positions[i]);
    }

    // now do MC sampling that randomly moves positions between segments and scores
    mstreal score = scoreAlignmentSegments(M, segments);
    if (verbose) {
      cout << "------------" << endl;
      cout << "segments[0] = " << MstUtils::vecToString(segments[0]) << endl;
      cout << "segments[1] = " << MstUtils::vecToString(segments[1]) << endl;
      cout << "initial score is " << score << endl;
    }
    for (int i = 0; i < 1E4; i++) {
      // pick a random position and re-assign
      int fromSide, toSide;
      if (MstUtils::randUnit() < 0.5) { fromSide = 0; toSide = 1; }
      else { fromSide = 1; toSide = 0; }
      int idx = MstUtils::randInt(segments[fromSide].size());
      segments[toSide].push_back(segments[fromSide][idx]);
      segments[fromSide].erase(segments[fromSide].begin() + idx);

      // score
      mstreal newScore = scoreAlignmentSegments(M, segments);

      // accept or reject
      if (newScore > score) {
        score = newScore;
        if (verbose) {
          cout << "score became " << score << ", with a " << segments[0].size() << " - " << segments[1].size() << " split" << endl;
          cout << "segments[0] = " << MstUtils::vecToString(MstUtils::keys(MstUtils::contents(segments[0]))) << endl;
          cout << "segments[1] = " << MstUtils::vecToString(MstUtils::keys(MstUtils::contents(segments[1]))) << endl;
        }
      } else {
        segments[fromSide].push_back(segments[toSide].back());
        segments[toSide].pop_back();
      }
    }
    if ((c == 0) || (score > bestScore)) {
      bestScore = score;
      bestSegments = segments;
    }
  }

  return bestSegments;
}

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Given an MSA, splits it into sub-alignments, such that the total diversity of variants within all sub-alignment is minimized. After such a split, the MSA can essentially be well approximated as a convolution of multiple MSA, each with as few variants as possible. Options:");
  op.addOption("f", "FASTA file with the MSA.");
  op.addOption("s", "a flat MSA file (one single-letter sequence per line).");
  op.addOption("e", "energy table. If specified, sequences will be generated by performing MC on this table.");
  op.addOption("n", "native sequence. If given, will treat sequences as perturbation of the native.");
  op.addOption("kT", "kT for the MC sampling.");
  op.addOption("nc", "number of MC cycles. One sequence will be taken from each cycle (the one sampled last).");
  op.addOption("ni", "number of iteration per cycle");
  op.addOption("kTf", "if specified, will anneal from the temperature specified by --kT to this temperature during sampling.");
  op.setOptions(argc, argv);
  if (!op.isGiven("f") && !op.isGiven("s") && !op.isGiven("e")) { cout << op.usage(); MstUtils::error("neither --f, --s, nor --e were given"); }
  vector<Sequence> M;
  if (op.isGiven("f")) SeqTools::readFasta(op.getString("f"), M);
  if (op.isGiven("s")) SeqTools::readSequences(op.getString("s"), M);

  long int x = std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
  srand(x);

  if (op.isGiven("e")) {
    MstUtils::assertCond(op.isGiven("nc") && op.isGiven("ni") && op.isGiven("kT"), "not all MC-related parameters are specifieed");
    EnergyTable Etab(op.getString("e"));
    for (int c = 0; c < op.getInt("nc"); c++) {
      vector<int> bestSol = Etab.mc(1, op.getInt("ni"), op.getReal("kT"), op.getReal("kTf", op.getReal("kT")));
      M.push_back(Etab.solutionToSequence(bestSol));
      cout << M.back().toString() << endl;
    }
  }
  if (M.empty()) MstUtils::error("no sequences were read/generated!");

  // an approach based on encoding each sequence as a delta from the native sequence
  if (op.isGiven("n")) {
    Sequence N(op.getString("n"));
    MstUtils::assertCond(N.length() == M[0].length(), "native sequence of different length than sequences in the MSA");
    int L = N.length();

    // build the alphabet at each position from the MSA and native
    int W = 0; // total number of non-natives allowed at all positions
    vector<map<res_t, int>> alpha(L);
    for (int i = 0; i < L; i++) {
      alpha[i][N[i]] = alpha[i].size();
      set<res_t> inMSA;
      for (int k = 0; k < M.size(); k++) {
        if (inMSA.find(M[k][i]) == inMSA.end()) {
          inMSA.insert(M[k][i]);
          alpha[i][M[k][i]] = alpha[i].size();
          W++;
        }
      }
    }

    // encode each matrix as a delta
    Matrix D(M.size(), W, 0);
    int off = 0;
    for (int i = 0; i < L; i++) {
      for (int k = 0; k < M.size(); k++) {
        if (alpha[i][M[k][i]]) { // if it is native, skip
          D(k, off + alpha[i][M[k][i]] - 1) = 1;
        }
      }
      off += alpha[i].size() - 1;
    }

    // dump each row
    for (int i = 0; i < D.numRows(); i++) {
      for (int j = 0; j < D.numCols(); j++) {
        cout << (int) D(i, j) << " ";
      }
      cout << endl;
    }
  } else {
    // an approach based on random splitting and minimizing the effective number
    // of sequences (essentially maximizing information in each segment or mini-
    // mizing the entropy)

    // -- optimally split the MSA into two segments
    cout << "before splitting:" << endl;
    vector<vector<int>> segments(1, MstUtils::range(0, (int) M[0].size()));
    scoreAlignmentSegment(M, segments[0], true);
    for (int c = 0; c < 2; c++) {
      cout << "cycle " << c << endl;
      vector<vector<int>> newSegments;
      for (int i = 0; i < segments.size(); i++) {
        cout << "\tsplitting segment " << i << ": " << MstUtils::vecToString(MstUtils::keys(MstUtils::contents(segments[i]))) << endl;
        vector<vector<int>> split = optimalSplit(M, segments[i]);
        cout << "\n\t\tsub-segment 0: " << MstUtils::vecToString(MstUtils::keys(MstUtils::contents(split[0]))) << endl;
        scoreAlignmentSegment(M, split[0], true);
        cout << "\n\t\tsub-segment 1: " << MstUtils::vecToString(MstUtils::keys(MstUtils::contents(split[1]))) << endl;
        scoreAlignmentSegment(M, split[1], true);
        newSegments.push_back(split[0]);
        newSegments.push_back(split[1]);
      }
      segments = newSegments;
    }
    // cout << endl << "***" << endl;
    // cout << "segment 0: " << MstUtils::vecToString(MstUtils::keys(MstUtils::contents(segments[0]))) << endl;
    // mstreal s1 = scoreAlignmentSegment(M, segments[0], true);
    // cout << "segment 1: " << MstUtils::vecToString(MstUtils::keys(MstUtils::contents(segments[1]))) << endl;
    // scoreAlignmentSegment(M, segments[1], true);
    // cout << endl << "***" << endl;
    //
    // cout << "before splitting:" << endl;
    // scoreAlignmentSegment(M, MstUtils::range(0, (int) M[0].size()), true);
    // vector<vector<int>> segments = optimalSplit(M, MstUtils::range(0, (int) M[0].size()));
    // cout << endl << "***" << endl;
    // cout << "segment 0: " << MstUtils::vecToString(MstUtils::keys(MstUtils::contents(segments[0]))) << endl;
    // mstreal s1 = scoreAlignmentSegment(M, segments[0], true);
    // cout << "segment 1: " << MstUtils::vecToString(MstUtils::keys(MstUtils::contents(segments[1]))) << endl;
    // scoreAlignmentSegment(M, segments[1], true);
    // cout << endl << "***" << endl;

    // vector<vector<int>> segmentsA = optimalSplit(M, segments[0]);
    // vector<vector<int>> segmentsB = optimalSplit(M, segments[1]);
    // cout << endl << "***" << endl;
    // cout << "segment 0: " << MstUtils::vecToString(MstUtils::keys(MstUtils::contents(segmentsA[0]))) << endl;
    // cout << "segment 1: " << MstUtils::vecToString(MstUtils::keys(MstUtils::contents(segmentsA[1]))) << endl;
    // cout << "segment 3: " << MstUtils::vecToString(MstUtils::keys(MstUtils::contents(segmentsB[0]))) << endl;
    // cout << "segment 4: " << MstUtils::vecToString(MstUtils::keys(MstUtils::contents(segmentsB[1]))) << endl;
    // cout << endl << "***" << endl;
  }
}
