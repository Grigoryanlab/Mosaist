#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <typeinfo>
#include <stdio.h>
#include "msttypes.h"
#include "mstoptions.h"
#include "mstsequence.h"

using namespace std;
using namespace MST;

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Tests MST::Sequence and MST::SeqTools. Options:");
  op.addOption("r", "do a rSearch test.");
  op.addOption("n", "number of sequences to generate (required if --c is given).");
  op.addOption("L", "sequence length (required if --c is given).");
  op.addOption("id", "sequence identity cutoff (default is 0.5).");
  op.addOption("b", "how much to bias the sequence selection (>= 1; 1 is default and means unbiased).");
  op.setOptions(argc, argv);
  srand(time(NULL) + (int) getpid());
  chrono::high_resolution_clock::time_point begin, end;

  // uSearch test
  if (op.isGiven("r")) {
    if (!op.isGiven("n") || !op.isInt("n")) MstUtils::error("--n must be given and must be integer");
    if (!op.isGiven("L") || !op.isInt("L")) MstUtils::error("--n must be given and must be integer");
    int N = op.getInt("n");
    int L = op.getInt("L");
    mstreal idCut = op.getReal("id", 0.5);
    mstreal b = op.getReal("b", 1.5);

    // build some biased random amino-acid sequences
    vector<Sequence> seqs(N, Sequence(vector<res_t>(L)));
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < L; j++) {
        mstreal r = pow(MstUtils::randUnit(), b);
        seqs[i][j] = (int) ceil(r*20) - 1;
      }
    }

    // apply rSearch
    cout << "doing rSearch..." << endl;
    begin = chrono::high_resolution_clock::now();
    vector<vector<int> > result = SeqTools::rSearch(seqs, idCut);
    end = chrono::high_resolution_clock::now();
    for (int i = 0; i < N; i++) {
      if (result[i].empty()) continue;
      cout << i << ": " << seqs[i].toString() << endl;
      for (int k = 0; k < result[i].size(); k++) {
        cout << "\t" << seqs[result[i][k]].toString() << endl;
      }
    }
    int searchTime = chrono::duration_cast<std::chrono::seconds>(end-begin).count();
    cout << "\trSearch took " << searchTime << " s" << std::endl;

    // a brute-force search
    cout << "doing a brute-force search..." << endl;
    begin = chrono::high_resolution_clock::now();
    int m = 0, nB = 0, nU;
    int nID = (int) ceil(1.0*L*idCut); // number of identities necessary to pass the cutoff
    for (int i = 0; i < N; i++) {
      vector<int> resultBrute;
      for (int j = 0; j < N; j++) {
        if (i == j) continue;
        int n = 0;
        for (int k = 0; k < L; k++) {
          if (seqs[i][k] == seqs[j][k]) n++;
        }
        if (n >= nID) resultBrute.push_back(j);
      }

      // compare
      nB += resultBrute.size();
      nU += result[i].size();
      if ((result[i].size() != resultBrute.size()) || (!MstUtils::setdiff(result[i], resultBrute).empty())) {
        cout << "rSearch for " << i << "    : " << MstUtils::vecToString(result[i]) << endl;
        cout << "brute-force for " << i << ": " << MstUtils::vecToString(resultBrute) << endl << endl;
        m++;
      }
    }
    end = chrono::high_resolution_clock::now();
    searchTime = chrono::duration_cast<std::chrono::seconds>(end-begin).count();
    cout << "\brute-force took " << searchTime << " s" << std::endl;

    cout << "rSearch found " << nU << " neighbors" << endl;
    cout << "brute-force found " << nB << " neighbors" << endl;
    cout << "DONE with " << m << " mistakes" << endl;
  } else {
    cout << op.usage() << endl;
  }
}
