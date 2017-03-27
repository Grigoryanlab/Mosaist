#ifndef _MSTOPTIM_H
#define _MSTOPTIM_H

#include <vector>
#include "mstlinalg.h"

using namespace std;

namespace MST {

class optimizerEvaluator {
  public:
    virtual vector<double> guessPoint() { return vector<double>(1, 0); }
    virtual double eval(vector<double>& point) { return 0.0; }
};

class Optim {

  public:

    // --- Nelder-Mead simplex algorithm as in Matlab's fminsearch
    static double fminsearch(optimizerEvaluator& E, int numIters, vector<double>& solution, bool verbose = false);

};

}
#endif
