#ifndef _MSTOPTIM_H
#define _MSTOPTIM_H

#include <vector>
#include <algorithm>    // for sort to work on linux
#include "mstlinalg.h"

using namespace std;

namespace MST {

class optimizerEvaluator {
  public:
    virtual vector<mstreal> guessPoint() { return vector<mstreal>(1, 0); }
    virtual mstreal eval(const vector<mstreal>& point) { return 0.0; }
    virtual mstreal eval(const vector<mstreal>& point, Vector& grad) { grad = finiteDifferenceGradient(point); return eval(point); }
    virtual Vector finiteDifferenceGradient(const vector<mstreal>& point, vector<mstreal> eps = vector<mstreal>(0));
};

class Optim {

  public:

    // --- Nelder-Mead simplex algorithm as in Matlab's fminsearch
    static mstreal fminsearch(optimizerEvaluator& E, int numIters, vector<mstreal>& solution, bool verbose = false);
    static mstreal gradDescent(optimizerEvaluator& E, vector<mstreal>& solution, int numIters = 1000, mstreal tol = 10E-8, bool verbose = false);

};

}
#endif
