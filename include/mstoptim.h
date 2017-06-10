#ifndef _MSTOPTIM_H
#define _MSTOPTIM_H

#include <vector>
#include <algorithm>    // for sort to work on linux
#include "mstlinalg.h"

using namespace std;

namespace MST {

class optimizerEvaluator {
  public:
    virtual vector<double> guessPoint() { return vector<double>(1, 0); }
    virtual double eval(const vector<double>& point) { return 0.0; }
    virtual double eval(const vector<double>& point, Vector& grad, vector<double> eps = vector<double>(0)) { grad = finiteDifferenceGradient(point, eps); return eval(point); }
    virtual Vector finiteDifferenceGradient(const vector<double>& point, vector<double> eps = vector<double>(0));
};

class Optim {

  public:

    // --- Nelder-Mead simplex algorithm as in Matlab's fminsearch
    static double fminsearch(optimizerEvaluator& E, int numIters, vector<double>& solution, bool verbose = false);
    static double gradDescent(optimizerEvaluator& E, vector<double>& solution, int numIters = 1000, double tol = 10E-6, bool verbose = false);

};

}
#endif
