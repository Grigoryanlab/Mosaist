#include "mstoptim.h"

double Optim::fminsearch(optimizerEvaluator& E, int numIters, vector<double>& solution, bool verbose) {
  double tol = 10E-6;
  Vector x0 = E.guessPoint();
  int n = x0.size();
  Matrix simplex(n + 1, x0.length());
  Vector values(n + 1);
  vector<int> order(n + 1);
  for (int i = 0; i < order.size(); i++) order[i] = i;

  // construct the initial simples
  for (int i = 0; i < simplex.size(); i++) {
    // the i-th point in the simplex is obtained by adjusting the i-th coordinate of the guess point
    simplex[i] = x0;
    if (i < n) {
      if (simplex[i][i] == 0) {
        simplex[i][i] = 0.00025; // default value
      } else {
        simplex[i][i] = 1.05 * simplex[i][i]; // add 5%
      }
    }
    values[i] = E.eval(simplex[i]);
  }

  // iterate until convergence of limit of iterations
  for (int k = 0; k < numIters; k++) {
    // sort points in the current simplex by function value
    sort(order.begin(), order.end(), [&values] (int i1, int i2) {return values[i1] < values[i2];});
    int b = order[0];   // best simplex
    int w = order[n-1]; // worst simplex to be kept
    int o = order[n];   // the very worst simplex (will be replaced)

    if (verbose) {
      cout << "iteration " << k+1 << ", the current " << order.size() << " simplex points:" << endl;
      for (int i = 0; i < order.size(); i++) {
        cout << MstUtils::vecToString(simplex[order[i]]) << " -> " << values[order[i]] << endl;
      }
    }

    // center of the current simplex
    vector<double> m(n, 0);
    for (int i = 0; i < simplex.size(); i++) Util::increment(m, simplex[i]);
    Util::scale(m, 1.0/simplex.size());

    // quit if the current simplex is too small
    bool tooSmall = true;
    for (int i = 0; i < simplex.size(); i++) {
      vector<double> diff = Util::subtract(m, simplex[i]);
      if (Util::norm(diff) > tol) { tooSmall = false; break; }
    }
    if (tooSmall) break;

    // calculate the reflected point
    vector<double> r(m);
    Util::scale(r, 2.0);
    Util::decrement(r, simplex[o]);
    double fr = E.eval(r);
    if ((values[b] <= fr) && (fr < values[w])) {
      // accept the reflection
      simplex[o] = r;
      values[o] = fr;
      if (verbose) cout << "Reflect\n";
      continue;
    } else if (fr < values[b]) {
      // calculate the expansion point
      vector<double> s = Util::subtract(m, simplex[o]);
      Util::scale(s, 2.0);
      Util::increment(s, m);
      double fs = E.eval(s);
      if (fs < fr) {
        // accept the expansion point
        simplex[o] = s;
        values[o] = fs;
        if (verbose) cout << "Expand\n";
        continue;
      } else {
        // accept the reflection
        simplex[o] = r;
        values[o] = fr;
        if (verbose) cout << "Reflect\n";
        continue;
      }
    } else if (fr >= values[w]) {
      // perform a contraction
      if (fr < values[o]) {
        vector<double> c = Util::subtract(r, m);
        Util::scale(c, 0.5);
        Util::increment(c, m);
        double fc = E.eval(c);
        if (fc < fr) {
          // accept contraction point
          simplex[o] = c;
          values[o] = fc;
          if (verbose) cout << "Contract outside\n";
          continue;
        }
      } else {
        vector<double> cc = Util::subtract(simplex[o], m);
        Util::scale(cc, 0.5);
        Util::increment(cc, m);
        double fcc = E.eval(cc);
        if (fcc < values[o]) {
          // accept contraction point
          simplex[o] = cc;
          values[o] = fcc;
          if (verbose) cout << "Contract inside\n";
          continue;
        }
      }
    }
    // if nothing else worked, shrink
    for (int j = 1; j < order.size(); j++) {
      int i = order[j];
      vector<double> v = Util::subtract(simplex[i], simplex[b]);
      Util::scale(v, 0.5);
      Util::increment(v, simplex[b]);
      simplex[i] = v;
      values[i] = E.eval(simplex[i]);
    }
    if (verbose) cout << "Shrink\n";
  }

  // find best solution
  int bestPoint;
  double bestScore = Util::min(values, -1, -1, &bestPoint);
  solution = simplex[bestPoint];
  return bestScore;
}
