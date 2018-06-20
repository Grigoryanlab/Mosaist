#include "mstoptim.h"

using namespace MST;

Vector optimizerEvaluator::finiteDifferenceGradient(const vector<mstreal>& point, vector<mstreal> eps) {
  Vector grad(point.size(), 0);
  vector<mstreal> p = point;
  if (eps.empty()) eps = vector<mstreal>(point.size(), 10E-6);
  mstreal a, b;
  for (int i = 0; i < grad.length(); i++) {
    p[i] = point[i] - eps[i];
    a = eval(p);
    p[i] = point[i] + eps[i];
    b = eval(p);
    p[i] = point[i];
    grad[i] = (b - a)/(2*eps[i]);
  }
  return grad;
}

mstreal Optim::fminsearch(optimizerEvaluator& E, int numIters, vector<mstreal>& solution, bool verbose) {
  mstreal tol = 10E-6;
  Matrix x0(E.guessPoint()); // row vector
  int n = x0.length();
  Matrix simplex(n + 1, n);
  Matrix values(n + 1, 1);
  vector<int> order(n + 1);
  for (int i = 0; i < order.size(); i++) order[i] = i;

  // construct the initial simples
  for (int i = 0; i < simplex.size(); i++) {
    // the i-th point in the simplex is obtained by adjusting the i-th coordinate of the guess point
    simplex.row(i) = x0;
    if (i < n) {
      if (simplex(i, i) == 0) {
        simplex(i, i) = 0.00025; // default value
      } else {
        simplex(i, i) = 1.1 * simplex(i, i); // add 10%
      }
    }
    values(i) = E.eval(simplex.row(i));
  }

  // iterate until convergence of limit of iterations
  for (int k = 0; k < numIters; k++) {
    // sort points in the current simplex by function value
    sort(order.begin(), order.end(), [&values] (int i1, int i2) {return values(i1) < values(i2);});
    int b = order[0];   // best simplex
    int w = order[n-1]; // worst simplex to be kept
    int o = order[n];   // the very worst simplex (will be replaced)

    if (verbose) {
      cout << "iteration " << k+1 << ", the best simplex has value " << values(b) << " ";
    }

    // center of the current simplex
    Matrix m = simplex.mean(1);

    // quit if the current simplex is too small
    bool tooSmall = true;
    for (int i = 0; i < simplex.size(); i++) {
      if ((m - simplex.row(i)).norm() > tol) { tooSmall = false; break; }
    }
    if (tooSmall) break;

    // calculate the reflected point
    Matrix r = 2.0*m - simplex.row(o);
    mstreal fr = E.eval(r);
    if ((values(b) <= fr) && (fr < values(w))) {
      // accept the reflection
      simplex.row(o) = r;
      values(o) = fr;
      if (verbose) cout << "Reflect\n";
      continue;
    } else if (fr < values(b)) {
      // calculate the expansion point
      Matrix s = 2*(m - simplex.row(o)) + m;
      mstreal fs = E.eval(s);
      if (fs < fr) {
        // accept the expansion point
        simplex.row(o) = s;
        values(o) = fs;
        if (verbose) cout << "Expand\n";
        continue;
      } else {
        // accept the reflection
        simplex.row(o) = r;
        values(o) = fr;
        if (verbose) cout << "Reflect\n";
        continue;
      }
    } else if (fr >= values(w)) {
      // perform a contraction
      if (fr < values(o)) {
        Matrix c = 0.5*(r - m) + m;
        mstreal fc = E.eval(c);
        if (fc < fr) {
          // accept contraction point
          simplex.row(o) = c;
          values(o) = fc;
          if (verbose) cout << "Contract outside\n";
          continue;
        }
      } else {
        Matrix cc = 0.5*(simplex.row(o) - m) + m;
        mstreal fcc = E.eval(cc);
        if (fcc < values(o)) {
          // accept contraction point
          simplex.row(o) = cc;
          values(o) = fcc;
          if (verbose) cout << "Contract inside\n";
          continue;
        }
      }
    }
    // if nothing else worked, shrink
    for (int j = 1; j < order.size(); j++) {
      int i = order[j];
      simplex.row(i) = 0.5*(simplex.row(i) - simplex.row(b)) + simplex.row(b);
      values(i) = E.eval(simplex.row(i));
    }
    if (verbose) cout << "Shrink\n";
  }

  // find best solution
  int bestPoint;
  mstreal bestScore = MstUtils::min((vector<mstreal>) values, -1, -1, &bestPoint);
  solution = simplex.row(bestPoint);
  return bestScore;
}

mstreal Optim::gradDescent(optimizerEvaluator& E, vector<mstreal>& solution, int numIters, mstreal tol, bool verbose) {
  Vector x0(E.guessPoint());
  Vector x1(x0);
  mstreal gamma = 0.001;      // initial step size (gradient descent rate)
  mstreal v0, v1;
  Vector g0(x0.length()), g1(x0.length());

  v0 = E.eval(x0, g0);
  for (int i = 0; i < numIters; i++) {
    if (verbose) printf("gradient descent %d: %e (gamma = %e)\n", i+1, v0, gamma);
    x1 = x0 - gamma * g0;
    v1 = E.eval(x1, g1);
    if (fabs(v0 - v1) < tol) break;
    // adaptively change the learning rate
    if (v1 < v0) {
      gamma *= 1.1;
      x0 = x1; v0 = v1; g0 = g1;
    } else {
      gamma /= 1.2;
      if (MstUtils::closeEnough(gamma, 0.0)) break;
    }
  }
  solution = x0;
  return v0;
}

mstreal Optim::conjGradMin(optimizerEvaluator& E, vector<mstreal>& solution, int numIters, mstreal tol, bool verbose) {
  vector<mstreal> curr = E.guessPoint();
  vector<mstreal> next = curr;
  Vector h(curr.size()), g(curr.size()), g1(curr.size());
  mstreal v = E.eval(curr, g); g = -g; g1 = g; h = g;

  for (int i = 0; i < numIters; i++) {
    mstreal vNext = Optim::lineSearch(E, curr, next, h);
    if (v - vNext < tol) break;
    curr = next;

    // update direction
    v = E.eval(curr, g1); g1 = -g1;
    mstreal gamma = g1.norm2()/g.norm2(); // Fletcher–Reeves
    // mstreal gamma = g1.dot(g1 - g)/g.norm2(); // Polak–Ribière
    h = g1 + gamma*h;
    g = g1;
  }

  solution = curr;
  return v;
}

mstreal Optim::lineSearch(optimizerEvaluator& E, const vector<mstreal>& point, vector<mstreal>& solution, const Vector& _dir, mstreal startStepSize, bool verbose) {
  Vector x(point);
  Vector midGrad(x.length());
  Vector dir = _dir;
  Vector curr = x, low = x, mid = x, high = x;
  mstreal v;
  if (dir.length() == 0) {
    v = E.eval(x, dir); // set the line-search direction as the gradient
  } else {
    v = E.eval(x);
  }
  if (startStepSize <= 0) {
    // startStepSize = MstUtils::min(0.1/(dir.abs()).max(), 1.0);
    startStepSize = 0.1/(dir.abs()).max();
  }
  mstreal gamma = startStepSize;

  // bracket the minimum between low and high
  while (true) {
    // make step
    curr = x + gamma*dir;
    if (E.eval(curr) < v) {
      low = mid;
      mid = curr;
      gamma += startStepSize;
    } else {
      high = curr;
      break;
    }
  }

  // seaerch in the bracketed range
  for (int i = 0; i < 10; i++) {
    v = E.eval(mid, midGrad);
    if (midGrad.dot(dir) < 0) {
      // sub-divide the right side
      curr = (mid + high)/2;
      if (E.eval(curr) < v) {
        low = mid;
        mid = curr;
      } else {
        high = curr;
      }
    } else {
      // sub-divide the left side
      curr = (low + mid)/2;
      if (E.eval(curr) < v) {
        high = mid;
        mid = curr;
      } else {
        low = curr;
      }
    }
  }

  solution = mid;
  return E.eval(mid);
}
