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
    mstreal vNext = Optim::lineSearch(E, curr, next, h, 0.01, verbose);
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

mstreal Optim::lineSearch(optimizerEvaluator& E, const vector<mstreal>& point, vector<mstreal>& solution, const Vector& dir, mstreal startStepSize, bool verbose) {
  Vector x0 = point, g0 = x0, x1 = x0, g1 = x0;
  mstreal v0 = E.eval(x0, g0), v1; g0 = g0.getUnit();
  mstreal gamma = startStepSize;      // initial step size
  mstreal tol = 10E-8;
  for (int i = 0; 100; i++) {
    if (verbose) printf("line search %d: %e (gamma = %e)\n", i+1, v0, gamma);
    x1 = x0 - gamma * dir.dot(g0) * dir;
    v1 = E.eval(x1, g1);
    if (fabs(v0 - v1) < tol) break;
    // adaptively change the learning rate
    if (v1 < v0) {
      gamma *= 1.5;
      x0 = x1; v0 = v1; g0 = g1.getUnit();
    } else {
      gamma /= 2;
      if (MstUtils::closeEnough(gamma, 0.0)) break;
    }
  }
  solution = x0;
  return v0;
}

vector<mstreal> Optim::langevinDynamics(optimizerEvaluator& E, const vector<mstreal>& masses, mstreal timeStep, mstreal gamma, mstreal kT, int numIters, vector<vector<mstreal> >& trajectory, int saveInterval, bool verbose) {
  // Integrator type:
  // 0 -- Brooks-Brunger-Karplus (BBK) integration scheme, appropriate for small-gamma regime
  // 1 -- Bussi and Parrinello (May, 2018), works well for large gamma (https://arxiv.org/pdf/0803.4083.pdf)
  // 2 -- a second-order Langevin integrator (see doc/molecular_dynamics_2015.pdf)
  int integrator = 1;

  if (saveInterval < 0) saveInterval = numIters;
  vector<mstreal> ener; // snapshot energies

  // some necessary vectors (needed by all integrators)
  Vector x0(E.guessPoint());    // starting positions (row vector)
  int df = x0.size();           // number of degrees of freedom
  int dim = df / masses.size(); // number of degrees of freedom per particle
  int np = df / dim;            // number of particles
  Vector V(df);                 // velocities
  Vector f(df);                 // forces divided by mass, at current point
  Vector f1(df);                // forces divided by mass, at next point
  Vector C(df);                 // means different things depending on the integrator, see below
  Vector x(x0);                 // current positions
  Vector m(df);                 // masses
  Vector R(df);                 // normally-distributed random numbers with zero mean and unit variance
  Vector R1(df);                // normally-distributed random numbers with zero mean and unit variance
  if (df % masses.size() != 0) {
    MstUtils::error(MstUtils::toString(df) + "-dimensional system, but received " + MstUtils::toString(masses.size()) + " masses", "Optim::langevinDynamics");
  }
  for (int i = 0; i < df; i += dim) {
    for (int k = 0; k < dim; k++) m[i + k] = masses[i / dim];
  }
  for (int i = 0; i < df; i++) V[i] = MstUtils::randNormal(0, sqrt(kT/m[i])); // assign random velocities
  mstreal timeStep2 = timeStep*timeStep;
  mstreal timeStep12 = sqrt(timeStep);
  mstreal timeStep32 = pow(timeStep12, 3);

  // arrays for Bussi and Parrinello integrator
  mstreal c1; Vector c2, sig;
  if (integrator == 1) {
    c1 = exp(-gamma*timeStep/2);  // equation 13a
    c2 = Vector(df);              // c2 from equation 13b divided by mass
    for (int i = 0; i < df; i++) c2[i] = sqrt((1 - c1*c1)*kT/m[i]);
  }
  if (integrator == 2) {
    sig = Vector(df);             // sqrt(2*kT*gamma/m)
    for (int i = 0; i < df; i++) sig[i] = sqrt(2*kT*gamma/m[i]);
  }

  mstreal en = E.eval(x, f); f = -f.div(m);
  mstreal mke = 0;
  if (verbose) cout << "Optim::langevinDynamics, initial potential energy = " << en << endl;
  for (int i = 0; i < numIters; i++) {
    switch (integrator) {
      case 0:
        // Brooks-Brunger-Karplus (BBK) integration scheme, appropriate for small-gamma regime
        for (int ii = 0; ii < df; ii++) {
          R[ii] =  MstUtils::randNormal(0, sqrt(2*gamma*kT/(timeStep*m[ii])));
          R1[ii] = MstUtils::randNormal(0, sqrt(2*gamma*kT/(timeStep*m[ii])));
        }

        C = V + (timeStep/2)*(f - gamma*V + R); // p(t + 1/2 t)
        x = x + timeStep*C;
        en = E.eval(x, f1); f1 = -f1.div(m);
        V = C + (timeStep/2)*(f1 - gamma*C + R1);
        break;

      case 1:
        // Bussi and Parrinello (May, 2018)
        for (int ii = 0; ii < df; ii++) {
          R[ii] = MstUtils::randNormal(0, 1);
          R1[ii] = MstUtils::randNormal(0, 1);
        }

        C = c1*V + c2.mult(R); // p(t+) from eq 12a
        x = x + timeStep*C + f*timeStep2/2;
        en = E.eval(x, f1); f1 = -f1.div(m);
        V = c1*(C + timeStep*(f + f1)/2) + c2.mult(R1);
        break;

      case 2:
        // a second-order Langevin integrator
        for (int ii = 0; ii < x0.size(); ii++) {
          R[ii] = MstUtils::randNormal(0, 1);
          R1[ii] = MstUtils::randNormal(0, 1);
        }

        C = (timeStep2/2)*(f - gamma*V) + timeStep32*sig.mult(0.5*R + (0.5/sqrt(3))*R1);
        x = x + timeStep * V + C;
        en = E.eval(x, f1); f1 = -f1.div(m);
        V = V + timeStep*(f + f1)/2 - timeStep*gamma*V + timeStep12*sig.mult(R) - gamma*C;
        break;

      default:
        MstUtils::error("unrecognized integrator type!", "Optim::langevinDynamics");
    }
    if (verbose) {
      mstreal ke = (m.mult(V.mult(V))/2).sum().norm(); // the last norm is to make a scalar
      mke = (mke*i + ke/df)/(i + 1); // per degree-of-freedom mean
      if (i % 100 == 0) cout << "Optim::langevinDynamics potential energy = " << en << ", kinetic energy = " << ke << " (running <kT> = " << 2*mke << ") total = " << en + ke << endl;
    }
    f = f1;
    if ((i + 1) % saveInterval == 0) {
      // save snapshot
      trajectory.push_back(x);
      ener.push_back(en);
    }
  }

  // a second-order Langevin integrator
  // Vector x0(E.guessPoint()); // starting positions (row vector)
  // Vector V(x0.size());       // velocities
  // Vector f(x0.size());       // forces divided by mass, at current point
  // Vector f1(x0.size());      // forces divided by mass, at next point
  // Vector C(x0.size());       // C(t) from equation 88
  // Vector x(x0);              // current positions
  // Vector m(x0.size());       // masses
  // Vector chi(x0.size());     // normally-distributed random numbers with zero mean and unit variance
  // Vector R1(x0.size());   // normally-distributed random numbers with zero mean and unit variance
  // Vector sig(x0.size());     // sqrt(2*kT*gamma/m)
  // mstreal timeStep2 = timeStep*timeStep;
  // mstreal timeStep12 = sqrt(timeStep);
  // mstreal timeStep32 = pow(timeStep12, 3);
  // if (x0.size() % masses.size() != 0) {
  //   MstUtils::error(MstUtils::toString(x0.size()) + "-dimensional system, but received " + MstUtils::toString(masses.size()) + " masses", "Optim::langevinDynamics");
  // }
  // int dim = x0.size() / masses.size(); // number of degrees of freedom per particle
  // int np = x0.size() / dim;
  // for (int i = 0; i < x0.size(); i += dim) {
  //   for (int k = 0; k < dim; k++) m[i + k] = masses[i / dim];
  // }
  // for (int i = 0; i < x0.size(); i++) sig[i] = sqrt(2*kT*gamma/m[i]);
  //
  // mstreal en = E.eval(x, f); f = -f.div(m);
  // if (verbose) cout << "Optim::langevinDynamics, initial potential energy = " << en << endl;
  // for (int i = 0; i < numIters; i++) {
  //   // generate two sets of independent normally-distributed random numbers
  //   for (int ii = 0; ii < x0.size(); ii++) {
  //     chi[ii] = MstUtils::randNormal(0, 1);
  //     theta[ii] = MstUtils::randNormal(0, 1);
  //   }
  //
  //   C = (timeStep2/2)*(f - gamma*V) + timeStep32*sig.mult(0.5*chi + (0.5/sqrt(3))*theta);
  //   x = x + timeStep * V + C;
  //   mstreal en = E.eval(x, f1); f1 = -f1.div(m);
  //   V = V + timeStep*(f + f1)/2 - timeStep*gamma*V + timeStep12*sig.mult(chi) - gamma*C;
  //   if (verbose) {
  //     mstreal ke = (m.mult(V.mult(V))/2).sum().norm(); // the last norm is to make a scalar
  //     cout << "Optim::langevinDynamics, potential energy = " << en << ", kinetic energy = " << ke << " (mean = " << ke/np << "), total = " << en + ke << endl;
  //   }
  //   f = f1;
  //   if ((i + 1) % saveInterval == 0) {
  //     // save snapshot
  //     trajectory.push_back(x);
  //     ener.push_back(en);
  //   }
  // }

  return ener;
}

// mstreal Optim::lineSearch(optimizerEvaluator& E, const vector<mstreal>& point, vector<mstreal>& solution, const Vector& _dir, mstreal startStepSize, bool verbose) {
//   Vector x(point);
//   Vector midGrad(x.length());
//   Vector dir = _dir;
//   Vector curr = x, low = x, mid = x, high = x;
//   mstreal v;
//   if (dir.length() == 0) {
//     v = E.eval(x, dir); // set the line-search direction as the gradient
//   } else {
//     v = E.eval(x);
//   }
//   if (startStepSize <= 0) {
//     // startStepSize = MstUtils::min(0.1/(dir.abs()).max(), 1.0);
//     startStepSize = 0.1/(dir.abs()).max();
//   }
//   mstreal gamma = startStepSize;
//
//   // bracket the minimum between low and high
//   while (true) {
//     // make step
//     curr = x + gamma*dir;
//     if (E.eval(curr) < v) {
//       low = mid;
//       mid = curr;
//       gamma += startStepSize;
//     } else {
//       high = curr;
//       break;
//     }
//   }
//
//   // seaerch in the bracketed range
//   for (int i = 0; i < 10; i++) {
//     v = E.eval(mid, midGrad);
//     if (midGrad.dot(dir) < 0) {
//       // sub-divide the right side
//       curr = (mid + high)/2;
//       if (E.eval(curr) < v) {
//         low = mid;
//         mid = curr;
//       } else {
//         high = curr;
//       }
//     } else {
//       // sub-divide the left side
//       curr = (low + mid)/2;
//       if (E.eval(curr) < v) {
//         high = mid;
//         mid = curr;
//       } else {
//         low = curr;
//       }
//     }
//   }
//   solution = mid;
//   return E.eval(mid);
// }
