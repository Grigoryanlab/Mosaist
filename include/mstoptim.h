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

    // --- Gradient descent
    static mstreal gradDescent(optimizerEvaluator& E, vector<mstreal>& solution, int numIters = 1000, mstreal tol = 10E-8, bool verbose = false);

    // --- Conjugate-gradient with line search minimization
    static mstreal conjGradMin(optimizerEvaluator& E, vector<mstreal>& solution, int numIters = 1000, mstreal tol = 10E-8, bool verbose = false);
    static mstreal lineSearch(optimizerEvaluator& E, const vector<mstreal>& point, vector<mstreal>& solution, const Vector& dir = Vector(), mstreal startStepSize = 0.01, bool verbose = false);

    /* --- Langevin integrator (second-order), implemented as in eq. 98 of molecular_dynamics_2015.pdf
     *  + masses  -- the mass of every particle. This vector must be of length
     *               k times the number of dimensions (i.e., the size of E.guessPoint()),
     *               where k (integer) is the number of dimensions per particle
     *               (e.g., often this will be three). It is assumed that the
     *               coordinate vector (e.g., those in E.guessPoint()) stores
     *               coordinates of a given particle adjacent to one anothher and
     *               in the same order for all particles (e.g., x1, y1, z2, x2,
     *               y2, z2, ..., xn, yn, zn).
     * + saveInterval -- how frequently to save the trajectory. If negative,
     *                   or if not specified, will save only the last snapshot.
     *                   Returns a vector if snapshot energies and trajectory
     *                   will be appended with coordinate vectors, one for each snapshot.
     * + timeStep -- is the integration time step.
     * + gamma    -- is the friction coefficient in units of 1/time, where time
                     is expressed in whatever units timeStep is expressed in.
     * + kT       -- k x T, in units of energy that are: [units of masses] x
     *               [units of coordinate vector]^2 / [units of timeStep]^2.
     */
    static vector<mstreal> langevinDynamics(optimizerEvaluator& E, const vector<mstreal>& masses, mstreal timeStep, mstreal gamma, mstreal kT, int numIters, vector<vector<mstreal> >& trajectory, int saveInterval = -1, bool verbose = false);
};

}
#endif
