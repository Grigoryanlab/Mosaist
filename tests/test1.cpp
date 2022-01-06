#include "msttypes.h"
#include "mstoptions.h"
#include "mstoptim.h"
#include "msttransforms.h"
#include <chrono>

//debugger
using namespace std;
using namespace MST;

class C2symmetryEvaluator : public optimizerEvaluator {
  public:
    C2symmetryEvaluator(const Structure& prot) {
      alignProtomer(prot);
    }

    void setReference(const Structure& protA, const Structure& protB) {
      ref = protA;
      for (int ci = 0; ci < protB.chainSize(); ci++) {
        ref.appendChain(new Chain(protB[ci]));
      }
      refAtoms = ref.getAtoms();
    }

    void unsetReference() { ref = Structure(); refAtoms = AtomPointerVector(); }

    vector<mstreal> guessPoint() {
      if (guess.empty()) return {MstUtils::randUnit()*360, MstUtils::randUnit()*5, MstUtils::randUnit()*5, MstUtils::randNormal()*10};
      return guess;
    }

    void setGuessPoint(const vector<mstreal>& _guess) { guess = _guess; }

    mstreal eval(const vector<mstreal>& point) {
      // construct B1 and B2 that are C2 symmetric
      buildImages(point);

      if (!refAtoms.empty()) {
        AtomPointerVector curr = atoms1;
        curr.push_back(atoms2);
        return rc.bestRMSD(curr, refAtoms);
      }

      // compute complementarity between B1 and B2 as roughly the van der Waals contact energy
      ProximitySearch ps1(atoms1, 15.0);
      mstreal vdw = 0.0, rmin = 3.0;
      mstreal rmin2 = rmin*rmin;
      for (Atom* a : atoms2) {
        vector<int> closeIndices = ps1.getPointsWithin(a, 0.0, 12.0);
        for (int i : closeIndices) {
          mstreal d2 = a->distance2(atoms1[i]);
          vdw += min(pow(rmin2/d2, 6) - 2*pow(rmin2/d2, 3), 100.0);
        }
      }

      // return -1.0*close + 10.0*tooClose + (atoms1.getGeometricCenter().distance(atoms2.getGeometricCenter()));
      return vdw;
    }

    mstreal eval(const vector<mstreal>& point, Vector& grad) { grad = finiteDifferenceGradient(point, {1.0, 1.0, 1.0, 0.01}); return eval(point); }

    Structure buildAssembly(const vector<mstreal>& point) {
      buildImages(point);
      Structure assembly = B1;
      for (int ci = 0; ci < B2.chainSize(); ci++) {
        assembly.appendChain(new Chain(B2[ci]));
      }
      return assembly;
    }

  protected:
    void buildImages(const vector<mstreal>& point) {
      atoms1.copyCoordinates(atoms);
      Transform T = TransformFactory::translate(0, point[3], 0) *
                    TransformFactory::rotateAroundZ(point[2]) *
                    TransformFactory::rotateAroundY(point[1]) *
                    TransformFactory::rotateAroundX(point[0]);
      T.apply(atoms1);
      atoms2.copyCoordinates(atoms1);
      TransformFactory::rotateAroundX(180).apply(atoms2);
    }

    Transform alignProtomer(const Structure& prot) {
      A = prot;
      // move to origin and align with principal axes (main along Z)
      CartesianPoint O = AtomPointerVector(A.getAtoms()).getGeometricCenter();
      Matrix axes = MstLinAlg::getPrincipalAxes(A.getAtoms());
      vector<int> axisOrder = {1, 2, 3};
      CartesianPoint X = CartesianPoint((vector<mstreal>) axes.column(axisOrder[0] - 1)).getUnit();
      CartesianPoint Y = CartesianPoint((vector<mstreal>) axes.column(axisOrder[1] - 1)).getUnit();
      CartesianPoint Z = CartesianPoint((vector<mstreal>) axes.column(axisOrder[2] - 1)).getUnit();
      // make sure system is right handed
      if ((X.cross(Y)).dot(Z) < 0) X = -X;
      // align along axes by switching frames
      CartesianPoint OO(0, 0, 0);
      Transform T = TransformFactory::switchFrames(Frame(OO, CartesianPoint(1, 0, 0), CartesianPoint(0, 1, 0), CartesianPoint(0, 0, 1)), Frame(OO, X, Y, Z)) *
                    TransformFactory::translate(-O);
      T.apply(A);
      // T = TransformFactory::translate(-AtomPointerVector(A.getAtoms()).getGeometricCenter()) * T;
      B1 = A;
      B2 = A;
      atoms = A.getAtoms();
      atoms1 = B1.getAtoms();
      atoms2 = B2.getAtoms();

      return T;
    }

  private:
    Structure A;
    Structure B1, B2; // for temporarily storing C2-symmetric images
    AtomPointerVector atoms, atoms1, atoms2;
    vector<mstreal> guess;
    Structure ref;
    AtomPointerVector refAtoms; // reference starting structure to fit to
    RMSDCalculator rc;
};

bool incrementCounter(vector<int>& count, const vector<int>& base) {
  if (count.empty()) {
    count = vector<int>(base.size(), 0);
    return true;
  }

  for (int i = count.size() - 1; i >= 0; i--) {
    count[i]++;
    if (count[i] < base[i]) return true;
    count[i] = 0;
  }
  return false;
}

// mstreal gridSearch(C2symmetryEvaluator& E, const vector<mstreal>& _lo, const vector<mstreal>& _hi, const vector<int>& num, vector<mstreal>& sol, int numFocus = 4) {
//   vector<mstreal> bestPoint;
//   mstreal bestScore = INFINITY;
//   mstreal f = 0.3; // focusing factor
//
//   vector<mstreal> lo = _lo, hi = _hi;
//   int numPars = lo.size();
//   for (int cyc = 0; cyc < numFocus; cyc++) {
//     cout << "cycle low = " << MstUtils::vecToString(lo) << endl;
//     cout << "cycle high = " << MstUtils::vecToString(hi) << endl;
//     // define grid points along each dimension
//     vector<vector<mstreal>> grid(numPars);
//     for (int i = 0; i < numPars; i++) {
//       mstreal bw = (hi[i] - lo[i])/num[i];
//       grid[i].resize(num[i]);
//       for (int k = 0; k < num[i]; k++) grid[i][k] = lo[i] + bw*(k + 0.5);
//       cout << "gird for paramter " << i << ": " << MstUtils::vecToString(grid[i]) << endl;
//     }
//
//     // exhaustively visit the grid and find best point
//     vector<int> gridPoint;
//     vector<mstreal> point(numPars);
//     int c = 0;
//     while (incrementCounter(gridPoint, num)) {
//       for (int i = 0; i < numPars; i++) point[i] = grid[i][gridPoint[i]];
//       mstreal score = E.eval(point);
//       if (score < bestScore) {
//         bestScore = score;
//         bestPoint = point;
//       }
//       if (c % 1000 == 0) {
//         cout << "point index " << MstUtils::vecToString(gridPoint) << ", value " << MstUtils::vecToString(point) << " got score " << score << "; best so far is " << bestScore << endl;
//       }
//       c++;
//     }
//     cout << "best point after cycle was " << MstUtils::vecToString(bestPoint) << " with score " << bestScore << endl;
//
//     // redefine low/high ranges arund the best point
//     for (int i = 0; i < numPars; i++) {
//       mstreal W = (hi[i] - lo[i]);
//       lo[i] = bestPoint[i] - W*f/2;
//       hi[i] = bestPoint[i] + W*f/2;
//     }
//   }
//
//   sol = bestPoint;
//   return bestScore;
// }
//
// mstreal gridMinSearch(C2symmetryEvaluator& E, const vector<mstreal>& _lo, const vector<mstreal>& _hi, const vector<int>& num, vector<mstreal>& sol) {
//   vector<mstreal> bestPoint;
//   mstreal bestScore = INFINITY;
//   mstreal f = 0.3; // focusing factor
//   int numFocus = 1;
//
//   vector<mstreal> lo = _lo, hi = _hi;
//   int numPars = lo.size();
//   for (int cyc = 0; cyc < numFocus; cyc++) {
//     cout << "cycle low = " << MstUtils::vecToString(lo) << endl;
//     cout << "cycle high = " << MstUtils::vecToString(hi) << endl;
//     // define grid points along each dimension
//     vector<vector<mstreal>> grid(numPars);
//     for (int i = 0; i < numPars; i++) {
//       mstreal bw = (hi[i] - lo[i])/num[i];
//       grid[i].resize(num[i]);
//       for (int k = 0; k < num[i]; k++) grid[i][k] = lo[i] + bw*(k + 0.5);
//       cout << "gird for paramter " << i << ": " << MstUtils::vecToString(grid[i]) << endl;
//     }
//
//     // exhaustively visit the grid; for each promising point minimize
//     vector<int> gridPoint;
//     vector<mstreal> point(numPars);
//     int c = 0;
//     while (incrementCounter(gridPoint, num)) {
//       for (int i = 0; i < numPars; i++) point[i] = grid[i][gridPoint[i]];
//       mstreal score = E.eval(point);
//       if (score < 0) {
//         cout << "point index " << MstUtils::vecToString(gridPoint) << ", value " << MstUtils::vecToString(point) << " has a promising score of " << score << "; best so far is " << bestScore << endl;;
//         E.setGuessPoint(point);
//         score = Optim::conjGradMin(E, point, 10000, 10E-8);
//         if (score < bestScore) {
//           cout << "\tafter minimization became point " << MstUtils::vecToString(point) << " with score " << score << " ==> new best" << endl;
//           bestScore = score;
//           bestPoint = point;
// E.buildAssembly(bestPoint).writePDB("/tmp/bestsofar.pdb");
//         }
//       }
//       c++;
//     }
//     cout << "best point after cycle was " << MstUtils::vecToString(bestPoint) << " with score " << bestScore << endl;
//
//     // redefine low/high ranges arund the best point
//     for (int i = 0; i < numPars; i++) {
//       mstreal W = (hi[i] - lo[i]);
//       lo[i] = bestPoint[i] - W*f/2;
//       hi[i] = bestPoint[i] + W*f/2;
//     }
//   }
//
//   sol = bestPoint;
//   return bestScore;
// }

int main(int argc, char** argv) {
  MstOptions op;
  op.setTitle("Optimizes a C2 assembly for complementarity, given an individual protomer structure. Options:");
  op.addOption("p", "protomer PDB file.", true);
  op.addOption("s", "a starting point for the second protomer (PDB file).");
  op.addOption("o", "output base.", true);
  op.setOptions(argc, argv);
  MstUtils::seedRandEngine();

  Structure unit(op.getString("p"));
  C2symmetryEvaluator E(unit);
  if (op.isGiven("s")) {
    cout << "setting an initial guess point based on provided reference..." << endl;
    E.setReference(unit, Structure(op.getString("s")));
    vector<mstreal> sol;
    mstreal score;
    for (int i = 0; i < 20; i++) {
      score = Optim::fminsearch(E, 10000, sol);
      E.setGuessPoint(sol);
      score = Optim::conjGradMin(E, sol, 10000, 10E-8);
      E.setGuessPoint(sol);
      cout << "improved to score " << score << endl;
    }
    E.setGuessPoint(sol);
    E.unsetReference();
  }

  vector<mstreal> sol;
  mstreal score;
  if (true) {
    for (int i = 0; i < 10; i++) {
      score = Optim::fminsearch(E, 10000, sol);
      E.setGuessPoint(sol);
      score = Optim::conjGradMin(E, sol, 10000, 10E-8);
      E.setGuessPoint(sol);
      cout << "improved to score " << score << endl;
    }
  // } else if (true) {
  //   score = gridMinSearch(E, {0, 0, 0, 5}, {360, 360, 360, 20}, {10, 10, 10, 10}, sol);
  } else {
    // do an initial search for a good starting point
    mstreal bestScore = INFINITY; vector<mstreal> bestStartPoint;
    for (int i = 0; i < 1000; i++) {
      score = Optim::fminsearch(E, 10000, sol);
      if (score < bestScore) {
        bestStartPoint = sol;
        bestScore = score;
      }
      cout << "found score " << score << ", best starting score so far is " << bestScore << endl;
    }

    // now iterate to find the best solution
    E.setGuessPoint(bestStartPoint);
    for (int i = 0; i < 100; i++) {
      score = Optim::conjGradMin(E, sol, 10000, 10E-8);
      E.setGuessPoint(sol);
      score = Optim::fminsearch(E, 10000, sol);
      E.setGuessPoint(sol);
      cout << "improved to score " << score << endl;
    }
  }


  cout << "best score is " << score << endl;
  cout << "best solution is " << MstUtils::vecToString(sol) << endl;
  Structure assembly = E.buildAssembly(sol);
  assembly.writePDB(op.getString("o") + ".ass.pdb");
}
