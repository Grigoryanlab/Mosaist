#include "mstgeometry.h"

bool MstGeometry::testPrimitiveGradients() {
  long int x = std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
  srand(x);

  // test a bunch of times
  CartesianPoint bondGrad(6, 0.0), angleGrad(9, 0.0), diheGrad(12, 0.0);
  CartesianPoint bondGradFD(6, 0.0), angleGradFD(9, 0.0), diheGradFD(12, 0.0);
  for (int k = 0; k < 100; k++) {
    mstreal L = 1 + MstUtils::randUnit()*10; // length scale
    mstreal del = 0.0001; // finite-difference coordinate step
    mstreal tol = 0.01;   // tollerated norm difference between analytical and FD

    // make some random points
    vector<CartesianPoint> P;
    for (int i = 0; i < 4; i++) {
      P.push_back(CartesianPoint(L*MstUtils::randUnit(), L*MstUtils::randUnit(), L*MstUtils::randUnit()));
    }

    // test bond derivatives
    MstGeometry::distance(P[0], P[1], bondGrad);
    bondGradFD[0] = P[1].distance(P[0] + CartesianPoint(del, 0, 0)) - P[1].distance(P[0] - CartesianPoint(del, 0, 0));
    bondGradFD[1] = P[1].distance(P[0] + CartesianPoint(0, del, 0)) - P[1].distance(P[0] - CartesianPoint(0, del, 0));
    bondGradFD[2] = P[1].distance(P[0] + CartesianPoint(0, 0, del)) - P[1].distance(P[0] - CartesianPoint(0, 0, del));
    bondGradFD[3] = P[0].distance(P[1] + CartesianPoint(del, 0, 0)) - P[0].distance(P[1] - CartesianPoint(del, 0, 0));
    bondGradFD[4] = P[0].distance(P[1] + CartesianPoint(0, del, 0)) - P[0].distance(P[1] - CartesianPoint(0, del, 0));
    bondGradFD[5] = P[0].distance(P[1] + CartesianPoint(0, 0, del)) - P[0].distance(P[1] - CartesianPoint(0, 0, del));
    bondGradFD /= 2*del;
    if ((bondGrad - bondGradFD).norm() > tol) {
      cout << "test FAILED for distance gradient, point set: " << MstUtils::vecToString(P, " | ") << endl;
      cout << "analytical: " << MstUtils::vecToString(bondGrad) << endl;
      cout << "numerical: " << MstUtils::vecToString(bondGradFD) << endl;
      return false;
    }

    // test angle derivatives
    MstGeometry::angle(P[0], P[1], P[2], angleGrad);
    angleGradFD[0] = CartesianGeometry::angle(P[0] + CartesianPoint(del, 0, 0), P[1], P[2], true) - CartesianGeometry::angle(P[0] - CartesianPoint(del, 0, 0), P[1], P[2], true);
    angleGradFD[1] = CartesianGeometry::angle(P[0] + CartesianPoint(0, del, 0), P[1], P[2], true) - CartesianGeometry::angle(P[0] - CartesianPoint(0, del, 0), P[1], P[2], true);
    angleGradFD[2] = CartesianGeometry::angle(P[0] + CartesianPoint(0, 0, del), P[1], P[2], true) - CartesianGeometry::angle(P[0] - CartesianPoint(0, 0, del), P[1], P[2], true);
    angleGradFD[3] = CartesianGeometry::angle(P[0], P[1] + CartesianPoint(del, 0, 0), P[2], true) - CartesianGeometry::angle(P[0], P[1] - CartesianPoint(del, 0, 0), P[2], true);
    angleGradFD[4] = CartesianGeometry::angle(P[0], P[1] + CartesianPoint(0, del, 0), P[2], true) - CartesianGeometry::angle(P[0], P[1] - CartesianPoint(0, del, 0), P[2], true);
    angleGradFD[5] = CartesianGeometry::angle(P[0], P[1] + CartesianPoint(0, 0, del), P[2], true) - CartesianGeometry::angle(P[0], P[1] - CartesianPoint(0, 0, del), P[2], true);
    angleGradFD[6] = CartesianGeometry::angle(P[0], P[1], P[2] + CartesianPoint(del, 0, 0), true) - CartesianGeometry::angle(P[0], P[1], P[2] - CartesianPoint(del, 0, 0), true);
    angleGradFD[7] = CartesianGeometry::angle(P[0], P[1], P[2] + CartesianPoint(0, del, 0), true) - CartesianGeometry::angle(P[0], P[1], P[2] - CartesianPoint(0, del, 0), true);
    angleGradFD[8] = CartesianGeometry::angle(P[0], P[1], P[2] + CartesianPoint(0, 0, del), true) - CartesianGeometry::angle(P[0], P[1], P[2] - CartesianPoint(0, 0, del), true);

    angleGradFD /= 2*del;
    if ((angleGrad - angleGradFD).norm() > tol) {
      cout << "test FAILED for angle gradient, point set: " << MstUtils::vecToString(P, " | ") << endl;
      cout << "analytical: " << MstUtils::vecToString(angleGrad) << endl;
      cout << "numerical: " << MstUtils::vecToString(angleGradFD) << endl;
      return false;
    }

    // test dihedral derivatives
    MstGeometry::dihedral(P[0], P[1], P[2], P[3], diheGrad);
    diheGradFD[0] = CartesianGeometry::dihedral(P[0] + CartesianPoint(del, 0, 0), P[1], P[2], P[3], true) - CartesianGeometry::dihedral(P[0] - CartesianPoint(del, 0, 0), P[1], P[2], P[3], true);
    diheGradFD[1] = CartesianGeometry::dihedral(P[0] + CartesianPoint(0, del, 0), P[1], P[2], P[3], true) - CartesianGeometry::dihedral(P[0] - CartesianPoint(0, del, 0), P[1], P[2], P[3], true);
    diheGradFD[2] = CartesianGeometry::dihedral(P[0] + CartesianPoint(0, 0, del), P[1], P[2], P[3], true) - CartesianGeometry::dihedral(P[0] - CartesianPoint(0, 0, del), P[1], P[2], P[3], true);
    diheGradFD[3] = CartesianGeometry::dihedral(P[0], P[1] + CartesianPoint(del, 0, 0), P[2], P[3], true) - CartesianGeometry::dihedral(P[0], P[1] - CartesianPoint(del, 0, 0), P[2], P[3], true);
    diheGradFD[4] = CartesianGeometry::dihedral(P[0], P[1] + CartesianPoint(0, del, 0), P[2], P[3], true) - CartesianGeometry::dihedral(P[0], P[1] - CartesianPoint(0, del, 0), P[2], P[3], true);
    diheGradFD[5] = CartesianGeometry::dihedral(P[0], P[1] + CartesianPoint(0, 0, del), P[2], P[3], true) - CartesianGeometry::dihedral(P[0], P[1] - CartesianPoint(0, 0, del), P[2], P[3], true);
    diheGradFD[6] = CartesianGeometry::dihedral(P[0], P[1], P[2] + CartesianPoint(del, 0, 0), P[3], true) - CartesianGeometry::dihedral(P[0], P[1], P[2] - CartesianPoint(del, 0, 0), P[3], true);
    diheGradFD[7] = CartesianGeometry::dihedral(P[0], P[1], P[2] + CartesianPoint(0, del, 0), P[3], true) - CartesianGeometry::dihedral(P[0], P[1], P[2] - CartesianPoint(0, del, 0), P[3], true);
    diheGradFD[8] = CartesianGeometry::dihedral(P[0], P[1], P[2] + CartesianPoint(0, 0, del), P[3], true) - CartesianGeometry::dihedral(P[0], P[1], P[2] - CartesianPoint(0, 0, del), P[3], true);
    diheGradFD[9] = CartesianGeometry::dihedral(P[0], P[1], P[2], P[3] + CartesianPoint(del, 0, 0), true) - CartesianGeometry::dihedral(P[0], P[1], P[2], P[3] - CartesianPoint(del, 0, 0), true);
    diheGradFD[10] = CartesianGeometry::dihedral(P[0], P[1], P[2], P[3] + CartesianPoint(0, del, 0), true) - CartesianGeometry::dihedral(P[0], P[1], P[2], P[3] - CartesianPoint(0, del, 0), true);
    diheGradFD[11] = CartesianGeometry::dihedral(P[0], P[1], P[2], P[3] + CartesianPoint(0, 0, del), true) - CartesianGeometry::dihedral(P[0], P[1], P[2], P[3] - CartesianPoint(0, 0, del), true);
    diheGradFD /= 2*del;
    if ((diheGrad - diheGradFD).norm() > tol) {
      cout << "test FAILED for dihedral gradient, point set: " << MstUtils::vecToString(P, " | ") << endl;
      cout << "analytical: " << MstUtils::vecToString(diheGrad) << endl;
      cout << "numerical: " << MstUtils::vecToString(diheGradFD) << endl;
      cout << "difference: " << (diheGrad - diheGradFD) << endl;
      cout << "norm difference: " << (diheGrad - diheGradFD).norm() << endl;
      return false;
    }
  }
  return true;
}

bool MstGeometry::testQCP() {
  long int x = std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
  srand(x);
  RMSDCalculator rc;

  // test a bunch of times
  long int Tkabsch = 0, Tqcp = 0;
  bool failed = false;
  for (int k = 0; k < 100; k++) {
    int N = MstUtils::randInt(100, 10); // number of atoms
    mstreal L = MstUtils::randUnit()*50; // length scale

    // create random atoms
    AtomPointerVector A(N), B(N);
    for (int i = 0; i < N; i++) {
      A[i] = new Atom(0, "X", MstUtils::randUnit()*L, MstUtils::randUnit()*L, MstUtils::randUnit()*L, 0, 0, false);
      B[i] = new Atom(0, "X", MstUtils::randUnit()*L, MstUtils::randUnit()*L, MstUtils::randUnit()*L, 0, 0, false);
    }

    // superimpose different ways
    Tkabsch -= std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
    mstreal rmsd = rc.bestRMSD(A, B);
    Tkabsch += std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
    Tqcp -= std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
    mstreal rmsdQCP = MstGeometry::qcpRMSD(A, B);
    Tqcp += std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();

    // test
    if (fabs(rmsd - rmsdQCP) > 10E-8) {
      failed = true;
      cout << "test FAILED for RMSD calculation, point sets:\n" << A << endl << B << endl;
      cout << "Kabsch: " << rmsd << endl;
      cout << "QCP: " << rmsdQCP << endl;
      cout << "difference: " << (rmsd - rmsdQCP) << endl;
    }
    A.deletePointers();
    B.deletePointers();
    if (failed) return false;
  }
  cout << "Kabsch: " << (10E9/Tkabsch) << " per second" << endl;
  cout << "KQcp:   " << (10E9/Tqcp) << " per second" << endl;
  return true;
}
