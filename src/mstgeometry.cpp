#include "mstgeometry.h"

bool MstGeometry::testPrimitiveGradients() {
  long int x = std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
  srand(x);

  // test a bunch of times
  CartesianPoint bondGrad(6, 0.0), angleGrad(9, 0.0), diheGrad(12, 0.0);
  CartesianPoint bondGradFD(6, 0.0), angleGradFD(9, 0.0), diheGradFD(12, 0.0);
  for (int k = 0; k < 100; k++) {
    mstreal L = 1 + MstUtils::randUnit()*10; // length scale
    mstreal del = 0.001; // finite-difference coordinate step
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

// QCP algorithm from: http://onlinelibrary.wiley.com/doi/10.1002/jcc.21439/epdf
mstreal MstGeometry::qcpRMSD(const vector<CartesianPoint> &A, const vector<CartesianPoint> &B) {
  int i, j;
  int N = A.size();
  if (N != B.size()) MstUtils::error("structures are of different length", "MstGeometry::qcpRMSD");

  //compute centers for vector sets x, y
  vector<mstreal> cA(3, 0.0), cB(3, 0.0);
  for (i = 0; i < N; i++){
    for (j = 0; j < 3; j++) {
      cA[j] += A[i][j];
      cB[j] += B[i][j];
    }
  }
  for (i = 0; i < 3; i++){
    cA[i] = cA[i]/N;
    cB[i] = cB[i]/N;
  }

  // compute the correlation matrix S and the inner products of the two structures
  vector<vector<mstreal> > S(3, vector<mstreal>(3, 0.0));
  mstreal ax, ay, az, bx, by, bz;
  mstreal GA = 0, GB = 0;
  for (i = 0; i < N; i++) {
    ax = A[i][0] - cA[0];
    ay = A[i][1] - cA[1];
    az = A[i][2] - cA[2];
    bx = B[i][0] - cB[0];
    by = B[i][1] - cB[1];
    bz = B[i][2] - cB[2];
    GA += ax*ax + ay*ay + az*az;
    GB += bx*bx + by*by + bz*bz;
    S[0][0] += bx * ax;
    S[0][1] += bx * ay;
    S[0][2] += bx * az;
    S[1][0] += by * ax;
    S[1][1] += by * ay;
    S[1][2] += by * az;
    S[2][0] += bz * ax;
    S[2][1] += bz * ay;
    S[2][2] += bz * az;
  }

  // square of S
  vector<vector<mstreal> > S2(3, vector<mstreal>(3, 0.0));
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) S2[i][j] = S[i][j]*S[i][j];
  }

  // calculate characteristic polynomial coefficients
  // NOTE: F, G, H, I could be made more efficient!
  mstreal C2 = -2*(S2[0][0] + S2[0][1] + S2[0][2] + S2[1][0] + S2[1][1] + S2[1][2] + S2[2][0] + S2[2][1] + S2[2][2]);
  mstreal C1 = 8*(S[0][0]*S[1][2]*S[2][1] + S[1][1]*S[2][0]*S[0][2] + S[2][2]*S[0][1]*S[1][0] -
                  S[0][0]*S[1][1]*S[2][2] - S[1][2]*S[2][0]*S[0][1] + S[2][0]*S[1][0]*S[0][2]);
  mstreal D = (S2[0][1] + S2[0][2] - S2[1][0] - S2[2][0]); D = D*D;
  mstreal E1 = -S2[0][0] + S2[1][1] + S2[2][2] + S2[1][2] + S2[2][1];
  mstreal E2 = 2*(S[1][1]*S[2][2] - S[1][2]*S[2][1]);
  mstreal E = (E1 - E2) * (E1 + E2);
  mstreal F = (-(S[0][2] + S[2][0])*(S[1][2] - S[2][1]) + (S[0][1] - S[1][0])*(S[0][0] - S[1][1] - S[2][2])) *
              (-(S[0][2] - S[2][0])*(S[1][2] + S[2][1]) + (S[0][1] - S[1][0])*(S[0][0] - S[1][1] + S[2][2]));
  mstreal G = (-(S[0][2] + S[2][0])*(S[1][2] + S[2][1]) - (S[0][1] + S[1][0])*(S[0][0] + S[1][1] - S[2][2])) *
              (-(S[0][2] - S[2][0])*(S[1][2] - S[2][1]) - (S[0][1] + S[1][0])*(S[0][0] + S[1][1] + S[2][2]));
  mstreal H = ( (S[0][1] + S[1][0])*(S[1][2] + S[2][1]) + (S[0][2] + S[2][0])*(S[0][0] - S[1][1] + S[2][2])) *
              (-(S[0][1] - S[1][0])*(S[1][2] - S[2][1]) + (S[0][2] + S[2][0])*(S[0][0] + S[1][1] + S[2][2]));
  mstreal I = ( (S[0][1] + S[1][0])*(S[1][2] - S[2][1]) + (S[0][2] - S[2][0])*(S[0][0] - S[1][1] - S[2][2])) *
              (-(S[0][1] - S[1][0])*(S[1][2] + S[2][1]) + (S[0][2] - S[2][0])*(S[0][0] + S[1][1] - S[2][2]));
  mstreal C0 = D + E + F + G + H + I;

  // now iterate Newton–Raphson method to find max eigenvalue
  mstreal tol = 10E-8;
  mstreal L = (GA + GB)/2; // this initial guess is key (see http://journals.iucr.org/a/issues/2005/04/00/sh5029/sh5029.pdf)
  mstreal Lold, L2, L3, L4, C22 = C2*2;
  do {
    Lold = L;
    L2 = L*L;
    L3 = L2*L;
    L4 = L3*L;
    L = L - (L4 + C2*L2 + C1*L + C0)/(4*L3 + C22*L + C1);
  } while (fabs(L - Lold) > tol);

  return sqrt((GA + GB - 2*L)/N);
}
