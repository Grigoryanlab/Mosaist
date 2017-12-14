#ifndef _MSTGEOMETRY_H
#define _MSTGEOMETRY_H

/* More complex geometry than in "msttypes.h", which often includes some linear
 * algebra and requires additional #includes and capabilities. */

#include "mstlinalg.h"
#include <chrono>

using namespace MST;

class MstGeometry {
  public:
    MstGeometry();
    /* gradient vector must be of length 6, and will be filled with partial
     * derivatives d(distance)/dc, where c runs over x, y, and z of the first
     * atom, then second atom. */
    template <class T>
    static mstreal distance(const CartesianPoint& atom1, const CartesianPoint& atom2, T& grad);

    /* gradient vector must be of length 9, and will be filled with partial
     * derivatives d(distance)/dc, where c runs over x, y, and z of the first
     * atom, then second atom, then third atom. */
    template <class T>
    static mstreal angle(const CartesianPoint& atom1, const CartesianPoint& atom2, const CartesianPoint& atom3, T& grad);

    /* gradient vector must be of length 12, and will be filled with partial
     * derivatives d(distance)/dc, where c runs over x, y, and z of the first
     * atom, then second atom, then third atom, then fourth atom. */
    template <class T>
    static mstreal dihedral(const CartesianPoint& atom1, const CartesianPoint& atom2, const CartesianPoint& atom3, const CartesianPoint& atom4, T& grad);

    template <class T>
    mstreal qcpRMSD(const T& A, const T& B, bool setTransform = false);

    /* Tests implementation of analytical gradients of bond, angle, and dihedral
     * using finite difference for comparison. */
    static bool testPrimitiveGradients();

    /* Tests QCP implementation for RMSD and RMSD gradient calculation. */
    static bool testQCP();

  protected:
    template <class T>
    static T& ref(T& obj) { return obj; }
    template <class T>
    static T& ref(T* obj) { return *obj; }

  private:
    mstreal rot[3][3];
    mstreal trans[3];
};

// see derivation in http://grigoryanlab.org/docs/dynamics_derivatives.pdf
template <class T>
mstreal MstGeometry::distance(const CartesianPoint& atom1, const CartesianPoint& atom2, T& grad) {
  mstreal x12 = atom1.getX() - atom2.getX();
  mstreal y12 = atom1.getY() - atom2.getY();
  mstreal z12 = atom1.getZ() - atom2.getZ();
  mstreal d = sqrt(x12*x12 + y12*y12 + z12*z12);
  if (MstUtils::closeEnough(d, 0.0)) {
    // by convention, set gradient to unity when the two atoms are on top of each other
    for (int i = 0; i < grad.size(); i++) grad[i] = 1.0;
  } else {
    grad[0] = x12/d;
    grad[1] = y12/d;
    grad[2] = z12/d;
    grad[3] = -grad[0];
    grad[4] = -grad[1];
    grad[5] = -grad[2];
  }
  return d;
}

// see derivation in http://grigoryanlab.org/docs/dynamics_derivatives.pdf
template <class T>
mstreal MstGeometry::angle(const CartesianPoint& atom1, const CartesianPoint& atom2, const CartesianPoint& atom3, T& grad) {
  mstreal x12 = atom1.getX() - atom2.getX();
  mstreal y12 = atom1.getY() - atom2.getY();
  mstreal z12 = atom1.getZ() - atom2.getZ();
  mstreal x32 = atom3.getX() - atom2.getX();
  mstreal y32 = atom3.getY() - atom2.getY();
  mstreal z32 = atom3.getZ() - atom2.getZ();
  mstreal L1 = sqrt(x12*x12 + y12*y12 + z12*z12);
  mstreal L2 = sqrt(x32*x32 + y32*y32 + z32*z32);
  if (MstUtils::closeEnough(L1, 0.0) || MstUtils::closeEnough(L2, 0.0)) {
    // angle is not really defined, when two subsequent atoms coicide, so just
    // set the derivative to something for the sake of code stability
    for (int i = 0; i < grad.size(); i++) grad[i] = 0.0;
  }
  mstreal p = x12*x32 + y12*y32 + z12*z32;
  mstreal d = p/(L1*L2);
  if (MstUtils::closeEnough(fabs(d), 1.0)) {
    // signs are arbitrary in this case (set to + by convention)
    grad[0] = (sqrt(y12*y12 + z12*z12)/L1)/L1;
    grad[1] = (sqrt(x12*x12 + z12*z12)/L1)/L1;
    grad[2] = (sqrt(x12*x12 + y12*y12)/L1)/L1;
    grad[6] = (sqrt(y32*y32 + z32*z32)/L2)/L2;
    grad[7] = (sqrt(x32*x32 + z32*z32)/L2)/L2;
    grad[8] = (sqrt(x32*x32 + y32*y32)/L2)/L2;
  } else {
    mstreal L12 = L1*L2;
    mstreal C = -1/(sqrt(1 - d*d)*L12*L12);
    mstreal pL12i = p*L1/L2;
    mstreal pL21i = p*L2/L1;
    grad[0] = C*(x32*L12 - pL21i*x12);
    grad[1] = C*(y32*L12 - pL21i*y12);
    grad[2] = C*(z32*L12 - pL21i*z12);
    grad[6] = C*(x12*L12 - pL12i*x32);
    grad[7] = C*(y12*L12 - pL12i*y32);
    grad[8] = C*(z12*L12 - pL12i*z32);
  }

  // gardient with respect to the middle point is always the minus sum of that
  // of the first and the third points
  grad[3] = -(grad[0] + grad[6]);
  grad[4] = -(grad[1] + grad[7]);
  grad[5] = -(grad[2] + grad[8]);
  return acos(d);
}

// see derivation in http://grigoryanlab.org/docs/dynamics_derivatives.pdf
template <class T>
mstreal MstGeometry::dihedral(const CartesianPoint& atom1, const CartesianPoint& atom2, const CartesianPoint& atom3, const CartesianPoint& atom4, T& grad) {
  mstreal x21 = atom2.getX() - atom1.getX();
  mstreal y21 = atom2.getY() - atom1.getY();
  mstreal z21 = atom2.getZ() - atom1.getZ();
  mstreal x32 = atom3.getX() - atom2.getX();
  mstreal y32 = atom3.getY() - atom2.getY();
  mstreal z32 = atom3.getZ() - atom2.getZ();
  mstreal x43 = atom4.getX() - atom3.getX();
  mstreal y43 = atom4.getY() - atom3.getY();
  mstreal z43 = atom4.getZ() - atom3.getZ();
  mstreal x31 = atom3.getX() - atom1.getX();
  mstreal y31 = atom3.getY() - atom1.getY();
  mstreal z31 = atom3.getZ() - atom1.getZ();
  mstreal x42 = atom4.getX() - atom2.getX();
  mstreal y42 = atom4.getY() - atom2.getY();
  mstreal z42 = atom4.getZ() - atom2.getZ();

  CartesianPoint N1(z21*y32 - y21*z32, x21*z32 - z21*x32, y21*x32 - x21*y32);
  CartesianPoint N2(y43*z32 - z43*y32, z43*x32 - x43*z32, x43*y32 - y43*x32);
  CartesianPoint angleGrad(9);
  mstreal th = MstGeometry::angle(N1, CartesianPoint(0, 0, 0), N2, angleGrad);

  grad[0]  = -angleGrad[1]*z32 + angleGrad[2]*y32;
  grad[1]  =  angleGrad[0]*z32 - angleGrad[2]*x32;
  grad[2]  = -angleGrad[0]*y32 + angleGrad[1]*x32;

  grad[3]  =  angleGrad[1]*z31 - angleGrad[2]*y31 - angleGrad[7]*z43 + angleGrad[8]*y43;
  grad[4]  = -angleGrad[0]*z31 + angleGrad[2]*x31 + angleGrad[6]*z43 - angleGrad[8]*x43;
  grad[5]  =  angleGrad[0]*y31 - angleGrad[1]*x31 - angleGrad[6]*y43 + angleGrad[7]*x43;

  grad[6]  = -angleGrad[1]*z21 + angleGrad[2]*y21 + angleGrad[7]*z42 - angleGrad[8]*y42;
  grad[7]  =  angleGrad[0]*z21 - angleGrad[2]*x21 - angleGrad[6]*z42 + angleGrad[8]*x42;
  grad[8]  = -angleGrad[0]*y21 + angleGrad[1]*x21 + angleGrad[6]*y42 - angleGrad[7]*x42;

  grad[9]  = -angleGrad[7]*z32 + angleGrad[8]*y32;
  grad[10] =  angleGrad[6]*z32 - angleGrad[8]*x32;
  grad[11] = -angleGrad[6]*y32 + angleGrad[7]*x32;

  if (N1 * CartesianPoint(x43, y43, z43) > 0) {
    for (int i = 0; i < grad.size(); i++) grad[i] = -grad[i];
    th = -th;
  }

  return th;
}

// QCP algorithm from: http://onlinelibrary.wiley.com/doi/10.1002/jcc.21439/epdf
template <class T>
mstreal MstGeometry::qcpRMSD(const T& A, const T& B, bool setTransform) {
  int i, j;
  int N = A.size();
  if (N != B.size()) MstUtils::error("structures are of different length", "MstGeometry::qcpRMSD");

  //compute centers for vector sets x, y
  mstreal cA[3], cB[3];
  for (i = 0; i < 3; i++) { cA[i] = 0.0; cB[i] = 0.0; }
  for (i = 0; i < N; i++) {
    cA[0] += A[i]->getX();
    cA[1] += A[i]->getY();
    cA[2] += A[i]->getZ();
    cB[0] += B[i]->getX();
    cB[1] += B[i]->getY();
    cB[2] += B[i]->getZ();
  }
  for (i = 0; i < 3; i++){
    cA[i] = cA[i]/N;
    cB[i] = cB[i]/N;
  }

  // compute the correlation matrix S and the inner products of the two structures
  mstreal S[3][3];
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) S[i][j] = 0.0;
  }
  mstreal ax, ay, az, bx, by, bz;
  mstreal GA = 0, GB = 0;
  for (i = 0; i < N; i++) {
    ax = A[i]->getX() - cA[0];
    ay = A[i]->getY() - cA[1];
    az = A[i]->getZ() - cA[2];
    bx = B[i]->getX() - cB[0];
    by = B[i]->getY() - cB[1];
    bz = B[i]->getZ() - cB[2];
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
  mstreal S2[3][3];
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) S2[i][j] = S[i][j]*S[i][j];
  }

  // calculate characteristic polynomial coefficients
  // NOTE: though there are repeated terms in F, G, H, and I, it actually turns
  // out to be better to let the compiler optimize this rathre than declare temp variables
  mstreal C2 = -2*(S2[0][0] + S2[0][1] + S2[0][2] + S2[1][0] + S2[1][1] + S2[1][2] + S2[2][0] + S2[2][1] + S2[2][2]);
  mstreal C1 = 8*(S[0][0]*S[1][2]*S[2][1] + S[1][1]*S[2][0]*S[0][2] + S[2][2]*S[0][1]*S[1][0] -
                  S[0][0]*S[1][1]*S[2][2] - S[1][2]*S[2][0]*S[0][1] - S[2][1]*S[1][0]*S[0][2]);
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
  mstreal tol = 10E-11;
  mstreal L = (GA + GB)/2; // this initial guess is key (see http://journals.iucr.org/a/issues/2005/04/00/sh5029/sh5029.pdf)
  mstreal Lold, L2, L3, L4, C22 = 2*C2;
  do {
    Lold = L;
    L2 = L*L;
    L3 = L2*L;
    L4 = L3*L;
    L = L - (L4 + C2*L2 + C1*L + C0)/(4*L3 + C22*L + C1);
  } while (fabs(L - Lold) > tol*L);

  // compute optimal rotation matrix
  if (setTransform) {
    mstreal K[4][4];
    K[0][0] =  S[0][0] + S[1][1] + S[2][2] - L;
    K[1][1] =  S[0][0] - S[1][1] - S[2][2] - L;
    K[2][2] = -S[0][0] + S[1][1] - S[2][2] - L;
    K[3][3] = -S[0][0] - S[1][1] + S[2][2] - L;
    K[0][1] = K[1][0] = S[1][2] - S[2][1];
    K[0][2] = K[2][0] = S[2][0] - S[0][2];
    K[0][3] = K[3][0] = S[0][1] - S[1][0];
    K[1][2] = K[2][1] = S[0][1] + S[1][0];
    K[1][3] = K[3][1] = S[0][2] + S[2][0];
    K[2][3] = K[3][2] = S[1][2] + S[2][1];

    mstreal evectol = 10E-11;
    mstreal q0, q1, q2, q3; // a column from the adjoint matrix
    q0 = K[1][1]*K[2][2]*K[3][3] - K[1][1]*K[2][3]*K[3][2] - K[1][2]*K[2][1]*K[3][3] + K[1][2]*K[2][3]*K[3][1] + K[1][3]*K[2][1]*K[3][2] - K[1][3]*K[2][2]*K[3][1];
    q1 = K[1][0]*K[2][3]*K[3][2] - K[1][0]*K[2][2]*K[3][3] + K[1][2]*K[2][0]*K[3][3] - K[1][2]*K[2][3]*K[3][0] - K[1][3]*K[2][0]*K[3][2] + K[1][3]*K[2][2]*K[3][0];
    q2 = K[1][0]*K[2][1]*K[3][3] - K[1][0]*K[2][3]*K[3][1] - K[1][1]*K[2][0]*K[3][3] + K[1][1]*K[2][3]*K[3][0] + K[1][3]*K[2][0]*K[3][1] - K[1][3]*K[2][1]*K[3][0];
    q3 = K[1][0]*K[2][2]*K[3][1] - K[1][0]*K[2][1]*K[3][2] + K[1][1]*K[2][0]*K[3][2] - K[1][1]*K[2][2]*K[3][0] - K[1][2]*K[2][0]*K[3][1] + K[1][2]*K[2][1]*K[3][0];
    mstreal qsqr = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
    if (qsqr < evectol) {
      q0 = K[0][1]*K[2][3]*K[3][2] - K[0][1]*K[2][2]*K[3][3] + K[0][2]*K[2][1]*K[3][3] - K[0][2]*K[2][3]*K[3][1] - K[0][3]*K[2][1]*K[3][2] + K[0][3]*K[2][2]*K[3][1];
      q1 = K[0][0]*K[2][2]*K[3][3] - K[0][0]*K[2][3]*K[3][2] - K[0][2]*K[2][0]*K[3][3] + K[0][2]*K[2][3]*K[3][0] + K[0][3]*K[2][0]*K[3][2] - K[0][3]*K[2][2]*K[3][0];
      q2 = K[0][0]*K[2][3]*K[3][1] - K[0][0]*K[2][1]*K[3][3] + K[0][1]*K[2][0]*K[3][3] - K[0][1]*K[2][3]*K[3][0] - K[0][3]*K[2][0]*K[3][1] + K[0][3]*K[2][1]*K[3][0];
      q3 = K[0][0]*K[2][1]*K[3][2] - K[0][0]*K[2][2]*K[3][1] - K[0][1]*K[2][0]*K[3][2] + K[0][1]*K[2][2]*K[3][0] + K[0][2]*K[2][0]*K[3][1] - K[0][2]*K[2][1]*K[3][0];
      qsqr = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
      if (qsqr < evectol) {
        q0 = K[0][1]*K[1][2]*K[3][3] - K[0][1]*K[1][3]*K[3][2] - K[0][2]*K[1][1]*K[3][3] + K[0][2]*K[1][3]*K[3][1] + K[0][3]*K[1][1]*K[3][2] - K[0][3]*K[1][2]*K[3][1];
        q1 = K[0][0]*K[1][3]*K[3][2] - K[0][0]*K[1][2]*K[3][3] + K[0][2]*K[1][0]*K[3][3] - K[0][2]*K[1][3]*K[3][0] - K[0][3]*K[1][0]*K[3][2] + K[0][3]*K[1][2]*K[3][0];
        q2 = K[0][0]*K[1][1]*K[3][3] - K[0][0]*K[1][3]*K[3][1] - K[0][1]*K[1][0]*K[3][3] + K[0][1]*K[1][3]*K[3][0] + K[0][3]*K[1][0]*K[3][1] - K[0][3]*K[1][1]*K[3][0];
        q3 = K[0][0]*K[1][2]*K[3][1] - K[0][0]*K[1][1]*K[3][2] + K[0][1]*K[1][0]*K[3][2] - K[0][1]*K[1][2]*K[3][0] - K[0][2]*K[1][0]*K[3][1] + K[0][2]*K[1][1]*K[3][0];
        qsqr = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
        if (qsqr < evectol) {
          q0 = K[0][1]*K[1][3]*K[2][2] - K[0][1]*K[1][2]*K[2][3] + K[0][2]*K[1][1]*K[2][3] - K[0][2]*K[1][3]*K[2][1] - K[0][3]*K[1][1]*K[2][2] + K[0][3]*K[1][2]*K[2][1];
          q1 = K[0][0]*K[1][2]*K[2][3] - K[0][0]*K[1][3]*K[2][2] - K[0][2]*K[1][0]*K[2][3] + K[0][2]*K[1][3]*K[2][0] + K[0][3]*K[1][0]*K[2][2] - K[0][3]*K[1][2]*K[2][0];
          q2 = K[0][0]*K[1][3]*K[2][1] - K[0][0]*K[1][1]*K[2][3] + K[0][1]*K[1][0]*K[2][3] - K[0][1]*K[1][3]*K[2][0] - K[0][3]*K[1][0]*K[2][1] + K[0][3]*K[1][1]*K[2][0];
          q3 = K[0][0]*K[1][1]*K[2][2] - K[0][0]*K[1][2]*K[2][1] - K[0][1]*K[1][0]*K[2][2] + K[0][1]*K[1][2]*K[2][0] + K[0][2]*K[1][0]*K[2][1] - K[0][2]*K[1][1]*K[2][0];
          qsqr = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
        }
      }
    }

    if (qsqr < evectol) {
      // if eigenvector still too small, then the optimal rotation is likely just the identity matrix
      rot[0][0] = 1; rot[0][1] = 0; rot[0][2] = 0;
      rot[1][0] = 0; rot[1][1] = 1; rot[1][2] = 0;
      rot[2][0] = 0; rot[2][1] = 0; rot[2][2] = 1;
    } else {
      mstreal normq = sqrt(qsqr);
      q0 /= normq; q1 /= normq; q2 /= normq; q3 /= normq;
      rot[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
      rot[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
      rot[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
      rot[0][1] = 2*(q1*q2 - q0*q3); rot[0][2] = 2*(q1*q3 + q0*q2);
      rot[1][0] = 2*(q1*q2 + q0*q3); rot[2][0] = 2*(q1*q3 - q0*q2);
      rot[1][2] = 2*(q2*q3 - q0*q1);
      rot[2][1] = 2*(q2*q3 + q0*q1);
    }
  }

  return sqrt((GA + GB - 2*L)/N);
}


#endif
