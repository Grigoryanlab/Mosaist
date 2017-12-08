#ifndef _MSTGEOMETRY_H
#define _MSTGEOMETRY_H

/* More complex geometry than in "msttypes.h", which often includes some linear
 * algebra and requires additional #includes and capabilities. */

#include "mstlinalg.h"

class MstGeometry {
  public:
    /* gradient vector must be of length 6, and will be filled with partial
     * derivatives d(distance)/dc, where c runs over x, y, and z of the first
     * atom, then second atom. */
    template <class T>
    static mstreal distnace(const CartesianPoint& atom1, const CartesianPoint& atom2, T& grad);

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

};

// see derivation in http://grigoryanlab.org/docs/dynamics_derivatives.pdf
template <class T>
mstreal MstGeometry::distnace(const CartesianPoint& atom1, const CartesianPoint& atom2, T& grad) {
  mstreal x12 = atom1.getX() - atom2.getX();
  mstreal y12 = atom1.getY() - atom2.getY();
  mstreal z12 = atom1.getZ() - atom2.getZ();
  mstreal d = sqrt(x12*x12 + y12*y12 + z12*z12);
  if (MstUtils::closeEnough(d, 0)) {
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
  if (MstUtils::closeEnough(L1, 0) || MstUtils::closeEnough(L2, 0)) {
    // angle is not really defined, when two subsequent atoms coicide, so just
    // set the derivative to something for the sake of code stability
    for (int i = 0; i < grad.size(); i++) grad[i] = 0.0;
  }
  msteral p = x12*x32 + y12*y32 + z12*z32;
  mstreal d = p/(L1*L2);
  if (MstUtils::closeEnough(fabs(d), 1)) {
    // signs are arbitrary in this case (set to + by convention)
    grad[0] = (sqrt(y12*y12 + z12*z12)/L1)/L1;
    grad[1] = (sqrt(x12*x12 + z12*z12)/L1)/L1;
    grad[2] = (sqrt(x12*x12 + y12*y12)/L1)/L1;
    grad[6] = (sqrt(y32*y32 + z32*z32)/L2)/L2;
    grad[7] = (sqrt(x32*x32 + z32*z32)/L2)/L2;
    grad[8] = (sqrt(x32*x32 + y32*y32)/L2)/L2;
  } else {
    mstreal L12 = L1*L2;
    mstreal C = 1/(sqrt(1 - d*d)*L12*L12);
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
mstreal MstGeometry::dihedral(const CartesianPoint& atom1, const CartesianPoint& atom2, const CartesianPoint& atom3, const CartesianPoint& atom3, T& grad) {
  mstreal x12 = atom1.getX() - atom2.getX();
  mstreal y12 = atom1.getY() - atom2.getY();
  mstreal z12 = atom1.getZ() - atom2.getZ();
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

  CartesianPoint N1(y12*z32 - z12*z32, z12*x32 - x12*z32, x12*y32 - y12*x32);
  CartesianPoint O(0, 0, 0);
  CartesianPoint N2(y43*z32 - z43*y32, z43*x32 - x43*z32, x43*y32 - y43*x32);
  Vector angleGrad(9);
  mstreal th = MstGeometry::angle(N1, 0, N2, angleGrad);

  grad[0]  = -angleGrad[1]*z32 + angleGrad[2]*y32;
  grad[1]  =  angleGrad[0]*z32 - angleGrad[2]*x32;
  grad[2]  = -angleGrad[0]*y32 + angleGrad[1]*x32;

  grad[3]  =  angleGrad[1]*z31 - angleGrad[2]*y31 - angleGrad[4]*y43 + angleGrad[5]*y43;
  grad[4]  = -angleGrad[0]*z31 + angleGrad[2]*x31 + angleGrad[3]*z43 - angleGrad[5]*x43;
  grad[5]  =  angleGrad[0]*y31 - angleGrad[1]*x31 - angleGrad[3]*y43 + angleGrad[4]*x43;

  grad[6]  = -angleGrad[1]*z21 + angleGrad[2]*y21 + angleGrad[4]*z42 - angleGrad[5]*y42;
  grad[7]  =  angleGrad[0]*z21 - angleGrad[2]*x21 - angleGrad[3]*z42 + angleGrad[5]*x42;
  grad[8]  = -angleGrad[0]*y21 + angleGrad[1]*x21 + angleGrad[3]*y42 - angleGrad[4]*x42;

  grad[9]  = -angleGrad[4]*z32 + angleGrad[5]*y32;
  grad[10] =  angleGrad[3]*z32 - angleGrad[5]*x32;
  grad[11] = -angleGrad[3]*y32 + angleGrad[4]*x32;

  if (N1 * CartesianPoint(x43, y43, z43) > 0) {
    for (int i = 0; i < grad.size(); i++) grad[i] = -grad[i];
    th = -th;
  }

  return th;
}

#endif
