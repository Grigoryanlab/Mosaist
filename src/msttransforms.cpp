#include "msttransforms.h"

using namespace MST;

/* --------- Frame --------- */

Frame::Frame() {
  O[0] = 0; O[1] = 0; O[2] = 0;
  X[0] = 1; X[1] = 0; X[2] = 0;
  Y[0] = 0; Y[1] = 1; Y[2] = 0;
  Z[0] = 0; Z[1] = 0; Z[2] = 1;
}

Frame::Frame(CartesianPoint& _O, CartesianPoint& _X, CartesianPoint& _Y, CartesianPoint& _Z) {
  MstUtils::assert((_O.size() == 3) && (_X.size() == 3) && (_Y.size() == 3) && (_Z.size() == 3),
      "Frame class currently supports only 3D coordinate frames; specified origin and axes must be 3D vectors", "Frame::Frame(CartesianPoint&, CartesianPoint&, CartesianPoint&, CartesianPoint&)");
  real xn = _X.norm(); real yn = _Y.norm(); real zn = _Z.norm();
  O[0] = _O[0]; O[1] = _O[1]; O[2] = _O[2];
  X[0] = _X[0]/xn; X[1] = _X[1]/xn; X[2] = _X[2]/xn;
  Y[0] = _Y[0]/yn; Y[1] = _Y[1]/yn; Y[2] = _Y[2]/yn;
  Z[0] = _Z[0]/zn; Z[1] = _Z[1]/zn; Z[2] = _Z[2]/zn;
}

Frame::Frame(real _ox, real _oy, real _oz, real _xx, real _xy, real _xz, real _yx, real _yy, real _yz, real _zx, real _zy, real _zz) {
  real xn = sqrt(_xx*_xx + _xy*_xy + _xz*_xz);
  real yn = sqrt(_yx*_yx + _yy*_yy + _yz*_yz);
  real zn = sqrt(_zx*_zx + _zy*_zy + _zz*_zz);
  O[0] = _ox; O[1] = _oy; O[2] = _oz;
  X[0] = _xx/xn; X[1] = _xy/xn; X[2] = _xz/xn;
  Y[0] = _yx/yn; Y[1] = _yy/yn; Y[2] = _yz/yn;
  Z[0] = _zx/zn; Z[1] = _zy/zn; Z[2] = _zz/zn;
}

Frame::Frame(Frame& other) {
  O[0] = other.O[0]; O[1] = other.O[1]; O[2] = other.O[2];
  X[0] = other.X[0]; X[1] = other.X[1]; X[2] = other.X[2];
  Y[0] = other.Y[0]; Y[1] = other.Y[1]; Y[2] = other.Y[2];
  Z[0] = other.Z[0]; Z[1] = other.Z[1]; Z[2] = other.Z[2];
}

/* --------- Transform --------- */
Transform::Transform() {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      M[i][j] = (i == j) ? 1 : 0;
    }
  }
}

Transform::Transform(const Transform& other) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      M[i][j] = other.M[i][j];
    }
  }
}

Transform::Transform(CartesianPoint A, CartesianPoint B, CartesianPoint C, fillOrder order) {
  Transform(); // called default constructor to make the identity matrix
  if ((A.size() != 3) || (B.size() != 3) || (C.size() != 3)) {
    MstUtils::error("expected 3D vectors!", "Transform::Transform(CartesianPoint, CartesianPoint, CartesianPoint, fillOrder)");
  }
  CartesianPoint* points[3];
  points[0] = &A;
  points[1] = &B;
  points[2] = &C;
  for (int i = 0; i < 3; i++) {
    CartesianPoint* P = points[i];
    for (int j = 0; j < 3; j++) {
      if (order == fillOrder::byColumn) {
        (*this)(j, i) = (*P)[j];
      } else {
        (*this)(i, j) = (*P)[j];
      }
    }
  }
  // now fill out the last row and column
  for (int i = 0; i < 4; i++) {
    if (i == 3) {
      (*this)(i, 3) = 1;
    } else {
      (*this)(i, 3) = 0;
      (*this)(3, i) = 0;
    }
  }
}

Transform::Transform(CartesianPoint A, CartesianPoint B, CartesianPoint C, CartesianPoint D, fillOrder order) {
  if ((A.size() != 4) || (B.size() != 4) || (C.size() != 4) || (D.size() != 4)) {
    MstUtils::error("expected 4D vectors!", "Transform::Transform(CartesianPoint, CartesianPoint, CartesianPoint, CartesianPoint, fillOrder)");
  }
  CartesianPoint* points[4];
  points[0] = &A;
  points[1] = &B;
  points[2] = &C;
  points[3] = &D;
  for (int i = 0; i < 4; i++) {
    CartesianPoint* P = points[i];
    for (int j = 0; j < 4; j++) {
      if (order == fillOrder::byColumn) {
        (*this)(i, j) = (*P)[j];
      } else {
        (*this)(j, i) = (*P)[j];
      }
    }
  }
}

Transform::Transform(vector<real> trans) {
  Transform(); // called default constructor to make the identity matrix
  if (trans.size() != 3) MstUtils::error("expected a 3D vector!", "Transform::Transform(vector<real>)");
  for (int i = 0; i < 3; i++) {
    (*this)(i, 3) = trans[i];
  }
}

Transform::Transform(vector<vector<real> > rot) {
  if (rot.size() != 3) MstUtils::error("expected a 3x3 matrix!", "Transform::Transform(vector<vector<real> >)");
  Transform(CartesianPoint(rot[0]), CartesianPoint(rot[1]), CartesianPoint(rot[2]), byRow);
}

Transform::Transform(vector<vector<real> > rot, vector<real> trans) {
  if (rot.size() != 3) MstUtils::error("expected a 3x3 matrix!", "Transform::Transform(vector<vector<real> >, vector<real>)");
  if (trans.size() != 3) MstUtils::error("expected a 3D vector!", "Transform::Transform(vector<vector<real> >, vector<real>)");
  Transform(CartesianPoint(rot[0]), CartesianPoint(rot[1]), CartesianPoint(rot[2]), byRow);
  for (int i = 0; i < 3; i++) {
    (*this)(i, 3) = trans[i];
  }
}

real& Transform::operator()(int i, int j) {
  if ((i < 0) || (j < 0) || (i > 3) || (j > 3)) MstUtils::error("element " + MstUtils::toString(i) + " x " + MstUtils::toString(j) + " out of range", "Transform::perator()");
  return M[i][j];
}

real Transform::operator()(int i, int j) const {
  if ((i < 0) || (j < 0) || (i > 3) || (j > 3)) MstUtils::error("element " + MstUtils::toString(i) + " x " + MstUtils::toString(j) + " out of range", "Transform::perator()");
  return M[i][j];
}

const Transform Transform::operator*(const Transform& rhs) const {
  Transform T;
  // matrix multiply
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      T(i, j) = 0;
      for (int k = 0; k < 4; k++) {
        T(i, j) += (*this)(i, k) * rhs(k, j);
      }
    }
  }
  return T;
}

Transform& Transform::operator*=(const Transform& rhs) {
  Transform T = (*this) * rhs;
  *this = T;
  return *this;
}

const CartesianPoint Transform::operator*(const CartesianPoint& rhs) const {
  if (rhs.size() != 3) {
    MstUtils::error("Transform currently supports only 3D transforms, whereas a point of dimensionality " + MstUtils::toString(rhs.size()) + " was passed", "Transform::operator*(const CartesianPoint&&)");
  }
  CartesianPoint p(3, 0);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 4; j++) {
      p[i] += (*this)(i, j) * (j == 3 ? 1 : rhs[j]); // add W (homogeneous coordinate) of 1 to the point
    }
  }
  return p;
}

Transform Transform::inverse() {
  Transform& T = *this;
  Transform Ti;

  // hard-code a 4x4 inverse (taken from MATLAB)
  real det = (T(0,0)*T(1,1)*T(2,2)*T(3,3) - T(0,0)*T(1,1)*T(2,3)*T(3,2) - T(0,0)*T(1,2)*T(2,1)*T(3,3) + T(0,0)*T(1,2)*T(2,3)*T(3,1) + T(0,0)*T(1,3)*T(2,1)*T(3,2) - T(0,0)*T(1,3)*T(2,2)*T(3,1) - T(0,1)*T(1,0)*T(2,2)*T(3,3) + T(0,1)*T(1,0)*T(2,3)*T(3,2) + T(0,1)*T(1,2)*T(2,0)*T(3,3) - T(0,1)*T(1,2)*T(2,3)*T(3,0) - T(0,1)*T(1,3)*T(2,0)*T(3,2) + T(0,1)*T(1,3)*T(2,2)*T(3,0) + T(0,2)*T(1,0)*T(2,1)*T(3,3) - T(0,2)*T(1,0)*T(2,3)*T(3,1) - T(0,2)*T(1,1)*T(2,0)*T(3,3) + T(0,2)*T(1,1)*T(2,3)*T(3,0) + T(0,2)*T(1,3)*T(2,0)*T(3,1) - T(0,2)*T(1,3)*T(2,1)*T(3,0) - T(0,3)*T(1,0)*T(2,1)*T(3,2) + T(0,3)*T(1,0)*T(2,2)*T(3,1) + T(0,3)*T(1,1)*T(2,0)*T(3,2) - T(0,3)*T(1,1)*T(2,2)*T(3,0) - T(0,3)*T(1,2)*T(2,0)*T(3,1) + T(0,3)*T(1,2)*T(2,1)*T(3,0));
  if (fabs(det) < 10E-15) {
    MstUtils::warn("matrix too close to singular, determinant near zero", "Transform::inverse()");
  }
  Ti(0, 0) =  (T(1,1)*T(2,2)*T(3,3) - T(1,1)*T(2,3)*T(3,2) - T(1,2)*T(2,1)*T(3,3) + T(1,2)*T(2,3)*T(3,1) + T(1,3)*T(2,1)*T(3,2) - T(1,3)*T(2,2)*T(3,1))/det;
  Ti(0, 1) = -(T(0,1)*T(2,2)*T(3,3) - T(0,1)*T(2,3)*T(3,2) - T(0,2)*T(2,1)*T(3,3) + T(0,2)*T(2,3)*T(3,1) + T(0,3)*T(2,1)*T(3,2) - T(0,3)*T(2,2)*T(3,1))/det;
  Ti(0, 2) =  (T(0,1)*T(1,2)*T(3,3) - T(0,1)*T(1,3)*T(3,2) - T(0,2)*T(1,1)*T(3,3) + T(0,2)*T(1,3)*T(3,1) + T(0,3)*T(1,1)*T(3,2) - T(0,3)*T(1,2)*T(3,1))/det;
  Ti(0, 3) = -(T(0,1)*T(1,2)*T(2,3) - T(0,1)*T(1,3)*T(2,2) - T(0,2)*T(1,1)*T(2,3) + T(0,2)*T(1,3)*T(2,1) + T(0,3)*T(1,1)*T(2,2) - T(0,3)*T(1,2)*T(2,1))/det;
  Ti(1, 0) = -(T(1,0)*T(2,2)*T(3,3) - T(1,0)*T(2,3)*T(3,2) - T(1,2)*T(2,0)*T(3,3) + T(1,2)*T(2,3)*T(3,0) + T(1,3)*T(2,0)*T(3,2) - T(1,3)*T(2,2)*T(3,0))/det;
  Ti(1, 1) =  (T(0,0)*T(2,2)*T(3,3) - T(0,0)*T(2,3)*T(3,2) - T(0,2)*T(2,0)*T(3,3) + T(0,2)*T(2,3)*T(3,0) + T(0,3)*T(2,0)*T(3,2) - T(0,3)*T(2,2)*T(3,0))/det;
  Ti(1, 2) = -(T(0,0)*T(1,2)*T(3,3) - T(0,0)*T(1,3)*T(3,2) - T(0,2)*T(1,0)*T(3,3) + T(0,2)*T(1,3)*T(3,0) + T(0,3)*T(1,0)*T(3,2) - T(0,3)*T(1,2)*T(3,0))/det;
  Ti(1, 3) =  (T(0,0)*T(1,2)*T(2,3) - T(0,0)*T(1,3)*T(2,2) - T(0,2)*T(1,0)*T(2,3) + T(0,2)*T(1,3)*T(2,0) + T(0,3)*T(1,0)*T(2,2) - T(0,3)*T(1,2)*T(2,0))/det;
  Ti(2, 0) =  (T(1,0)*T(2,1)*T(3,3) - T(1,0)*T(2,3)*T(3,1) - T(1,1)*T(2,0)*T(3,3) + T(1,1)*T(2,3)*T(3,0) + T(1,3)*T(2,0)*T(3,1) - T(1,3)*T(2,1)*T(3,0))/det;
  Ti(2, 1) = -(T(0,0)*T(2,1)*T(3,3) - T(0,0)*T(2,3)*T(3,1) - T(0,1)*T(2,0)*T(3,3) + T(0,1)*T(2,3)*T(3,0) + T(0,3)*T(2,0)*T(3,1) - T(0,3)*T(2,1)*T(3,0))/det;
  Ti(2, 2) =  (T(0,0)*T(1,1)*T(3,3) - T(0,0)*T(1,3)*T(3,1) - T(0,1)*T(1,0)*T(3,3) + T(0,1)*T(1,3)*T(3,0) + T(0,3)*T(1,0)*T(3,1) - T(0,3)*T(1,1)*T(3,0))/det;
  Ti(2, 3) = -(T(0,0)*T(1,1)*T(2,3) - T(0,0)*T(1,3)*T(2,1) - T(0,1)*T(1,0)*T(2,3) + T(0,1)*T(1,3)*T(2,0) + T(0,3)*T(1,0)*T(2,1) - T(0,3)*T(1,1)*T(2,0))/det;
  Ti(3, 0) = -(T(1,0)*T(2,1)*T(3,2) - T(1,0)*T(2,2)*T(3,1) - T(1,1)*T(2,0)*T(3,2) + T(1,1)*T(2,2)*T(3,0) + T(1,2)*T(2,0)*T(3,1) - T(1,2)*T(2,1)*T(3,0))/det;
  Ti(3, 1) =  (T(0,0)*T(2,1)*T(3,2) - T(0,0)*T(2,2)*T(3,1) - T(0,1)*T(2,0)*T(3,2) + T(0,1)*T(2,2)*T(3,0) + T(0,2)*T(2,0)*T(3,1) - T(0,2)*T(2,1)*T(3,0))/det;
  Ti(3, 2) = -(T(0,0)*T(1,1)*T(3,2) - T(0,0)*T(1,2)*T(3,1) - T(0,1)*T(1,0)*T(3,2) + T(0,1)*T(1,2)*T(3,0) + T(0,2)*T(1,0)*T(3,1) - T(0,2)*T(1,1)*T(3,0))/det;
  Ti(3, 3) =  (T(0,0)*T(1,1)*T(2,2) - T(0,0)*T(1,2)*T(2,1) - T(0,1)*T(1,0)*T(2,2) + T(0,1)*T(1,2)*T(2,0) + T(0,2)*T(1,0)*T(2,1) - T(0,2)*T(1,1)*T(2,0))/det;

  return Ti;
}

CartesianPoint Transform::applyToCopy(CartesianPoint& p) {
  return (*this) * p;
}

void Transform::apply(CartesianPoint& p) {
  p = (*this) * p;
}

void Transform::apply(Atom* a) {
  CartesianPoint point(*a);
  point = (*this) * point;
  a->setCoor(point);
}

void Transform::apply(real& x, real& y, real& z) {
  CartesianPoint point(x, y, z);
  point = (*this) * point;
  x = point[0]; y = point[1]; z = point[2];
}

void Transform::apply(Residue* res) {
  for (int i = 0; i < res->atomSize(); i++) {
    apply((*res)[i]);
  }
}

void Transform::apply(Chain* chain) {
  for (int j = 0; j < chain->residueSize(); j++) {
    Residue& res = (*chain)[j];
    for (int i = 0; i < res.atomSize(); i++) {
      apply(res[i]);
    }
  }
}

void Transform::apply(Structure* S) {
  for (int k = 0; k < S->chainSize(); k++) {
    Chain& chain = (*S)[k];
    for (int j = 0; j < chain.residueSize(); j++) {
      Residue& res = chain[j];
      for (int i = 0; i < res.atomSize(); i++) {
        apply(res[i]);
      }
    }
  }
}


/* --------- TransformFactory --------- */
/* lots of valuable information was taken from http://web.cs.iastate.edu/~cs577/handouts/homogeneous-transform.pdf */

const real TransformFactory::degreesToRadians = M_PI/180.0;
const real TransformFactory::radiansToDegrees = 180.0/M_PI;

Transform TransformFactory::translate(real x, real y, real z) {
  Transform T;
  T(0, 3) = x; T(1, 3) = y; T(2, 3) = z;
  return T;
}

Transform TransformFactory::translate(CartesianPoint& p) {
  if (p.size() != 3) MstUtils::error("Translation vector of unexpected dimension " + MstUtils::toString(p.size()), "TransformFactory::translate");
  return translate(p[0], p[1], p[2]);
}

Transform TransformFactory::rotateAroundX(real angle) {
  return rotateAroundX(sin(angle * degreesToRadians), cos(angle * degreesToRadians));
}

Transform TransformFactory::rotateAroundY(real angle) {
  return rotateAroundY(sin(angle * degreesToRadians), cos(angle * degreesToRadians));
}

Transform TransformFactory::rotateAroundZ(real angle) {
  return rotateAroundZ(sin(angle * degreesToRadians), cos(angle * degreesToRadians));
}

Transform TransformFactory::rotateAroundX(real s, real c) {
  Transform T;
  T(1, 1) =  c; T(1, 2) = -s;
  T(2, 1) =  s; T(2, 2) =  c;
  return T;
}

Transform TransformFactory::rotateAroundY(real s, real c) {
  Transform T;
  T(0, 0) =  c; T(0, 2) =  s;
  T(2, 0) = -s; T(2, 2) =  c;
  return T;
}

Transform TransformFactory::rotateAroundZ(real s, real c) {
  Transform T;
  T(0, 0) =  c; T(0, 1) = -s;
  T(1, 0) =  s; T(1, 1) =  c;
  return T;
}

Transform TransformFactory::rotateAroundAxis(real u, real v, real w, real a) {
  if (sqrt(u*u + v*v) < 10E-10) {
    return rotateAroundZ(a);
  }

  // matrix for rotating vector (u,v,w) about Z-axis into XY-plane
  real uv = sqrt(u*u + v*v);
  Transform Rxz = rotateAroundZ(-v/uv, u/uv);

  // matrix for rotating vector (u,v,w) about Y-axis to the Z-axis
  real uvw = sqrt(u*u + v*v + w*w);
  Transform Rxz2z = rotateAroundY(-uv/uvw, w/uvw);

  // matrix for rotating around vector (u,v,w) by angle a
  return Rxz.inverse() * Rxz2z.inverse() * rotateAroundZ(a) * Rxz2z * Rxz;

}

Transform TransformFactory::rotateAroundAxis(CartesianPoint& p, real a) {
  if (p.size() != 3) MstUtils::error("Rotation axis of unexpected dimension " + MstUtils::toString(p.size()), "TransformFactory::rotateAroundOrigAxis");
  return rotateAroundAxis(p[0], p[1], p[2], a);
}

Transform TransformFactory::rotateAroundLine(real p1, real p2, real p3, real q1, real q2, real q3, real a) {
  real r1 = q1 - p1;
  real r2 = q2 - p2;
  real r3 = q3 - p3;

  real sinThx = r2/sqrt(r2*r2 + r3*r3);
  real cosThx = r3/sqrt(r2*r2 + r3*r3);
  real sinThy = r1/sqrt(r1*r1 + r2*r2 + r3*r3);
  real cosThy = sqrt((r2*r2 + r3*r3)/(r1*r1 + r2*r2 + r3*r3));

  return translate(p1, p2, p3) * rotateAroundX(-sinThx, cosThx) * rotateAroundY(sinThy, cosThy) *
         rotateAroundZ(sin(a * degreesToRadians), cos(a * degreesToRadians)) *
         rotateAroundY(-sinThy, cosThy) * rotateAroundX(sinThx, cosThx) * translate(-p1, -p2, -p3);
}

Transform TransformFactory::rotateAroundLine(CartesianPoint& p, CartesianPoint& q, real a) {
  if ((p.size() != 3) || (q.size() != 3)) {
    MstUtils::error("Rotation axes of unexpected dimension " + MstUtils::toString(p.size()) + " and " + MstUtils::toString(q.size()), "TransformFactory::rotateAroundLine");
  }
  return rotateAroundLine(p[0], p[1], p[2], q[0], q[1], q[2], a);
}

Transform TransformFactory::alignVectorWithXAxis(real u, real v, real w) {
//  return matRotY(90) x matAlignVectorWithZAxis(u, v, w);
  real vw = sqrt(v*v + w*w);
  real uvw = sqrt(u*u + v*v + w*w);
  if (vw == 0) {
    return Transform();
  }

  return rotateAroundZ(u/uvw, -vw/uvw) * rotateAroundX(v/vw, -w/vw);
}

Transform TransformFactory::alignVectorWithYAxis(real u, real v, real w) {
  real uw = sqrt(u*u + w*w);
  real uvw = sqrt(u*u + v*v + w*w);
  if (uw == 0) {
    return Transform();
  }

  return rotateAroundX(v/uvw, -uw/uvw) * rotateAroundY(w/uw, -u/uw);
}

Transform TransformFactory::alignVectorWithZAxis(real u, real v, real w) {
  real uv = sqrt(u*u + v*v);
  real uvw = sqrt(u*u + v*v + w*w);
  if (uv == 0) {
    return Transform();
  }

  return rotateAroundY(w/uvw, -uv/uvw) * rotateAroundZ(u/uv, -v/uv);
}

Transform TransformFactory::alignVectorWithXAxis(CartesianPoint& p) {
  if (p.size() != 3) MstUtils::error("Specified axis of unexpected dimension " + MstUtils::toString(p.size()), "TransformFactory::alignVectorWithXAxis");
  return alignVectorWithXAxis(p[0], p[1], p[2]);
}

Transform TransformFactory::alignVectorWithYAxis(CartesianPoint& p) {
  if (p.size() != 3) MstUtils::error("Specified axis of unexpected dimension " + MstUtils::toString(p.size()), "TransformFactory::alignVectorWithYAxis");
  return alignVectorWithYAxis(p[0], p[1], p[2]);
}

Transform TransformFactory::alignVectorWithZAxis(CartesianPoint& p) {
  if (p.size() != 3) MstUtils::error("Specified axis of unexpected dimension " + MstUtils::toString(p.size()), "TransformFactory::alignVectorWithZAxis");
  return alignVectorWithZAxis(p[0], p[1], p[2]);
}

Transform TransformFactory::switchFrames(Frame& _from, Frame& _to) {
  // the procedure will be simple enough:
  // 1. change of basis from _from to the standard basis (laboratory basis)
  // 2. inverse change of basis from _to to the standard basis
  // 3. translate from the origin of _from to the origin of _to
  Transform T1(_from.getX(), _from.getY(), _from.getZ(), Transform::byColumn);
  Transform T2(_to.getX(), _to.getY(), _to.getZ(), Transform::byColumn);
  CartesianPoint ori = _from.getO() - _to.getO();

  return translate(ori) * (T2.inverse()) * T1;
}
