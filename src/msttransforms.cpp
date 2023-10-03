#include "msttransforms.h"

using namespace MST;

/* --------- Frame --------- */

Frame::Frame() {
  O[0] = 0; O[1] = 0; O[2] = 0;
  X[0] = 1; X[1] = 0; X[2] = 0;
  Y[0] = 0; Y[1] = 1; Y[2] = 0;
  Z[0] = 0; Z[1] = 0; Z[2] = 1;
}

Frame::Frame(const CartesianPoint& _O, const CartesianPoint& _X, const CartesianPoint& _Y, const CartesianPoint& _Z) {
  constructFrame(_O,_X,_Y,_Z);
}

Frame::Frame(mstreal _ox, mstreal _oy, mstreal _oz, mstreal _xx, mstreal _xy, mstreal _xz, mstreal _yx, mstreal _yy, mstreal _yz, mstreal _zx, mstreal _zy, mstreal _zz) {
  mstreal xn = sqrt(_xx*_xx + _xy*_xy + _xz*_xz);
  mstreal yn = sqrt(_yx*_yx + _yy*_yy + _yz*_yz);
  mstreal zn = sqrt(_zx*_zx + _zy*_zy + _zz*_zz);
  O[0] = _ox; O[1] = _oy; O[2] = _oz;
  X[0] = _xx/xn; X[1] = _xy/xn; X[2] = _xz/xn;
  Y[0] = _yx/yn; Y[1] = _yy/yn; Y[2] = _yz/yn;
  Z[0] = _zx/zn; Z[1] = _zy/zn; Z[2] = _zz/zn;
}

Frame::Frame(const Frame& other) {
  O[0] = other.O[0]; O[1] = other.O[1]; O[2] = other.O[2];
  X[0] = other.X[0]; X[1] = other.X[1]; X[2] = other.X[2];
  Y[0] = other.Y[0]; Y[1] = other.Y[1]; Y[2] = other.Y[2];
  Z[0] = other.Z[0]; Z[1] = other.Z[1]; Z[2] = other.Z[2];
}

void Frame::setX(const CartesianPoint& _X) {
  X[0] = _X[0]; X[1] = _X[1]; X[2] = _X[2];
}

void Frame::setY(const CartesianPoint& _Y) {
  Y[0] = _Y[0]; Y[1] = _Y[1]; Y[2] = _Y[2];
}

void Frame::setZ(const CartesianPoint& _Z) {
  Z[0] = _Z[0]; Z[1] = _Z[1]; Z[2] = _Z[2];
}

void Frame::setO(const CartesianPoint& _O) {
  O[0] = _O[0]; O[1] = _O[1]; O[2] = _O[2];
}

void Frame::constructFrame(const CartesianPoint& _O, const CartesianPoint& _X, const CartesianPoint& _Y, const CartesianPoint& _Z) {
  MstUtils::assertCond((_O.size() == 3) && (_X.size() == 3) && (_Y.size() == 3) && (_Z.size() == 3),
      "Frame class currently supports only 3D coordinate frames; specified origin and axes must be 3D vectors", "Frame::Frame(CartesianPoint&, CartesianPoint&, CartesianPoint&, CartesianPoint&)");
  mstreal xn = _X.norm(); mstreal yn = _Y.norm(); mstreal zn = _Z.norm();
  O[0] = _O[0]; O[1] = _O[1]; O[2] = _O[2];
  X[0] = _X[0]/xn; X[1] = _X[1]/xn; X[2] = _X[2]/xn;
  Y[0] = _Y[0]/yn; Y[1] = _Y[1]/yn; Y[2] = _Y[2]/yn;
  Z[0] = _Z[0]/zn; Z[1] = _Z[1]/zn; Z[2] = _Z[2]/zn;
}

/* --------- Transform --------- */
Transform::Transform(const Transform& other) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      M[i][j] = other.M[i][j];
    }
  }
}

void Transform::fill(vector<mstreal>& A, vector<mstreal>& B, vector<mstreal>& C, vector<mstreal>& D, fillOrder order) {
  if ((A.size() != 4) || (B.size() != 4) || (C.size() != 4) || (D.size() != 4)) {
    MstUtils::error("expected 4D vectors!", "Transform::Transform(CartesianPoint, CartesianPoint, CartesianPoint, CartesianPoint, fillOrder)");
  }
  vector<mstreal>* points[4];
  points[0] = &A;
  points[1] = &B;
  points[2] = &C;
  points[3] = &D;
  for (int i = 0; i < 4; i++) {
    vector<mstreal>* P = points[i];
    for (int j = 0; j < 4; j++) {
      if (order == fillOrder::byColumn) {
        (*this)(i, j) = (*P)[j];
      } else {
        (*this)(j, i) = (*P)[j];
      }
    }
  }
}

void Transform::fill(vector<mstreal>& A, vector<mstreal>& B, vector<mstreal>& C, fillOrder order) {
  makeIdentity();
  if ((A.size() != 3) || (B.size() != 3) || (C.size() != 3)) {
    MstUtils::error("expected 3D vectors!", "Transform::Transform(CartesianPoint, CartesianPoint, CartesianPoint, fillOrder)");
  }
  vector<mstreal>* points[3];
  points[0] = &A;
  points[1] = &B;
  points[2] = &C;
  for (int i = 0; i < 3; i++) {
    vector<mstreal>* P = points[i];
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

Transform::Transform(CartesianPoint A, CartesianPoint B, CartesianPoint C, fillOrder order) {
  fill(A, B, C, order);
}

Transform::Transform(CartesianPoint A, CartesianPoint B, CartesianPoint C, CartesianPoint D, fillOrder order) {
  fill(A, B, C, D, order);
}

Transform::Transform(vector<mstreal> trans) {
  makeIdentity();
  if (trans.size() != 3) MstUtils::error("expected a 3D vector!", "Transform::Transform(vector<real>)");
  for (int i = 0; i < 3; i++) {
    (*this)(i, 3) = trans[i];
  }
}

Transform::Transform(mstreal* trans) {
  makeIdentity();
  for (int i = 0; i < 3; i++) (*this)(i, 3) = trans[i];
}

Transform::Transform(vector<vector<mstreal> > rot) {
  if (rot.size() != 3) MstUtils::error("expected a 3x3 matrix!", "Transform::Transform(vector<vector<real> >)");
  fill(rot[0], rot[1], rot[2], byRow);
}

Transform::Transform(mstreal** rot) {
  makeIdentity();
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*this)(i, j) = rot[i][j];
    }
  }
}

Transform::Transform(vector<vector<mstreal> > rot, vector<mstreal> trans) {
  if (rot.size() != 3) MstUtils::error("expected a 3x3 matrix!", "Transform::Transform(vector<vector<real> >, vector<real>)");
  if (trans.size() != 3) MstUtils::error("expected a 3D vector!", "Transform::Transform(vector<vector<real> >, vector<real>)");
  fill(rot[0], rot[1], rot[2], byRow);
  for (int i = 0; i < 3; i++) {
    (*this)(i, 3) = trans[i];
  }
}

Transform::Transform(mstreal** rot, mstreal* trans) {
  for (int i = 0; i < 3; i++) {
    (*this)(i, 3) = trans[i];
    (*this)(3, i) = 0;
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*this)(i, j) = rot[i][j];
    }
  }
  (*this)(3, 3) = 1;
}

void Transform::makeIdentity() {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      M[i][j] = (i == j) ? 1 : 0;
    }
  }
}

mstreal& Transform::operator()(int i, int j) {
  if ((i < 0) || (j < 0) || (i > 3) || (j > 3)) MstUtils::error("element " + MstUtils::toString(i) + " x " + MstUtils::toString(j) + " out of range", "Transform::operator()");
  return M[i][j];
}

mstreal Transform::operator()(int i, int j) const {
  if ((i < 0) || (j < 0) || (i > 3) || (j > 3)) MstUtils::error("element " + MstUtils::toString(i) + " x " + MstUtils::toString(j) + " out of range", "Transform::operator()");
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

  // invert the rotation part by taking the transpose
  for (int i = 0; i < 3; i++) {
    for (int j = i; j < 3; j++) {
      Ti(i, j) = T(j, i);
      Ti(j, i) = T(i, j);
    }
  }

  // if this is a pure rotation, then this will be enough; otherwise, invert the translation too
  if ((T(0, 3) != 0) || (T(1, 3) != 0) || (T(2, 3) != 0)) {
    Transform Tti;
    for (int i = 0; i < 3; i++) Tti(i, 3) = -T(i, 3);
    Ti =  Ti * Tti;
  }

  return Ti;
}

Transform Transform::rotation() {
  Transform T = *this;
  for (int i = 0; i < 3; i++) T(i, 3) = 0; // zero out the translation part
  return T;
}

Transform Transform::translation() {
  Transform T; // initializes to identity
  for (int i = 0; i < 3; i++) T(i, 3) = (*this)(i, 3);
  return T;
}

void Transform::eulerAngles(mstreal& x, mstreal& y, mstreal& z) {
  // check for gimbal lock
  if (MstUtils::closeEnough((float) M[2][0], (float) -1.0)) {
    float z = 0; // gimbal lock, value of x doesn't matter
    float y = M_PI / 2;
    float x = z + atan2(M[1][0], M[2][0]);
  } else if (MstUtils::closeEnough((float) M[2][0], (float) 1.0)) {
    float z = 0;
    float y = -M_PI / 2;
    float x = -z + atan2(-M[1][0], -M[2][0]);
  } else { // two solutions exist
    float y1 = -asin(M[2][0]);
    float y2 = M_PI - y1;

    float x1 = atan2(M[2][1] / cos(y1), M[2][2] / cos(y1));
    float x2 = atan2(M[2][1] / cos(y2), M[2][2] / cos(y2));

    float z1 = atan2(M[1][0] / cos(y1), M[0][0] / cos(y1));
    float z2 = atan2(M[1][0] / cos(y2), M[0][0] / cos(y2));

    // choose one solution to return
    // for example the "shortest" rotation
    if ((fabsf(x1) + fabsf(y1) + fabsf(z1)) <= (fabsf(x2) + fabsf(y2) + fabsf(z2))) {
      x = x1; y = y1; z = z1;
    } else {
      x = x2; y = y2; z = z2;
    }
  }
  x *= TransformFactory::radiansToDegrees;
  y *= TransformFactory::radiansToDegrees;
  z *= TransformFactory::radiansToDegrees;
}

CartesianPoint Transform::applyToCopy(CartesianPoint& p) {
  return (*this) * p;
}

void Transform::apply(mstreal& x, mstreal& y, mstreal& z) {
  mstreal p[3];
  for (int i = 0; i < 3; i++) {
    p[i] += (*this)(i, 0) * x;
    p[i] += (*this)(i, 1) * y;
    p[i] += (*this)(i, 2) * z;
    p[i] += (*this)(i, 3); // add W (homogeneous coordinate) of 1 to the point
  }
  x = p[0];
  y = p[1];
  z = p[2];
}

void Transform::apply(Frame& f) {
  CartesianPoint O = f.getO();
  CartesianPoint X = f.getX() + O;
  CartesianPoint Y = f.getY() + O;
  CartesianPoint Z = f.getZ() + O;
  O = (*this) * O;
  X = (*this) * X;
  Y = (*this) * Y;
  Z = (*this) * Z;
  f.setO(O);
  f.setX(X-O);
  f.setY(Y-O);
  f.setZ(Z-O);
}

void Transform::apply(CartesianPoint& p) {
  this->apply(p[0], p[1], p[2]);
}

void Transform::apply(Atom* a) {
  this->apply((*a)[0], (*a)[1], (*a)[2]);
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

void Transform::apply(const AtomPointerVector& vec) {
  for (int i = 0; i < vec.size(); i++) apply(vec[i]);
}

void Transform::write(ostream& _os) const {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) MstUtils::writeBin(_os, M[i][j]);
  }
}

void Transform::read(istream& _is) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) MstUtils::readBin(_is, M[i][j]);
  }
}


/* --------- TransformFactory --------- */
/* lots of valuable information was taken from http://web.cs.iastate.edu/~cs577/handouts/homogeneous-transform.pdf */

const mstreal TransformFactory::degreesToRadians = M_PI/180.0;
const mstreal TransformFactory::radiansToDegrees = 180.0/M_PI;

Transform TransformFactory::translate(mstreal x, mstreal y, mstreal z) {
  Transform T;
  T(0, 3) = x; T(1, 3) = y; T(2, 3) = z;
  return T;
}

Transform TransformFactory::translate(const CartesianPoint& p) {
  if (p.size() != 3) MstUtils::error("Translation vector of unexpected dimension " + MstUtils::toString(p.size()), "TransformFactory::translate");
  return translate(p[0], p[1], p[2]);
}

Transform TransformFactory::rotateAroundX(mstreal angle) {
  return rotateAroundX(sin(angle * degreesToRadians), cos(angle * degreesToRadians));
}

Transform TransformFactory::rotateAroundY(mstreal angle) {
  return rotateAroundY(sin(angle * degreesToRadians), cos(angle * degreesToRadians));
}

Transform TransformFactory::rotateAroundZ(mstreal angle) {
  return rotateAroundZ(sin(angle * degreesToRadians), cos(angle * degreesToRadians));
}

Transform TransformFactory::rotateAroundX(mstreal s, mstreal c) {
  Transform T;
  T(1, 1) =  c; T(1, 2) = -s;
  T(2, 1) =  s; T(2, 2) =  c;
  return T;
}

Transform TransformFactory::rotateAroundY(mstreal s, mstreal c) {
  Transform T;
  T(0, 0) =  c; T(0, 2) =  s;
  T(2, 0) = -s; T(2, 2) =  c;
  return T;
}

Transform TransformFactory::rotateAroundZ(mstreal s, mstreal c) {
  Transform T;
  T(0, 0) =  c; T(0, 1) = -s;
  T(1, 0) =  s; T(1, 1) =  c;
  return T;
}

Transform TransformFactory::rotateAroundAxis(mstreal u, mstreal v, mstreal w, mstreal a) {
  if (sqrt(u*u + v*v) < 10E-10) {
    return rotateAroundZ(a);
  }

  // matrix for rotating vector (u,v,w) about Z-axis into XY-plane
  mstreal uv = sqrt(u*u + v*v);
  Transform Rxz = rotateAroundZ(-v/uv, u/uv);

  // matrix for rotating vector (u,v,w) about Y-axis to the Z-axis
  mstreal uvw = sqrt(u*u + v*v + w*w);
  Transform Rxz2z = rotateAroundY(-uv/uvw, w/uvw);

  // matrix for rotating around vector (u,v,w) by angle a
  return Rxz.inverse() * Rxz2z.inverse() * rotateAroundZ(a) * Rxz2z * Rxz;

}

Transform TransformFactory::rotateAroundAxis(CartesianPoint& p, mstreal a) {
  if (p.size() != 3) MstUtils::error("Rotation axis of unexpected dimension " + MstUtils::toString(p.size()), "TransformFactory::rotateAroundOrigAxis");
  return rotateAroundAxis(p[0], p[1], p[2], a);
}

Transform TransformFactory::rotateAroundLine(mstreal p1, mstreal p2, mstreal p3, mstreal q1, mstreal q2, mstreal q3, mstreal a) {
  mstreal r1 = q1 - p1;
  mstreal r2 = q2 - p2;
  mstreal r3 = q3 - p3;

  mstreal sinThx = r2/sqrt(r2*r2 + r3*r3);
  mstreal cosThx = r3/sqrt(r2*r2 + r3*r3);
  mstreal sinThy = r1/sqrt(r1*r1 + r2*r2 + r3*r3);
  mstreal cosThy = sqrt((r2*r2 + r3*r3)/(r1*r1 + r2*r2 + r3*r3));

  return translate(p1, p2, p3) * rotateAroundX(-sinThx, cosThx) * rotateAroundY(sinThy, cosThy) *
         rotateAroundZ(sin(a * degreesToRadians), cos(a * degreesToRadians)) *
         rotateAroundY(-sinThy, cosThy) * rotateAroundX(sinThx, cosThx) * translate(-p1, -p2, -p3);
}

Transform TransformFactory::rotateAroundLine(CartesianPoint& p, CartesianPoint& q, mstreal a) {
  if ((p.size() != 3) || (q.size() != 3)) {
    MstUtils::error("Rotation axes of unexpected dimension " + MstUtils::toString(p.size()) + " and " + MstUtils::toString(q.size()), "TransformFactory::rotateAroundLine");
  }
  return rotateAroundLine(p[0], p[1], p[2], q[0], q[1], q[2], a);
}

Transform TransformFactory::alignVectorWithXAxis(mstreal u, mstreal v, mstreal w) {
//  return matRotY(90) x matAlignVectorWithZAxis(u, v, w);
  mstreal vw = sqrt(v*v + w*w);
  mstreal uvw = sqrt(u*u + v*v + w*w);
  if (vw == 0) {
    return Transform();
  }

  return rotateAroundZ(-vw/uvw, u/uvw) * rotateAroundX(-w/vw, v/vw);
}

Transform TransformFactory::alignVectorWithYAxis(mstreal u, mstreal v, mstreal w) {
  mstreal uw = sqrt(u*u + w*w);
  mstreal uvw = sqrt(u*u + v*v + w*w);
  if (uw == 0) {
    return Transform();
  }

  return rotateAroundX(-uw/uvw, v/uvw) * rotateAroundY(-u/uw, w/uw);
}

Transform TransformFactory::alignVectorWithZAxis(mstreal u, mstreal v, mstreal w) {
  mstreal uv = sqrt(u*u + v*v);
  mstreal uvw = sqrt(u*u + v*v + w*w);
  if (uv == 0) {
    return Transform();
  }

  return rotateAroundY(-uv/uvw, w/uvw) * rotateAroundZ(-v/uv, u/uv);
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

Transform TransformFactory::switchFrames(const Frame& _from, const Frame& _to) {
  /* The procedure will be simple enough:
   * 1. change of basis from _from to the standard basis (laboratory basis)
   * 2. adjust origin
   * 3. inverse change of basis from _to to the standard basis
   * Two things always confuse me at first, to explaining for clarity:
   * a. it is fromOrigin minus toOrigin because coordinates should change in the
   *    opposite direction of the origin change.
   * b. the update of origins is between the two change of bases because the
   *    coordinates of the origins in both to and from bases are relative to the
   *    laboratory frame. */
  Transform T1(_from.getX(), _from.getY(), _from.getZ(), Transform::byColumn);
  Transform T2(_to.getX(), _to.getY(), _to.getZ(), Transform::byColumn);
  CartesianPoint ori = _from.getO() - _to.getO();

  return  (T2.inverse()) * translate(ori) * T1;
}


/* --------- TransformRMSD --------- */
TransformRMSD::TransformRMSD() {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) C[i][j] = 0;
  }
}

void TransformRMSD::init(const AtomPointerVector& atoms) {
  // compute origin
  mstreal Tr[3];
  Tr[0] = 0; Tr[1] = 0; Tr[2] = 0;
  for (int i = 0; i < atoms.size(); i++) {
    Atom& a = *(atoms[i]);
    Tr[0] += a[0];
    Tr[1] += a[1];
    Tr[2] += a[2];
  }
  Tr[0] /= atoms.size();
  Tr[1] /= atoms.size();
  Tr[2] /= atoms.size();

  // compute covariance matrix
  for (int i = 0; i < 3; i++) {
    for (int j = i; j < 3; j++) {
      C[i][j] = 0;
      for (int k = 0; k < atoms.size(); k++) {
        Atom& a = *(atoms[k]);
        C[i][j] += (a[i] - Tr[i]) * (a[j] - Tr[j]);
      }
      C[i][j] /= atoms.size();
      C[j][i] = C[i][j];
    }
  }
}

void TransformRMSD::init(const Structure& S) {
  AtomPointerVector atoms = S.getAtoms();
  init(atoms);
}

mstreal TransformRMSD::getRMSD(Transform& T1, Transform& T2) {
  mstreal rot = 0;
  // the translation component
  for (int i = 0; i < 3; i++) {
    rot += (T1(i, 3) - T2(i, 3)) * (T1(i, 3) - T2(i, 3));
  }

  // trace(R^{a,b} * C) (in the notation of Hildebrandt and co-workers, 10.1002/jcc.23513, eq. 3)
  Transform T = T1.rotation().inverse() * T2.rotation();
  mstreal trace = C[0][0] + C[1][1] + C[2][2];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      trace -= T(i, j)*C[j][i];
    }
  }
  trace *= 2;

  return sqrt(rot + trace);
}

mstreal TransformRMSD::getRMSD(Transform& T1) {
  Transform I;
  return getRMSD(T1, I);
}
