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
Transform::Transform(const Transform& other) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      M[i][j] = other.M[i][j];
    }
  }
}

void Transform::fill(vector<real>& A, vector<real>& B, vector<real>& C, vector<real>& D, fillOrder order) {
  if ((A.size() != 4) || (B.size() != 4) || (C.size() != 4) || (D.size() != 4)) {
    MstUtils::error("expected 4D vectors!", "Transform::Transform(CartesianPoint, CartesianPoint, CartesianPoint, CartesianPoint, fillOrder)");
  }
  vector<real>* points[4];
  points[0] = &A;
  points[1] = &B;
  points[2] = &C;
  points[3] = &D;
  for (int i = 0; i < 4; i++) {
    vector<real>* P = points[i];
    for (int j = 0; j < 4; j++) {
      if (order == fillOrder::byColumn) {
        (*this)(i, j) = (*P)[j];
      } else {
        (*this)(j, i) = (*P)[j];
      }
    }
  }
}

void Transform::fill(vector<real>& A, vector<real>& B, vector<real>& C, fillOrder order) {
  makeIdentity();
  if ((A.size() != 3) || (B.size() != 3) || (C.size() != 3)) {
    MstUtils::error("expected 3D vectors!", "Transform::Transform(CartesianPoint, CartesianPoint, CartesianPoint, fillOrder)");
  }
  vector<real>* points[3];
  points[0] = &A;
  points[1] = &B;
  points[2] = &C;
  for (int i = 0; i < 3; i++) {
    vector<real>* P = points[i];
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

Transform::Transform(vector<real> trans) {
  makeIdentity();
  if (trans.size() != 3) MstUtils::error("expected a 3D vector!", "Transform::Transform(vector<real>)");
  for (int i = 0; i < 3; i++) {
    (*this)(i, 3) = trans[i];
  }
}

Transform::Transform(real* trans) {
  makeIdentity();
  for (int i = 0; i < 3; i++) (*this)(i, 3) = trans[i];
}

Transform::Transform(vector<vector<real> > rot) {
  if (rot.size() != 3) MstUtils::error("expected a 3x3 matrix!", "Transform::Transform(vector<vector<real> >)");
  fill(rot[0], rot[1], rot[2], byRow);
}

Transform::Transform(real** rot) {
  makeIdentity();
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*this)(i, j) = rot[i][j];
    }
  }
}

Transform::Transform(vector<vector<real> > rot, vector<real> trans) {
  if (rot.size() != 3) MstUtils::error("expected a 3x3 matrix!", "Transform::Transform(vector<vector<real> >, vector<real>)");
  if (trans.size() != 3) MstUtils::error("expected a 3D vector!", "Transform::Transform(vector<vector<real> >, vector<real>)");
  fill(rot[0], rot[1], rot[2], byRow);
  for (int i = 0; i < 3; i++) {
    (*this)(i, 3) = trans[i];
  }
}

Transform::Transform(real** rot, real* trans) {
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

real& Transform::operator()(int i, int j) {
  if ((i < 0) || (j < 0) || (i > 3) || (j > 3)) MstUtils::error("element " + MstUtils::toString(i) + " x " + MstUtils::toString(j) + " out of range", "Transform::operator()");
  return M[i][j];
}

real Transform::operator()(int i, int j) const {
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

void Transform::eulerAngles(real& x, real& y, real& z) {
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

void Transform::apply(const AtomPointerVector& vec) {
  for (int i = 0; i < vec.size(); i++) apply(vec[i]);
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

Transform TransformFactory::translate(const CartesianPoint& p) {
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

Transform TransformFactory::switchFrames(const Frame& _from, const Frame& _to) {
  // the procedure will be simple enough:
  // 1. change of basis from _from to the standard basis (laboratory basis)
  // 2. inverse change of basis from _to to the standard basis
  // 3. translate from the origin of _from to the origin of _to
  Transform T1(_from.getX(), _from.getY(), _from.getZ(), Transform::byColumn);
  Transform T2(_to.getX(), _to.getY(), _to.getZ(), Transform::byColumn);
  CartesianPoint ori = _from.getO() - _to.getO();

  return translate(ori) * (T2.inverse()) * T1;
}


/* --------- TransformRMSD --------- */
TransformRMSD::TransformRMSD() {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) C[i][j] = 0;
  }
}

void TransformRMSD::init(const AtomPointerVector& atoms) {
  // compute origin
  real Tr[3];
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

real TransformRMSD::getRMSD(Transform& T1, Transform& T2) {
  real rot = 0;
  // the translation component
  for (int i = 0; i < 3; i++) {
    rot += (T1(i, 3) - T2(i, 3)) * (T1(i, 3) - T2(i, 3));
  }

  // trace(R^{a,b} * C) (in the notation of Hildebrandt and co-workers, 10.1002/jcc.23513, eq. 3)
  Transform T = T1.rotation().inverse() * T2.rotation();
  real trace = C[0][0] + C[1][1] + C[2][2];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      trace -= T(i, j)*C[j][i];
    }
  }
  trace *= 2;

  return sqrt(rot + trace);
}

real TransformRMSD::getRMSD(Transform& T1) {
  Transform I;
  return getRMSD(T1, I);
}

/* --------- Matrix --------- */

Matrix::Matrix(int rows, int cols, real val) {
  if ((rows < 0) || (cols < 0)) MstUtils::error("invalid dimensions specified: " + MstUtils::toString(rows) + " x " + MstUtils::toString(cols), "Matrix::Matrix");
  M.resize(rows, vector<real>(cols, val));
}

Matrix::Matrix(const vector<real>& p, bool col) {
  if (col) {
    M.resize(p.size(), vector<real>(1));
    for (int i = 0; i < p.size(); i++) M[i][0] = p[i];
  } else {
    M.resize(1, p);
  }
}

int Matrix::size(bool dim) const {
  if (dim) {
    if (M.size() == 0) return 0;
    return M[0].size();
  }
  return M.size();
}

Matrix& Matrix::operator/=(const real& s) {
  for (int i = 0; i < M.size(); i++) {
    for (int j = 0; j < M[i].size(); j++) {
      M[i][j] /= s;
    }
  }
  return *this;
}

Matrix& Matrix::operator*=(const real& s) {
  for (int i = 0; i < M.size(); i++) {
    for (int j = 0; j < M[i].size(); j++) {
      M[i][j] *= s;
    }
  }
  return *this;
}

const Matrix Matrix::operator/(const real& s) const {
  Matrix R = *this;
  R /= s;
  return R;
}

const Matrix Matrix::operator*(const real& s) const {
  Matrix R = *this;
  R *= s;
  return R;
}

Matrix& Matrix::operator*=(const Matrix& P) {
  Matrix R = *this * P;
  *this = R;
  return *this;
}

const Matrix Matrix::operator*(const Matrix& P) const {
  if (this->size(1) != P.size(0)) MstUtils::error("matrix dimensions do not agree", "Matrix::operator*(Matrix&)");
  int n = this->size(0);
  int m = P.size(1);
  int p = this->size(1);
  Matrix R(n, m, 0);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      for (int k = 0; k < p; k++) {
        R[i][j] += (*this)[i][k] * P[k][j];
      }
    }
  }
  return R;
}

Matrix& Matrix::operator+=(const Matrix& P) {
  if ((this->size(0) != P.size(0)) || (this->size(1) != P.size(1))) {
    MstUtils::error("matrix dimensions do not agree", "Matrix::operator+=(Matrix&)");
  }
  for (int i = 0; i < M.size(); i++) {
    for (int j = 0; j < M[i].size(); j++) {
      M[i][j] += P[i][j];
    }
  }
  return *this;
}

const Matrix Matrix::operator+(const Matrix& P) const {
  Matrix R = *this;
  R += P;
  return R;
}

Matrix& Matrix::operator-=(const Matrix& P) {
  if ((this->size(0) != P.size(0)) || (this->size(1) != P.size(1))) {
    MstUtils::error("matrix dimensions do not agree", "Matrix::operator-=(Matrix&)");
  }
  for (int i = 0; i < M.size(); i++) {
    for (int j = 0; j < M[i].size(); j++) {
      M[i][j] -= P[i][j];
    }
  }
  return *this;
}

const Matrix Matrix::operator-(const Matrix& P) const {
  Matrix R = *this;
  R -= P;
  return R;
}

const Matrix Matrix::operator-() const {
  Matrix R = *this;
  R *= -1;
  return R;
}

Matrix Matrix::inverse() {
  if (size(0) != size(1)) MstUtils::error("inverse of non-square matrix requested", "Matrix::inverse()");
  int N = size(0);
  Matrix Mi = *this;
  real det;

  switch(N) {
    case 0:
      MstUtils::error("inverse of an empty matrix requested", "Matrix::inverse()");
      break;

    case 1:
      Mi[0][0] = 1/M[0][0];
      break;

    case 2:
      det = (M[0][0]*M[1][1] - M[0][1]*M[1][0]);
      if (fabs(det) < 10E-15) MstUtils::warn("matrix too close to singular, determinant near zero", "Matrix::inverse()");
      Mi[0][0] =  M[1][1]/det;
      Mi[0][1] = -M[0][1]/det;
      Mi[1][0] = -M[1][0]/det;
      Mi[1][1] =  M[0][0]/det;
      break;

    case 3:
      det = (M[0][0]*M[1][1]*M[2][2] - M[0][0]*M[1][2]*M[2][1] - M[0][1]*M[1][0]*M[2][2] + M[0][1]*M[1][2]*M[2][0] + M[0][2]*M[1][0]*M[2][1] - M[0][2]*M[1][1]*M[2][0]);
      if (fabs(det) < 10E-15) MstUtils::warn("matrix too close to singular, determinant near zero", "Matrix::inverse()");
      Mi[0][0] =  (M[1][1]*M[2][2] - M[1][2]*M[2][1])/det;
      Mi[0][1] = -(M[0][1]*M[2][2] - M[0][2]*M[2][1])/det;
      Mi[0][2] =  (M[0][1]*M[1][2] - M[0][2]*M[1][1])/det;
      Mi[1][0] = -(M[1][0]*M[2][2] - M[1][2]*M[2][0])/det;
      Mi[1][1] =  (M[0][0]*M[2][2] - M[0][2]*M[2][0])/det;
      Mi[1][2] = -(M[0][0]*M[1][2] - M[0][2]*M[1][0])/det;
      Mi[2][0] =  (M[1][0]*M[2][1] - M[1][1]*M[2][0])/det;
      Mi[2][1] = -(M[0][0]*M[2][1] - M[0][1]*M[2][0])/det;
      Mi[2][2] =  (M[0][0]*M[1][1] - M[0][1]*M[1][0])/det;
      break;

    case 4:
      det = (M[0][0]*M[1][1]*M[2][2]*M[3][3] - M[0][0]*M[1][1]*M[2][3]*M[3][2] - M[0][0]*M[1][2]*M[2][1]*M[3][3] + M[0][0]*M[1][2]*M[2][3]*M[3][1] + M[0][0]*M[1][3]*M[2][1]*M[3][2] - M[0][0]*M[1][3]*M[2][2]*M[3][1] - M[0][1]*M[1][0]*M[2][2]*M[3][3] + M[0][1]*M[1][0]*M[2][3]*M[3][2] + M[0][1]*M[1][2]*M[2][0]*M[3][3] - M[0][1]*M[1][2]*M[2][3]*M[3][0] - M[0][1]*M[1][3]*M[2][0]*M[3][2] + M[0][1]*M[1][3]*M[2][2]*M[3][0] + M[0][2]*M[1][0]*M[2][1]*M[3][3] - M[0][2]*M[1][0]*M[2][3]*M[3][1] - M[0][2]*M[1][1]*M[2][0]*M[3][3] + M[0][2]*M[1][1]*M[2][3]*M[3][0] + M[0][2]*M[1][3]*M[2][0]*M[3][1] - M[0][2]*M[1][3]*M[2][1]*M[3][0] - M[0][3]*M[1][0]*M[2][1]*M[3][2] + M[0][3]*M[1][0]*M[2][2]*M[3][1] + M[0][3]*M[1][1]*M[2][0]*M[3][2] - M[0][3]*M[1][1]*M[2][2]*M[3][0] - M[0][3]*M[1][2]*M[2][0]*M[3][1] + M[0][3]*M[1][2]*M[2][1]*M[3][0]);
      if (fabs(det) < 10E-15) MstUtils::warn("matrix too close to singular, determinant near zero", "Matrix::inverse()");
      Mi[0][0] =  (M[1][1]*M[2][2]*M[3][3] - M[1][1]*M[2][3]*M[3][2] - M[1][2]*M[2][1]*M[3][3] + M[1][2]*M[2][3]*M[3][1] + M[1][3]*M[2][1]*M[3][2] - M[1][3]*M[2][2]*M[3][1])/det;
      Mi[0][1] = -(M[0][1]*M[2][2]*M[3][3] - M[0][1]*M[2][3]*M[3][2] - M[0][2]*M[2][1]*M[3][3] + M[0][2]*M[2][3]*M[3][1] + M[0][3]*M[2][1]*M[3][2] - M[0][3]*M[2][2]*M[3][1])/det;
      Mi[0][2] =  (M[0][1]*M[1][2]*M[3][3] - M[0][1]*M[1][3]*M[3][2] - M[0][2]*M[1][1]*M[3][3] + M[0][2]*M[1][3]*M[3][1] + M[0][3]*M[1][1]*M[3][2] - M[0][3]*M[1][2]*M[3][1])/det;
      Mi[0][3] = -(M[0][1]*M[1][2]*M[2][3] - M[0][1]*M[1][3]*M[2][2] - M[0][2]*M[1][1]*M[2][3] + M[0][2]*M[1][3]*M[2][1] + M[0][3]*M[1][1]*M[2][2] - M[0][3]*M[1][2]*M[2][1])/det;
      Mi[1][0] = -(M[1][0]*M[2][2]*M[3][3] - M[1][0]*M[2][3]*M[3][2] - M[1][2]*M[2][0]*M[3][3] + M[1][2]*M[2][3]*M[3][0] + M[1][3]*M[2][0]*M[3][2] - M[1][3]*M[2][2]*M[3][0])/det;
      Mi[1][1] =  (M[0][0]*M[2][2]*M[3][3] - M[0][0]*M[2][3]*M[3][2] - M[0][2]*M[2][0]*M[3][3] + M[0][2]*M[2][3]*M[3][0] + M[0][3]*M[2][0]*M[3][2] - M[0][3]*M[2][2]*M[3][0])/det;
      Mi[1][2] = -(M[0][0]*M[1][2]*M[3][3] - M[0][0]*M[1][3]*M[3][2] - M[0][2]*M[1][0]*M[3][3] + M[0][2]*M[1][3]*M[3][0] + M[0][3]*M[1][0]*M[3][2] - M[0][3]*M[1][2]*M[3][0])/det;
      Mi[1][3] =  (M[0][0]*M[1][2]*M[2][3] - M[0][0]*M[1][3]*M[2][2] - M[0][2]*M[1][0]*M[2][3] + M[0][2]*M[1][3]*M[2][0] + M[0][3]*M[1][0]*M[2][2] - M[0][3]*M[1][2]*M[2][0])/det;
      Mi[2][0] =  (M[1][0]*M[2][1]*M[3][3] - M[1][0]*M[2][3]*M[3][1] - M[1][1]*M[2][0]*M[3][3] + M[1][1]*M[2][3]*M[3][0] + M[1][3]*M[2][0]*M[3][1] - M[1][3]*M[2][1]*M[3][0])/det;
      Mi[2][1] = -(M[0][0]*M[2][1]*M[3][3] - M[0][0]*M[2][3]*M[3][1] - M[0][1]*M[2][0]*M[3][3] + M[0][1]*M[2][3]*M[3][0] + M[0][3]*M[2][0]*M[3][1] - M[0][3]*M[2][1]*M[3][0])/det;
      Mi[2][2] =  (M[0][0]*M[1][1]*M[3][3] - M[0][0]*M[1][3]*M[3][1] - M[0][1]*M[1][0]*M[3][3] + M[0][1]*M[1][3]*M[3][0] + M[0][3]*M[1][0]*M[3][1] - M[0][3]*M[1][1]*M[3][0])/det;
      Mi[2][3] = -(M[0][0]*M[1][1]*M[2][3] - M[0][0]*M[1][3]*M[2][1] - M[0][1]*M[1][0]*M[2][3] + M[0][1]*M[1][3]*M[2][0] + M[0][3]*M[1][0]*M[2][1] - M[0][3]*M[1][1]*M[2][0])/det;
      Mi[3][0] = -(M[1][0]*M[2][1]*M[3][2] - M[1][0]*M[2][2]*M[3][1] - M[1][1]*M[2][0]*M[3][2] + M[1][1]*M[2][2]*M[3][0] + M[1][2]*M[2][0]*M[3][1] - M[1][2]*M[2][1]*M[3][0])/det;
      Mi[3][1] =  (M[0][0]*M[2][1]*M[3][2] - M[0][0]*M[2][2]*M[3][1] - M[0][1]*M[2][0]*M[3][2] + M[0][1]*M[2][2]*M[3][0] + M[0][2]*M[2][0]*M[3][1] - M[0][2]*M[2][1]*M[3][0])/det;
      Mi[3][2] = -(M[0][0]*M[1][1]*M[3][2] - M[0][0]*M[1][2]*M[3][1] - M[0][1]*M[1][0]*M[3][2] + M[0][1]*M[1][2]*M[3][0] + M[0][2]*M[1][0]*M[3][1] - M[0][2]*M[1][1]*M[3][0])/det;
      Mi[3][3] =  (M[0][0]*M[1][1]*M[2][2] - M[0][0]*M[1][2]*M[2][1] - M[0][1]*M[1][0]*M[2][2] + M[0][1]*M[1][2]*M[2][0] + M[0][2]*M[1][0]*M[2][1] - M[0][2]*M[1][1]*M[2][0])/det;
      break;

    default:
      MstUtils::error("this is a basic matrix class that does not know how to invert matrices of dimension " + MstUtils::toString(N), "Matrix::inverse");
      break;
  }

  return Mi;
}

Matrix Matrix::transpose() {
  Matrix R = *this;
  for (int i = 0; i < M.size(); i++) {
    for (int j = 0; j < M[i].size(); j++) R[j][i] = M[i][j];
  }
  return R;
}

real Matrix::norm() {
  real n = 0;
  for (int i = 0; i < M.size(); i++) {
    for (int j = 0; j < M[i].size(); j++) n += M[i][j] * M[i][j];
  }
  return sqrt(n);
}

real Matrix::norm2() {
  real n = 0;
  for (int i = 0; i < M.size(); i++) {
    for (int j = 0; j < M[i].size(); j++) n += M[i][j] * M[i][j];
  }
  return n;
}
