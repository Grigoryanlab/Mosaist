#include "mstlinalg.h"

using namespace MST;

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
