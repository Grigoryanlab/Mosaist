#include "mstlinalg.h"

using namespace MST;

/* --------- Matrix --------- */

Matrix::~Matrix() {
  if (own) clear();
}

Matrix::Matrix(int rows, int cols, real val) {
  if ((rows < 0) || (cols < 0)) MstUtils::error("invalid dimensions specified: " + MstUtils::toString(rows) + " x " + MstUtils::toString(cols), "Matrix::Matrix");
  M.resize(rows, vector<real*>(cols, NULL));
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      M[i][j] = new real(val);
    }
  }
  own = true;
}

Matrix::Matrix(const vector<vector<real> >& _M) {
  int rows = _M.size();
  int cols = (_M.size() > 0) ? _M[0].size() : 0;
  M.resize(rows, vector<real*>(cols, NULL));
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      M[i][j] = new real(_M[i][j]);
    }
  }
  own = true;
}

Matrix::Matrix(const vector<vector<real*> >& _M, int rowBeg, int rowEnd, int colBeg, int colEnd) {
  if (rowEnd < 0) rowEnd = _M.size() - 1;
  if (colEnd < 0) colEnd = ((_M.size() > 0) ? _M[0].size() : 0) - 1;
  int rows = rowEnd - rowBeg + 1;
  int cols = colEnd - colBeg + 1;
  M.resize(rows, vector<real*>(cols, NULL));
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      M[i][j] = _M[i + rowBeg][j + colBeg];
    }
  }
  own = false;
}

void Matrix::clear() {
  for (int i = 0; i < M.size(); i++) {
    for (int j = 0; j < M[i].size(); j++) {
      delete(M[i][j]);
    }
  }
}

Matrix::Matrix(const Matrix& _M) {
  M.resize(_M.numRows(), vector<real*>(_M.numCols(), NULL));
  for (int i = 0; i < _M.numRows(); i++) {
    for (int j = 0; j < _M.numCols(); j++) {
      M[i][j] = new real(_M(i,j));
    }
  }
  own = true;
}

Matrix::Matrix(const vector<real>& p, bool col) {
  if (col) {
    M.resize(p.size(), vector<real*>(1));
    for (int i = 0; i < p.size(); i++) M[i][0] = new real(p[i]);
  } else {
    M.resize(1, vector<real*>(p.size()));
    for (int i = 0; i < p.size(); i++) M[0][i] = new real(p[i]);
  }
  own = true;
}

int Matrix::size(int dim) const {
  switch(dim) {
    case 1:
      return M.size();
    case 2:
      if (M.size() == 0) return 0;
      return M[0].size();
    default:
      MstUtils::error("out of range dimension " + MstUtils::toString(dim), "Matrix::size");
      return 0; // to make the compiler happy
  }
}

Matrix& Matrix::operator=(const Matrix& _M) {
  if ((numRows() != _M.numRows()) || (numCols() != _M.numCols())) {
    if (!getOwnFlag()) {
      MstUtils::error("dimensions must agree when assigning to a sub-Matrix", "Matrix::operator=");
    } else {
      clear();
      M.resize(_M.numRows(), vector<real*>(_M.numCols(), NULL));
      for (int i = 0; i < _M.numRows(); i++) {
        for (int j = 0; j < _M.numCols(); j++) {
          M[i][j] = new real(_M(i, j));
        }
      }
    }
  } else {
    for (int i = 0; i < this->numRows(); i++) {
      for (int j = 0; j < this->numCols(); j++) {
        *(M[i][j]) = _M(i, j);
      }
    }
  }
  return *this;
}

Matrix& Matrix::operator/=(const real& s) {
  for (int i = 0; i < M.size(); i++) {
    for (int j = 0; j < M[i].size(); j++) {
      *(M[i][j]) /= s;
    }
  }
  return *this;
}

Matrix& Matrix::operator*=(const real& s) {
  for (int i = 0; i < M.size(); i++) {
    for (int j = 0; j < M[i].size(); j++) {
      *(M[i][j]) *= s;
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
  if (this->size(2) != P.size(1)) MstUtils::error("matrix dimensions do not agree", "Matrix::operator*(Matrix&)");
  int n = this->size(1);
  int m = P.size(2);
  int p = this->size(2);
  Matrix R(n, m, 0);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      for (int k = 0; k < p; k++) {
        R(i,j) += (*this)(i,k) * P(k,j);
      }
    }
  }
  return R;
}

Matrix& Matrix::operator+=(const Matrix& P) {
  if ((this->size(1) != P.size(1)) || (this->size(2) != P.size(2))) {
    MstUtils::error("matrix dimensions do not agree", "Matrix::operator+=(Matrix&)");
  }
  for (int i = 0; i < M.size(); i++) {
    for (int j = 0; j < M[i].size(); j++) {
      *(M[i][j]) += *(P.M[i][j]);
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
  if ((this->size(1) != P.size(1)) || (this->size(2) != P.size(2))) {
    MstUtils::error("matrix dimensions do not agree", "Matrix::operator-=(Matrix&)");
  }
  for (int i = 0; i < M.size(); i++) {
    for (int j = 0; j < M[i].size(); j++) {
      *(M[i][j]) -= *(P.M[i][j]);
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

Matrix Matrix::row(int i) {
  Matrix subMat(M, i, i);
  return subMat;
}

Matrix Matrix::column(int i) {
  Matrix subMat(M, 0, -1, i, i);
  return subMat;
}

Matrix::operator vector<real>() const {
  vector<real> cat(numRows() * numCols());
  int k = 0;
  for (int i = 0; i < numRows(); i++) {
    for (int j = 0; j < numCols(); j++) {
      cat[k] = (*this)(i, j);
      k++;
    }
  }
  return cat;
}

Matrix Matrix::inverse() {
  if (size(1) != size(2)) MstUtils::error("inverse of non-square matrix requested", "Matrix::inverse()");
  int N = size(1);
  Matrix Mi = *this;
  Matrix& Mo = *this;
  real det;

  switch(N) {
    case 0:
      MstUtils::error("inverse of an empty matrix requested", "Matrix::inverse()");
      break;

    case 1:
      Mi(0,0) = 1/Mo(0,0);
      break;

    case 2:
      det = (Mo(0,0)*Mo(1,1) - Mo(0,1)*Mo(1,0));
      if (fabs(det) < 10E-15) MstUtils::warn("matrix too close to singular, determinant near zero", "Matrix::inverse()");
      Mi(0,0) =  Mo(1,1)/det;
      Mi(0,1) = -Mo(0,1)/det;
      Mi(1,0) = -Mo(1,0)/det;
      Mi(1,1) =  Mo(0,0)/det;
      break;

    case 3:
      det = (Mo(0,0)*Mo(1,1)*Mo(2,2) - Mo(0,0)*Mo(1,2)*Mo(2,1) - Mo(0,1)*Mo(1,0)*Mo(2,2) + Mo(0,1)*Mo(1,2)*Mo(2,0) + Mo(0,2)*Mo(1,0)*Mo(2,1) - Mo(0,2)*Mo(1,1)*Mo(2,0));
      if (fabs(det) < 10E-15) MstUtils::warn("matrix too close to singular, determinant near zero", "Matrix::inverse()");
      Mi(0,0) =  (Mo(1,1)*Mo(2,2) - Mo(1,2)*Mo(2,1))/det;
      Mi(0,1) = -(Mo(0,1)*Mo(2,2) - Mo(0,2)*Mo(2,1))/det;
      Mi(0,2) =  (Mo(0,1)*Mo(1,2) - Mo(0,2)*Mo(1,1))/det;
      Mi(1,0) = -(Mo(1,0)*Mo(2,2) - Mo(1,2)*Mo(2,0))/det;
      Mi(1,1) =  (Mo(0,0)*Mo(2,2) - Mo(0,2)*Mo(2,0))/det;
      Mi(1,2) = -(Mo(0,0)*Mo(1,2) - Mo(0,2)*Mo(1,0))/det;
      Mi(2,0) =  (Mo(1,0)*Mo(2,1) - Mo(1,1)*Mo(2,0))/det;
      Mi(2,1) = -(Mo(0,0)*Mo(2,1) - Mo(0,1)*Mo(2,0))/det;
      Mi(2,2) =  (Mo(0,0)*Mo(1,1) - Mo(0,1)*Mo(1,0))/det;
      break;

    case 4:
      det = (Mo(0,0)*Mo(1,1)*Mo(2,2)*Mo(3,3) - Mo(0,0)*Mo(1,1)*Mo(2,3)*Mo(3,2) - Mo(0,0)*Mo(1,2)*Mo(2,1)*Mo(3,3) + Mo(0,0)*Mo(1,2)*Mo(2,3)*Mo(3,1) + Mo(0,0)*Mo(1,3)*Mo(2,1)*Mo(3,2) - Mo(0,0)*Mo(1,3)*Mo(2,2)*Mo(3,1) - Mo(0,1)*Mo(1,0)*Mo(2,2)*Mo(3,3) + Mo(0,1)*Mo(1,0)*Mo(2,3)*Mo(3,2) + Mo(0,1)*Mo(1,2)*Mo(2,0)*Mo(3,3) - Mo(0,1)*Mo(1,2)*Mo(2,3)*Mo(3,0) - Mo(0,1)*Mo(1,3)*Mo(2,0)*Mo(3,2) + Mo(0,1)*Mo(1,3)*Mo(2,2)*Mo(3,0) + Mo(0,2)*Mo(1,0)*Mo(2,1)*Mo(3,3) - Mo(0,2)*Mo(1,0)*Mo(2,3)*Mo(3,1) - Mo(0,2)*Mo(1,1)*Mo(2,0)*Mo(3,3) + Mo(0,2)*Mo(1,1)*Mo(2,3)*Mo(3,0) + Mo(0,2)*Mo(1,3)*Mo(2,0)*Mo(3,1) - Mo(0,2)*Mo(1,3)*Mo(2,1)*Mo(3,0) - Mo(0,3)*Mo(1,0)*Mo(2,1)*Mo(3,2) + Mo(0,3)*Mo(1,0)*Mo(2,2)*Mo(3,1) + Mo(0,3)*Mo(1,1)*Mo(2,0)*Mo(3,2) - Mo(0,3)*Mo(1,1)*Mo(2,2)*Mo(3,0) - Mo(0,3)*Mo(1,2)*Mo(2,0)*Mo(3,1) + Mo(0,3)*Mo(1,2)*Mo(2,1)*Mo(3,0));
      if (fabs(det) < 10E-15) MstUtils::warn("matrix too close to singular, determinant near zero", "Matrix::inverse()");
      Mi(0,0) =  (Mo(1,1)*Mo(2,2)*Mo(3,3) - Mo(1,1)*Mo(2,3)*Mo(3,2) - Mo(1,2)*Mo(2,1)*Mo(3,3) + Mo(1,2)*Mo(2,3)*Mo(3,1) + Mo(1,3)*Mo(2,1)*Mo(3,2) - Mo(1,3)*Mo(2,2)*Mo(3,1))/det;
      Mi(0,1) = -(Mo(0,1)*Mo(2,2)*Mo(3,3) - Mo(0,1)*Mo(2,3)*Mo(3,2) - Mo(0,2)*Mo(2,1)*Mo(3,3) + Mo(0,2)*Mo(2,3)*Mo(3,1) + Mo(0,3)*Mo(2,1)*Mo(3,2) - Mo(0,3)*Mo(2,2)*Mo(3,1))/det;
      Mi(0,2) =  (Mo(0,1)*Mo(1,2)*Mo(3,3) - Mo(0,1)*Mo(1,3)*Mo(3,2) - Mo(0,2)*Mo(1,1)*Mo(3,3) + Mo(0,2)*Mo(1,3)*Mo(3,1) + Mo(0,3)*Mo(1,1)*Mo(3,2) - Mo(0,3)*Mo(1,2)*Mo(3,1))/det;
      Mi(0,3) = -(Mo(0,1)*Mo(1,2)*Mo(2,3) - Mo(0,1)*Mo(1,3)*Mo(2,2) - Mo(0,2)*Mo(1,1)*Mo(2,3) + Mo(0,2)*Mo(1,3)*Mo(2,1) + Mo(0,3)*Mo(1,1)*Mo(2,2) - Mo(0,3)*Mo(1,2)*Mo(2,1))/det;
      Mi(1,0) = -(Mo(1,0)*Mo(2,2)*Mo(3,3) - Mo(1,0)*Mo(2,3)*Mo(3,2) - Mo(1,2)*Mo(2,0)*Mo(3,3) + Mo(1,2)*Mo(2,3)*Mo(3,0) + Mo(1,3)*Mo(2,0)*Mo(3,2) - Mo(1,3)*Mo(2,2)*Mo(3,0))/det;
      Mi(1,1) =  (Mo(0,0)*Mo(2,2)*Mo(3,3) - Mo(0,0)*Mo(2,3)*Mo(3,2) - Mo(0,2)*Mo(2,0)*Mo(3,3) + Mo(0,2)*Mo(2,3)*Mo(3,0) + Mo(0,3)*Mo(2,0)*Mo(3,2) - Mo(0,3)*Mo(2,2)*Mo(3,0))/det;
      Mi(1,2) = -(Mo(0,0)*Mo(1,2)*Mo(3,3) - Mo(0,0)*Mo(1,3)*Mo(3,2) - Mo(0,2)*Mo(1,0)*Mo(3,3) + Mo(0,2)*Mo(1,3)*Mo(3,0) + Mo(0,3)*Mo(1,0)*Mo(3,2) - Mo(0,3)*Mo(1,2)*Mo(3,0))/det;
      Mi(1,3) =  (Mo(0,0)*Mo(1,2)*Mo(2,3) - Mo(0,0)*Mo(1,3)*Mo(2,2) - Mo(0,2)*Mo(1,0)*Mo(2,3) + Mo(0,2)*Mo(1,3)*Mo(2,0) + Mo(0,3)*Mo(1,0)*Mo(2,2) - Mo(0,3)*Mo(1,2)*Mo(2,0))/det;
      Mi(2,0) =  (Mo(1,0)*Mo(2,1)*Mo(3,3) - Mo(1,0)*Mo(2,3)*Mo(3,1) - Mo(1,1)*Mo(2,0)*Mo(3,3) + Mo(1,1)*Mo(2,3)*Mo(3,0) + Mo(1,3)*Mo(2,0)*Mo(3,1) - Mo(1,3)*Mo(2,1)*Mo(3,0))/det;
      Mi(2,1) = -(Mo(0,0)*Mo(2,1)*Mo(3,3) - Mo(0,0)*Mo(2,3)*Mo(3,1) - Mo(0,1)*Mo(2,0)*Mo(3,3) + Mo(0,1)*Mo(2,3)*Mo(3,0) + Mo(0,3)*Mo(2,0)*Mo(3,1) - Mo(0,3)*Mo(2,1)*Mo(3,0))/det;
      Mi(2,2) =  (Mo(0,0)*Mo(1,1)*Mo(3,3) - Mo(0,0)*Mo(1,3)*Mo(3,1) - Mo(0,1)*Mo(1,0)*Mo(3,3) + Mo(0,1)*Mo(1,3)*Mo(3,0) + Mo(0,3)*Mo(1,0)*Mo(3,1) - Mo(0,3)*Mo(1,1)*Mo(3,0))/det;
      Mi(2,3) = -(Mo(0,0)*Mo(1,1)*Mo(2,3) - Mo(0,0)*Mo(1,3)*Mo(2,1) - Mo(0,1)*Mo(1,0)*Mo(2,3) + Mo(0,1)*Mo(1,3)*Mo(2,0) + Mo(0,3)*Mo(1,0)*Mo(2,1) - Mo(0,3)*Mo(1,1)*Mo(2,0))/det;
      Mi(3,0) = -(Mo(1,0)*Mo(2,1)*Mo(3,2) - Mo(1,0)*Mo(2,2)*Mo(3,1) - Mo(1,1)*Mo(2,0)*Mo(3,2) + Mo(1,1)*Mo(2,2)*Mo(3,0) + Mo(1,2)*Mo(2,0)*Mo(3,1) - Mo(1,2)*Mo(2,1)*Mo(3,0))/det;
      Mi(3,1) =  (Mo(0,0)*Mo(2,1)*Mo(3,2) - Mo(0,0)*Mo(2,2)*Mo(3,1) - Mo(0,1)*Mo(2,0)*Mo(3,2) + Mo(0,1)*Mo(2,2)*Mo(3,0) + Mo(0,2)*Mo(2,0)*Mo(3,1) - Mo(0,2)*Mo(2,1)*Mo(3,0))/det;
      Mi(3,2) = -(Mo(0,0)*Mo(1,1)*Mo(3,2) - Mo(0,0)*Mo(1,2)*Mo(3,1) - Mo(0,1)*Mo(1,0)*Mo(3,2) + Mo(0,1)*Mo(1,2)*Mo(3,0) + Mo(0,2)*Mo(1,0)*Mo(3,1) - Mo(0,2)*Mo(1,1)*Mo(3,0))/det;
      Mi(3,3) =  (Mo(0,0)*Mo(1,1)*Mo(2,2) - Mo(0,0)*Mo(1,2)*Mo(2,1) - Mo(0,1)*Mo(1,0)*Mo(2,2) + Mo(0,1)*Mo(1,2)*Mo(2,0) + Mo(0,2)*Mo(1,0)*Mo(2,1) - Mo(0,2)*Mo(1,1)*Mo(2,0))/det;
      break;

    default:
      MstUtils::error("this is a basic matrix class that does not know how to invert matrices of dimension " + MstUtils::toString(N), "Matrix::inverse");
      break;
  }

  return Mi;
}

Matrix Matrix::transpose() {
  Matrix R(numCols(), numRows());
  for (int i = 0; i < M.size(); i++) {
    for (int j = 0; j < M[i].size(); j++) R(j,i) = *(M[i][j]);
  }
  return R;
}

Matrix Matrix::sum(int dim, bool norm) {
  switch(dim) {
    case 1: {
      Matrix S(1, numCols(), 0);
      for (int i = 0; i < numCols(); i++) {
        for (int j = 0; j < numRows(); j++) S(0, i) += (*this)(j, i);
        if (norm) S(0, i) /= numRows();
      }
      return S;
    }

    case 2: {
      Matrix S(numRows(), 1, 0);
      for (int i = 0; i < numRows(); i++) {
        for (int j = 0; j < numCols(); j++) S(i, 0) += (*this)(i, j);
        if (norm) S(0, i) /= numCols();
      }
      return S;
    }

    default:
      MstUtils::error("out of range dimension " + MstUtils::toString(dim));
      return Matrix(0, 0); // to make the compiler happy
  }
}

real Matrix::norm() const {
  real n = 0;
  for (int i = 0; i < M.size(); i++) {
    for (int j = 0; j < M[i].size(); j++) n += (*(M[i][j])) * (*(M[i][j]));
  }
  return sqrt(n);
}

real Matrix::norm2() const {
  real n = 0;
  for (int i = 0; i < M.size(); i++) {
    for (int j = 0; j < M[i].size(); j++) n += (*(M[i][j])) * (*(M[i][j]));
  }
  return n;
}
