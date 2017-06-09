#ifndef _MSTLINALG_H
#define _MSTLINALG_H

#include <vector>
#include "msttypes.h"

using namespace std;

namespace MST {

class Matrix {
  public:
    Matrix(int rows, int cols, real val = 0.0);
    Matrix(const vector<vector<real> >& _M);
    Matrix(const Matrix& _M);
    Matrix(const vector<real>& p, bool col = false); // by default makes row vectors
    ~Matrix();

    int size(int dim = 1) const;
    int length() { return max(size(1), size(2)); }
    int numRows() const { return size(1); }
    int numCols() const { return size(2); }
    real& operator()(int i, int j) { return *(M[i][j]); }
    real operator()(int i, int j) const { return *(M[i][j]); }

    Matrix& operator=(const Matrix& _M);
    Matrix& operator/=(const real& s);
    Matrix& operator*=(const real& s);
    const Matrix operator/(const real& s) const;
    const Matrix operator*(const real& s) const;
    Matrix& operator*=(const Matrix& P);
    const Matrix operator*(const Matrix& P) const;
    Matrix& operator+=(const Matrix& P);
    const Matrix operator+(const Matrix& P) const;
    Matrix& operator-=(const Matrix& P);
    const Matrix operator-(const Matrix& P) const;
    const Matrix operator-() const;
    // so that we can also do scalar * Matrix (this is a global operator, not a member operator)
    friend Matrix operator* (real s, const Matrix& M) { return M*s; }

    Matrix row(int i); // sub-matrix corresponding to the i-th row
    Matrix column(int i); // sub-matrix corresponding to the i-th column

    // when typecast as a vector<real>, appends all rows together
    operator vector<real>() const;

    Matrix inverse();
    Matrix transpose();
    Matrix sum(int dim = 1, bool norm = false);
    Matrix mean(int dim = 1) { return sum(dim, true); }
    real norm() const;   // Euclidean norm of the matrix
    real norm2() const;  // Euclidean norm squared

    friend ostream & operator<<(ostream &_os, const Matrix& _M) {
      for (int i = 0; i < _M.numRows(); i++) {
        for (int j = 0; j < _M.numCols(); j++) {
          _os << std::setprecision(6) << std::fixed << _M(i, j) << " ";
        }
        if (i != _M.numRows() - 1) _os << endl;
      }
      return _os;
    }

  protected:
    Matrix(const vector<vector<real*> >& _M, int rowBeg = 0, int rowEnd = -1, int colBeg = 0, int colEnd = -1);
    void setOwnFlag(bool _own) { own = _own; }
    bool getOwnFlag() { return own; }
    vector<vector<real*> > M;
    bool own;                  // do I own the data in the matrix, or is my matrix a sub-matrix from some other matrix?
};

class Vector : public Matrix {
  public:
    Vector(int numel, real val = 0.0, bool col = false) : Matrix(vector<real>(numel, val), col) {}
    Vector(const vector<real>& p, bool col = false) : Matrix(p, col) {}
    Vector(const Vector& V) : Matrix(V) {}

    bool isRow() const { return numRows() >= numCols(); }
    bool isColumn() const { return numRows() <= numCols(); }
    real& operator()(int i) { return isRow() ? *(M[0][i]) : *(M[i][0]); }
    real operator()(int i) const { return isRow() ? *(M[0][i]) : *(M[i][0]); }

};

}

#endif
