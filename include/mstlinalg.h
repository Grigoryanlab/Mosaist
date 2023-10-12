#ifndef _MSTLINALG_H
#define _MSTLINALG_H

#include <vector>
#include "msttypes.h"
#undef assert

using namespace std;

namespace MST {

class Matrix {
  public:
    Matrix(int rows, int cols, mstreal val = 0.0);
    Matrix(const vector<vector<mstreal> >& _M);
    Matrix(const Matrix& _M);
    Matrix(const vector<mstreal>& p, bool col = false); // by default makes row vectors
    ~Matrix();

    int size(int dim = 1) const;
    int length() const { return MstUtils::max(size(1), size(2)); }
    int numRows() const { return size(1); }
    int numCols() const { return size(2); }
    mstreal& operator()(int i, int j) { return *(M[i][j]); }
    mstreal operator()(int i, int j) const { return *(M[i][j]); }

    // single-subscript access (column-major order, as in Matlab)
    mstreal& operator()(int i) { int n = numRows(); return *(M[i % n][i / n]); }
    mstreal operator()(int i) const { int n = numRows(); return *(M[i % n][i / n]); }
    mstreal& operator[](int i) { int n = numRows(); return *(M[i % n][i / n]); }
    mstreal operator[](int i) const { int n = numRows(); return *(M[i % n][i / n]); }

    Matrix& operator=(const Matrix& _M);
    Matrix& operator/=(const mstreal& s);
    Matrix& operator*=(const mstreal& s);
    const Matrix operator/(const mstreal& s) const;
    const Matrix operator*(const mstreal& s) const;
    Matrix& operator*=(const Matrix& P);
    const Matrix operator*(const Matrix& P) const;
    Matrix& operator+=(const Matrix& P);
    const Matrix operator+(const Matrix& P) const;
    Matrix& operator-=(const Matrix& P);
    const Matrix operator-(const Matrix& P) const;
    const Matrix operator-() const;
    // so that we can also do scalar * Matrix (this is a global operator, not a member operator)
    friend Matrix operator* (mstreal s, const Matrix& M) { return M*s; }

    Matrix row(int i); // sub-matrix corresponding to the i-th row
    Matrix column(int i); // sub-matrix corresponding to the i-th column

    // when typecast as a vector<real>, appends all rows together
    operator vector<mstreal>() const;
    // TODO: typecast as a Vector (fail if neither dimension is 1)

    Matrix inverse();
    Matrix transpose();
    Matrix sum(int dim = -1, bool norm = false) const;
    Matrix mean(int dim = -1) const { return sum(dim, true); }
    mstreal norm() const;   // Euclidean norm of the matrix
    mstreal norm2() const;  // Euclidean norm squared
    mstreal max() const;
    mstreal min() const;
    Matrix abs() const;
    Matrix mult(const Matrix& other) const; // element-wise multiply
    Matrix div(const Matrix& other) const;  // element-wise divide

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
    Matrix(const vector<vector<mstreal*> >& _M, int rowBeg = 0, int rowEnd = -1, int colBeg = 0, int colEnd = -1);
    void setOwnFlag(bool _own) { own = _own; }
    bool getOwnFlag() { return own; }
    void clear();
    vector<vector<mstreal*> > M;
    bool own;                  // do I own the data in the matrix, or is my matrix a sub-matrix from some other matrix?
};

// A convenience class. A Vector is really just a Matrix, with one of the
// dimensions required to be unity.
class Vector : public Matrix {
  public:
    Vector(int numel = 0, mstreal val = 0.0, bool col = false) : Matrix(vector<mstreal>(numel, val), col) {}
    Vector(const vector<mstreal>& p, bool col = false) : Matrix(p, col) {}
    Vector(const Vector& V) : Matrix(V) {}
    Vector(const Matrix& _M) : Matrix(_M) { MstUtils::assertCond((_M.numRows() == 1) || (_M.numCols() == 1), "cannot construct a vector from a matrix with non-unitary dimensions", "Vector::Vector(const Matrix& _M)"); }
    int size() const { return length(); }
    mstreal dot(const Vector& v) const;
    Vector getUnit() const { return (*this)/this->norm(); }
};

}

/* If compiling with a proper linear algebra package, many more things are possible */
#ifdef ARMA

#include <armadillo>
#define armaVec arma::Col<MST::mstreal>
#define armaMat arma::Mat<MST::mstreal>

// collection of linear algebra functions
class MstLinAlg {
  public:
    static MST::Matrix getPrincipalAxes(const MST::AtomPointerVector& atoms);
};

#endif

#endif
