#ifndef _MSTLINALG_H
#define _MSTLINALG_H

#include <vector>
#include "msttypes.h"

using namespace std;

namespace MST {

/* This elemetary matrix class is provided for completeness and is not really
 * by other objects here. This class is not good for serious matrix algebra. */
class Matrix {
  public:
    Matrix(int rows, int cols, real val = 0.0);
    Matrix(const vector<vector<real> >& _M) { M = _M; }
    Matrix(const Matrix& _M) { M = _M.M; }
    Matrix(const vector<real>& p, bool col = false); // by default makes row vectors

    int size(bool dim) const;
    int numRows() { return size(0); }
    int numCols() { return size(1); }
    real& operator()(int i, int j) { return M[i][j]; }

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
    vector<real> operator[](int i) const { return M[i]; }
    vector<real>& operator[](int i) { return M[i]; }

    Matrix inverse();
    Matrix transpose();

    real norm();   // Euclidean norm of the matrix
    real norm2();  // Euclidean norm squared

    friend ostream & operator<<(ostream &_os, const Matrix& _M) {
      for (int i = 0; i < _M.size(0); i++) {
        _os << endl;
        for (int j = 0; j < _M.size(1); j++) {
          _os << std::setprecision(6) << std::fixed << _M[i][j] << " ";
        }
      }
      return _os;
    }

  private:
    vector<vector<real> > M;
};

}

#endif
