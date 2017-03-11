#ifndef _MSTTRANSFORMS_H
#define _MSTTRANSFORMS_H

#include "msttypes.h"
#include <math.h>

namespace MST {

class Frame {
  public:
    Frame(); // same as the laboratory frame
    Frame(CartesianPoint& _O, CartesianPoint& _X, CartesianPoint& _Y, CartesianPoint& _Z);
    Frame(real _ox, real _oy, real _oz, real _xx, real _xy, real _xz, real _yx, real _yy, real _yz, real _zx, real _zy, real _zz);
    Frame(Frame& other);

    CartesianPoint getX() { return CartesianPoint(X[0], X[1], X[2]); }
    CartesianPoint getY() { return CartesianPoint(Y[0], Y[1], Y[2]); }
    CartesianPoint getZ() { return CartesianPoint(Z[0], Z[1], Z[2]); }
    CartesianPoint getO() { return CartesianPoint(O[0], O[1], O[2]); }
    real getX(int i) { return X[i]; }
    real getY(int i) { return Y[i]; }
    real getZ(int i) { return Z[i]; }
    real getO(int i) { return O[i]; }

    friend ostream & operator<<(ostream &_os, Frame& _F) {
      _os << "O: " << _F.O[0] << " " << _F.O[1] << " " << _F.O[2] << endl;
      _os << "X: " << _F.X[0] << " " << _F.X[1] << " " << _F.X[2] << endl;
      _os << "Y: " << _F.Y[0] << " " << _F.Y[1] << " " << _F.Y[2] << endl;
      _os << "Z: " << _F.Z[0] << " " << _F.Z[1] << " " << _F.Z[2] << endl;
      return _os;
    }

  private:
    real O[3];             // coordinates of the origin in the laboratory frame
    real X[3], Y[3], Z[3]; // vectors corresponding to the X, Y, and Z axes in the laboratory frame
};

class Transform {
  public:
    enum fillOrder { byColumn = 1, byRow };
    Transform() { makeIdentity(); }               // defaults to the identity transform
    Transform(const Transform& other);
    Transform(CartesianPoint A, CartesianPoint B, CartesianPoint C, fillOrder order);
    Transform(CartesianPoint A, CartesianPoint B, CartesianPoint C, CartesianPoint D, fillOrder order);
    Transform(vector<real> trans);
    Transform(vector<vector<real> > rot);
    Transform(vector<vector<real> > rot, vector<real> trans);
    void makeIdentity();
    void fill(vector<real>& A, vector<real>& B, vector<real>& C, vector<real>& D, fillOrder order);
    void fill(vector<real>& A, vector<real>& B, vector<real>& C, fillOrder order);

    real& operator()(int i, int j);             // for access and setting
    real operator()(int i, int j) const;        // for access only
    Transform& operator*=(const Transform& rhs);
    const Transform operator*(const Transform& rhs) const;
    const CartesianPoint operator*(const CartesianPoint& rhs) const;
    Transform inverse();                   // return the einverse transform
    Transform rotation();                  // extract the rotation component
    Transform translation();               // extract the translation component


    CartesianPoint applyToCopy(CartesianPoint& p);
    void apply(CartesianPoint& p);
    void apply(Atom& a) { apply(&a); }
    void apply(Atom* a);
    void apply(real& x, real& y, real& z);
    void apply(Residue& res) { apply(&res); }
    void apply(Residue* res);
    void apply(Chain& chain) { apply(&chain); }
    void apply(Chain* chain);
    void apply(Structure& S) { apply(&S); }
    void apply(Structure* S);

    friend ostream & operator<<(ostream &_os, Transform& _T) {
      for (int i = 0; i < 4; i++) {
        _os << endl;
        for (int j = 0; j < 4; j++) {
          _os << std::setprecision(6) << std::fixed << _T.M[i][j] << "\t";
        }
      }
      return _os;
    }
  private:
    real M[4][4];         // transformation matrix in homogeneous coordinates
};


class TransformFactory {
  public:
    static Transform translate(real x, real y, real z);
    static Transform translate(CartesianPoint& p);

    // rotations based on angle in degrees
    static Transform rotateAroundX(real angle);
    static Transform rotateAroundY(real angle);
    static Transform rotateAroundZ(real angle);

    // rotations based on cosine and sine values (sometimes, it is convenient not to compute the angle itself)
    static Transform rotateAroundX(real s, real c);
    static Transform rotateAroundY(real s, real c);
    static Transform rotateAroundZ(real s, real c);

    // rotate around an arbitrary axis passing through the origin
    static Transform rotateAroundAxis(real u, real v, real w, real a);
    static Transform rotateAroundAxis(CartesianPoint& p, real a);

    // rotate around an arbitrary line passing through two given points
    static Transform rotateAroundLine(real p1, real p2, real p3, real q1, real q2, real q3, real a);
    static Transform rotateAroundLine(CartesianPoint& p, CartesianPoint& q, real a);

    // transformation matrices for aligning arbitrary axes with laboratory frame axes
    static Transform alignVectorWithXAxis(real u, real v, real w);
    static Transform alignVectorWithXAxis(CartesianPoint& p);

    static Transform alignVectorWithYAxis(real u, real v, real w);
    static Transform alignVectorWithYAxis(CartesianPoint& p);

    static Transform alignVectorWithZAxis(real u, real v, real w);
    static Transform alignVectorWithZAxis(CartesianPoint& p);

    // this transformation, when applied to points with coordinates in the from Frame
    // will produce coordinates of corresponding the same points but in the to Frame
    static Transform switchFrames(Frame& from, Frame& to);

    static const real degreesToRadians;
    static const real radiansToDegrees;
};

/* This class is for fast calculation of the RMSD between two transformations
 * of a given structure. */
class TransformRMSD {
  public:
    TransformRMSD();
    TransformRMSD(Structure& S) { init(S); }
    TransformRMSD(AtomPointerVector& atoms) { init(atoms); }
    void init(const AtomPointerVector& atoms);
    void init(const Structure& S);

    real getRMSD(Transform& T1, Transform& T2);   // RMSD between two transforms
    real getRMSD(Transform& T1);                  // RMSD between the given transform and the identity transform

  protected:

  private:
    real C[3][3];  // co-variance matrix of structure in question
    int N;         // number of atoms in the structure in question
};

}

#endif
