#ifndef _MSTTRANSFORMS_H
#define _MSTTRANSFORMS_H

#include "msttypes.h"
#include <math.h>

namespace MST {

class Frame {
  public:
    Frame(); // same as the laboratory frame
    Frame(const CartesianPoint& _O, const CartesianPoint& _X, const CartesianPoint& _Y, const CartesianPoint& _Z);
    Frame(const CartesianPoint& _X, const CartesianPoint& _Y, const CartesianPoint& _Z) : Frame(CartesianPoint(0, 0, 0), _X, _Y, _Z) {};
    Frame(mstreal _ox, mstreal _oy, mstreal _oz, mstreal _xx, mstreal _xy, mstreal _xz, mstreal _yx, mstreal _yy, mstreal _yz, mstreal _zx, mstreal _zy, mstreal _zz);
    Frame(const Frame& other);

    CartesianPoint getX() const { return CartesianPoint(X[0], X[1], X[2]); }
    CartesianPoint getY() const { return CartesianPoint(Y[0], Y[1], Y[2]); }
    CartesianPoint getZ() const { return CartesianPoint(Z[0], Z[1], Z[2]); }
    CartesianPoint getO() const { return CartesianPoint(O[0], O[1], O[2]); }
    mstreal getX(int i) const { return X[i]; }
    mstreal getY(int i) const { return Y[i]; }
    mstreal getZ(int i) const { return Z[i]; }
    mstreal getO(int i) const { return O[i]; }

    void setX(const CartesianPoint& _X);
    void setY(const CartesianPoint& _Y);
    void setZ(const CartesianPoint& _Z);
    void setO(const CartesianPoint& _O);

    friend ostream & operator<<(ostream &_os, Frame& _F) {
      _os << "O: " << _F.O[0] << " " << _F.O[1] << " " << _F.O[2] << endl;
      _os << "X: " << _F.X[0] << " " << _F.X[1] << " " << _F.X[2] << endl;
      _os << "Y: " << _F.Y[0] << " " << _F.Y[1] << " " << _F.Y[2] << endl;
      _os << "Z: " << _F.Z[0] << " " << _F.Z[1] << " " << _F.Z[2] << endl;
      return _os;
    }

  protected:
    void constructFrame(const CartesianPoint& _O, const CartesianPoint& _X, const CartesianPoint& _Y, const CartesianPoint& _Z);

  private:
    mstreal O[3];             // coordinates of the origin in the laboratory frame
    mstreal X[3], Y[3], Z[3]; // vectors corresponding to the X, Y, and Z axes in the laboratory frame
};

class Transform {
  public:
    enum fillOrder { byColumn = 1, byRow };
    Transform() { makeIdentity(); }               // defaults to the identity transform
    Transform(const Transform& other);
    Transform(CartesianPoint A, CartesianPoint B, CartesianPoint C, fillOrder order);
    Transform(CartesianPoint A, CartesianPoint B, CartesianPoint C, CartesianPoint D, fillOrder order);
    Transform(vector<mstreal> trans);
    Transform(mstreal* trans);
    Transform(vector<vector<mstreal> > rot);
    Transform(mstreal** rot);
    Transform(vector<vector<mstreal> > rot, vector<mstreal> trans);
    Transform(mstreal** rot, mstreal* trans);
    void makeIdentity();
    void fill(vector<mstreal>& A, vector<mstreal>& B, vector<mstreal>& C, vector<mstreal>& D, fillOrder order);
    void fill(vector<mstreal>& A, vector<mstreal>& B, vector<mstreal>& C, fillOrder order);

    mstreal& operator()(int i, int j);             // for access and setting
    mstreal operator()(int i, int j) const;        // for access only
    Transform& operator*=(const Transform& rhs);
    const Transform operator*(const Transform& rhs) const;
    const CartesianPoint operator*(const CartesianPoint& rhs) const;
    Transform inverse();                   // return the einverse transform
    Transform rotation();                  // extract the rotation component
    Transform translation();               // extract the translation component

    /* Computes Euler angles corresponding to the the rotation component of the
     * current transform. Assume the order of rotations is around X, then Y,
     * then Z. Code adopted from stackoverflow.org based on derivations in:
     * http://www.soi.city.ac.uk/~sbbh653/publications/euler.pdf */
    void eulerAngles(mstreal& x, mstreal& y, mstreal& z);


    CartesianPoint applyToCopy(CartesianPoint& p);
    void apply(CartesianPoint& p);
    void apply(Frame& f);
    void apply(Atom& a) { apply(&a); }
    void apply(Atom* a);
    void apply(mstreal& x, mstreal& y, mstreal& z);
    void apply(Residue& res) { apply(&res); }
    void apply(Residue* res);
    void apply(Chain& chain) { apply(&chain); }
    void apply(Chain* chain);
    void apply(Structure& S) { apply(&S); }
    void apply(Structure* S);
    void apply(const AtomPointerVector& vec);

    void write(ostream& _os) const; // write Tansform to a binary stream
    void read(istream& _is);  // read Tansform from a binary stream


    friend ostream & operator<<(ostream &_os, const Transform& _T) {
      for (int i = 0; i < 4; i++) {
        _os << endl;
        for (int j = 0; j < 4; j++) {
          _os << std::setprecision(6) << std::fixed << _T.M[i][j] << "\t";
        }
      }
      return _os;
    }

  private:
    mstreal M[4][4];         // transformation matrix in homogeneous coordinates
};


class TransformFactory {
  public:
    static Transform translate(mstreal x, mstreal y, mstreal z);
    static Transform translate(const CartesianPoint& p);

    // rotations based on angle in degrees
    static Transform rotateAroundX(mstreal angle);
    static Transform rotateAroundY(mstreal angle);
    static Transform rotateAroundZ(mstreal angle);

    // rotations based on cosine and sine values (sometimes, it is convenient not to compute the angle itself)
    static Transform rotateAroundX(mstreal s, mstreal c);
    static Transform rotateAroundY(mstreal s, mstreal c);
    static Transform rotateAroundZ(mstreal s, mstreal c);

    // rotate around an arbitrary axis passing through the origin
    static Transform rotateAroundAxis(mstreal u, mstreal v, mstreal w, mstreal a);
    static Transform rotateAroundAxis(CartesianPoint& p, mstreal a);

    // rotate around an arbitrary line passing through two given points
    static Transform rotateAroundLine(mstreal p1, mstreal p2, mstreal p3, mstreal q1, mstreal q2, mstreal q3, mstreal a);
    static Transform rotateAroundLine(CartesianPoint& p, CartesianPoint& q, mstreal a);

    // transformation matrices for aligning arbitrary axes with laboratory frame axes
    static Transform alignVectorWithXAxis(mstreal u, mstreal v, mstreal w);
    static Transform alignVectorWithXAxis(CartesianPoint& p);

    static Transform alignVectorWithYAxis(mstreal u, mstreal v, mstreal w);
    static Transform alignVectorWithYAxis(CartesianPoint& p);

    static Transform alignVectorWithZAxis(mstreal u, mstreal v, mstreal w);
    static Transform alignVectorWithZAxis(CartesianPoint& p);

    // this transformation, when applied to points with coordinates in the from Frame
    // will produce coordinates corresponding the same points but in the to Frame
    static Transform switchFrames(const Frame& from, const Frame& to);

    static const mstreal degreesToRadians;
    static const mstreal radiansToDegrees;
};

/* This class is for fast calculation of the RMSD between two transformations
 * of a given structure. IMPORTANT: the transformations are interpreted as being
 * applied to the structure after it is centered at the origin. */
class TransformRMSD {
  public:
    TransformRMSD();
    TransformRMSD(Structure& S) { init(S); }
    TransformRMSD(AtomPointerVector& atoms) { init(atoms); }
    void init(const AtomPointerVector& atoms);
    void init(const Structure& S);

    /* IMPORTANT: the transformations are interpreted as being
     * applied to the structure after it is centered at the origin. */
    mstreal getRMSD(Transform& T1, Transform& T2);   // RMSD between two transforms
    mstreal getRMSD(Transform& T1);                  // RMSD between the given transform and the identity transform

  protected:

  private:
    mstreal C[3][3];  // co-variance matrix of structure in question (normalized by the number of atoms)
};

}

#endif
