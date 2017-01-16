#ifndef _MSTTRANSFORMS_H
#define _MSTTRANSFORMS_H

#include "msttypes.h"

namespace MST {

class Frame3D { // TODO
  public:
    Frame3D(); // same as the laboratory frame
    Frame3D(CartesianPoint& _O, CartesianPoint& _X, CartesianPoint& _Y, CartesianPoint& _Z);
    Frame3D(real _ox, real _oy, real _oz, real _xx, real _xy, real _xz, real _yx, real _yy, real _yz, real _zx, real _zy, real _zz);
    Frame3D(Frame3D& other);

  private:
    real O[3];             // coordinates of the origin in the laboratory frame
    real X[3], Y[3], Z[3]; // vectors corresponding to the X, Y, and Z axes in the laboratory frame
};

class Transform3D { // TODO
  public:
    Transform3D();                   // defaults to the identity transform
    Transform3D(Transform3D& other);

    real& operator()(int i, int j);  // for access and setting
    operator*                        // multiplication of transforms
    CartesianPoint appyToCopy(CartesianPoint& p);
    void apply(CartesianPoint& p);
    void apply(Atom& a);
    void apply(Atom* a);
    void apply(real& x, real& y, real& z);

  private:
    real M[4][4];         // transformation matrix in homogeneous coordinates
};


class Transform3DFactory { // TODO
  public:
    Transform3D rotateAroundX(real angle);
    Transform3D rotateAroundY(real angle);
    Transform3D rotateAroundZ(real angle);

    Transform3D rotateAroundOrigAxis(real angle, real u, real v, real w);
    Transform3D rotateAroundOrigAxis(real angle, CartesianPoint& p);

    Transform3D translate(real x, real y, real z);
    Transform3D translate(CartesianPoint& p);

    Transform3D alignVectorWithXAxis(real x, real y, real z);
    Transform3D alignVectorWithXAxis(CartesianPoint& p);

    Transform3D alignVectorWithYAxis(real x, real y, real z);
    Transform3D alignVectorWithYAxis(CartesianPoint& p);

    Transform3D alignVectorWithZAxis(real x, real y, real z);
    Transform3D alignVectorWithZAxis(CartesianPoint& p);

    Transform3D switchFrames(Frame3D& from, Frame3D& to);

}

#endif
