#ifndef _MSTTRANSFORMS_H
#define _MSTTRANSFORMS_H

#include "msttypes.h"

namespace MST {

class Frame { // TODO
  public:
    Frame(); // same as the laboratory frame
    Frame(CartesianPoint& _O, CartesianPoint& _X, CartesianPoint& _Y, CartesianPoint& _Z);
    Frame(real _ox, real _oy, real _oz, real _xx, real _xy, real _xz, real _yx, real _yy, real _yz, real _zx, real _zy, real _zz);
    Frame(Frame& other);

  private:
    real O[3];             // coordinates of the origin in the laboratory frame
    real X[3], Y[3], Z[3]; // vectors corresponding to the X, Y, and Z axes in the laboratory frame
};

class Transform { // TODO
  public:
    Transform();                   // defaults to the identity transform
    Transform(Transform& other);

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


class TransformFactory { // TODO
  public:
    Transform rotateAroundX(real angle);
    Transform rotateAroundY(real angle);
    Transform rotateAroundZ(real angle);

    Transform rotateAroundOrigAxis(real angle, real u, real v, real w);
    Transform rotateAroundOrigAxis(real angle, CartesianPoint& p);

    Transform translate(real x, real y, real z);
    Transform translate(CartesianPoint& p);

    Transform alignVectorWithXAxis(real x, real y, real z);
    Transform alignVectorWithXAxis(CartesianPoint& p);

    Transform alignVectorWithYAxis(real x, real y, real z);
    Transform alignVectorWithYAxis(CartesianPoint& p);

    Transform alignVectorWithZAxis(real x, real y, real z);
    Transform alignVectorWithZAxis(CartesianPoint& p);

    Transform switchFrames(Frame& from, Frame& to);

}

#endif
