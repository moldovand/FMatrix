#ifndef _statcomp_h_
#define _statcomp_h_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

class VectorN;
class Vector;
class MatrixMN;
class Matrix;

// class Motion;

extern const Vector      Zv;
extern const Vector      Iv;
extern const Vector      Jv;
extern const Vector      Kv;

extern const Matrix      Im;
extern const Matrix      Zm;
extern const Matrix      Pk;

//extern const Motion      II;

#include "VectorN.h"
#include "Vector.h"
#include "MatrixMN.h"
#include "Matrix.h"

#include "ImagePosition.h"
#include "SpacePosition.h"

// #include <Motion.h>

#endif  // _statcomp_h_
