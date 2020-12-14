#include "Vector.h"

const Vector Zv(0.0, 0.0, 0.0);
const Vector Iv(1.0, 0.0, 0.0);
const Vector Jv(0.0, 1.0, 0.0);
const Vector Kv(0.0, 0.0, 1.0);

//
//     Class Vector
//

// constructor
Vector::Vector() : VectorN(3)
{
  // dummy
}

Vector::Vector(double m1, double m2, double m3) : VectorN(3)
{ 
  v[0] = m1; v[1] = m2; v[2] = m3;
}

Vector::Vector(double *a) : VectorN(a, 3)
{
  //  for (int i=0; i<3; ++i) v[i] = a[i];
}

// copy constructor
Vector::Vector(const Vector& N) : VectorN(N.v, 3)
{ 
  //  for (int i=0; i<3; ++i) v[i] = N.v[i];
}

// destructor
Vector::~Vector()
{
  // dummy
}

// assingment operator
Vector& Vector::operator=(const Vector& N)
{
  if (this != &N) 
    for (int i=0; i<3; ++i) v[i] = N.v[i];

  return *this;
}

// ベクトルの要素
//double& Vector::operator[](int x)
//{ 
//  return v[x];
//}

// 出力
int Vector::Print(FILE *fp, char *str) const 
{ 
  fprintf(fp, "%s(%13.7g,%13.7g,%13.7g )\n", str, v[0],v[1],v[2]);

  return true;
}

int Vector::LongPrint(FILE *fp, char *str) const 
{ 
  fprintf(fp, "%s(%22.16g,%22.16g,%22.16g )\n", str, v[0],v[1],v[2]);

  return true;
}

// 入力 
int Vector::Scan(FILE *fp)
{
  int    i;

  for(i=0; i<3; ++i)
    if (ScanForScan(fp) == DIGIT) fscanf(fp, "%lf", &v[i]);
    else {
      fprintf(stderr, "Scan error, in Vector.\n");
      return false;
    }
  
  return true;
}

Vector& Vector::operator+=(const Vector& N)
{
  for (int i=0; i<3; ++i) v[i] += N.v[i];

  return *this;
}
  
Vector& Vector::operator-=(const Vector& N)
{
  for (int i=0; i<3; ++i) v[i] -= N.v[i];

  return *this;
}
  
Vector& Vector::operator*=(double x)
{
  for (int i=0; i<3; ++i) v[i] *= x;

  return *this;
}
  
Vector& Vector::operator/=(double x)
{
  if (fabs(x) < eps)
    fprintf(stderr, "Division by zero, Vector/=.\n");
  else
    for (int i=0; i<3; ++i) v[i] /= x;

  return *this;
}

// 次元を1次元増やす
VectorN operator+(const Vector& N, double x)
{
  int i;
  VectorN ret(4);

  for (i=0; i<3; ++i) ret[i] = N.v[i];

  ret[i] = x;

  return ret;
}

VectorN operator+(double x, const Vector& N)
{
  VectorN ret(4);

  ret[0] = x;
  for (int i=1; i<4; ++i) ret[i] = N.v[i-1];

  return ret;
}

// ベクトルの単項マイナス
Vector operator-(const Vector& N)
{
  return Vector(-N.v[0], -N.v[1], -N.v[2]);
}

// ベクトル同士の差
Vector operator-(const Vector& N1, const Vector& N2)
{ 
  return Vector(N1.v[0] - N2.v[0], N1.v[1] - N2.v[1], N1.v[2] - N2.v[2]);
}

// ベクトル同士の和
Vector operator+(const Vector& N1, const Vector& N2)
{ 
  return Vector(N1.v[0] + N2.v[0], N1.v[1] + N2.v[1], N1.v[2] + N2.v[2]);
}

// ベクトルを定数倍
Vector operator*(double r, const Vector& N)
{
  return Vector(N.v[0] * r, N.v[1] * r, N.v[2] * r);
}

// ベクトルを定数倍
Vector operator*(const Vector& N, double r)
{ 
  return Vector(N.v[0] * r, N.v[1] * r, N.v[2] * r);
}

// ベクトルの外積
Vector operator%(const Vector& N1, const Vector& N2)
{
  return Vector(N1.v[1] * N2.v[2] - N1.v[2] * N2.v[1],
		N1.v[2] * N2.v[0] - N1.v[0] * N2.v[2],
		N1.v[0] * N2.v[1] - N1.v[1] * N2.v[0]);
}

// ベクトルの商
Vector operator/(const Vector& N, double r)
{
  if (fabs(r) < eps) {
    fprintf(stderr, "Error : division by zero, Vector/\n");
    return Zv;
  } 
  else return Vector(N.v[0] / r, N.v[1] / r, N.v[2] / r);
}

// ベクトル同士の比較
int operator==(const Vector& N1, const Vector& N2)
{ 
  if (fabs(Norm(N1 - N2)) < eps) return true;
  else return false;
}

// ベクトル同士の比較
int operator!=(const Vector& N1, const Vector& N2)
{ 
  if (N1 == N2) return false;
  else return true;
}

// 内積
double operator,(const Vector& N1, const Vector& N2)
{
  double x=0.0;

  for(int i=0; i<3; ++i) x += N1.v[i] * N2.v[i];

  return x;
}

// スカラ３重積
double ScalarTP(const Vector& N1, const Vector& N2, const Vector& N3)
{
  return ((N1 % N2), N3);
}

// ベクトルのノルム
double Norm(const Vector& N)
{
  return sqrt(sqr(N.v[0]) + sqr(N.v[1]) + sqr(N.v[2]));
}

// ベクトルの正規化
Vector Normalize(const Vector& N)
{ 
  double d;

  d = Norm(N);

  if (fabs(d) < eps) return Zv;
  else
    return Vector(N.v[0] / d, N.v[1] / d, N.v[2] / d);
}

// 透視作用素
Vector Proj(const Vector& N)
{
  if (fabs(N.v[2]) > eps)
    return Vector(N.v[0] / N.v[2], N.v[1] / N.v[2], 1.0);
  else {
    fprintf(stderr, "Can't projection.\n");
    return N;
  }
}
