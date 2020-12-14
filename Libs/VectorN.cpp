#include "VectorN.h"

double sqr(double x)
{
  return pow(x, 2.0);
}

int sgn(double x)
{
  if      (x > eps)  return 1;
  else if (x < -eps) return -1;
  else               return 0;
}

double cbrt(double x)
{
  if (x == 0.0) return 0.0;
  else
    if (x < 0.0) return -exp(log(-x) / 3);
    else         return  exp(log(x) / 3);
}

//
//     Class VectorN (N-Dimension Vector)
//

// constructor
VectorN::VectorN(int d)
{
  dim = d;
  v = new double[dim];
  for (int i=0; i<dim; ++i) v[i] = 0.0;
}

VectorN::VectorN(double *a, int d)
{
  dim = d;
  v = new double[dim];
  for (int i=0; i<dim; ++i) v[i] = a[i];
}

// copy constructor
VectorN::VectorN(const VectorN& N)
{
  dim = N.dim;
  v = new double[dim];
  for (int i=0; i<dim; ++i) v[i] = N.v[i];
}

// destructor
VectorN::~VectorN()
{
  delete [] v;
}

// assingment operator
VectorN& VectorN::operator=(const VectorN& N)
{
  if (this != &N)
    if (dim == N.dim) for (int i=0; i<dim; ++i) v[i] = N.v[i];
    else fprintf(stderr, "Dimension error, VectorN =.\n");
  
  return *this;
}

// ベクトルの要素
//double& VectorN::operator[](int x)
//{ 
//  return v[x];
//}

// 出力
int VectorN::Print(FILE *fp, char *str) const
{
  int i;
  fprintf(fp, "%s(", str);
  for (i=0; i<dim-1; ++i) fprintf(fp, "%13.7g,",v[i]);
  fprintf(fp, "%13.7g )\n", v[i]);

  return true;
}

int VectorN::LongPrint(FILE *fp, char *str) const
{ 
  int i;
  fprintf(fp, "%s(", str);
  for (i=0; i<dim-1; ++i) fprintf(fp, "%22.16g,",v[i]);
  fprintf(fp, "%22.16g )\n", v[i]);

  return true;
}

// 入力 
int VectorN::Scan(FILE *fp)
{
  int    i;

  for(i=0; i<dim; ++i)
    if (ScanForScan(fp) == DIGIT) fscanf(fp, "%lf", &v[i]);
    else {
      fprintf(stderr, "Scan error, in VectorN.\n");
      return false;
    }
  
  return true;
}

// 初期化
int VectorN::Init(int d)
{
  if (d != 0) {
    delete [] v;

    dim = d;
    v = new double[dim];
  }
  
  for (int i=0; i<dim; ++i) v[i] = 0.0;

  return true;
}
    
VectorN& VectorN::operator+=(const VectorN& N)
{
  if (dim == N.dim)
    for (int i=0; i<dim; ++i) v[i] += N.v[i];

  return *this;
}
  
VectorN& VectorN::operator-=(const VectorN& N)
{
  if (dim == N.dim)
    for (int i=0; i<dim; ++i) v[i] -= N.v[i];

  return *this;
}
  
VectorN& VectorN::operator*=(double x)
{
  for (int i=0; i<dim; ++i) v[i] *= x;

  return *this;
}
  
VectorN& VectorN::operator/=(double x)
{
  if (fabs(x) < eps)
    fprintf(stderr, "Division by zero, VectorN/=.\n");
  else
    for (int i=0; i<dim; ++i) v[i] /= x;

  return *this;
}

// 次元を1次元増やす
VectorN operator+(const VectorN& N, double x)
{
  int i;
  VectorN ret(N.dim+1);

  for (i=0; i<ret.dim-1; ++i) ret.v[i] = N.v[i];

  ret.v[i] = x;

  return ret;
}
     
VectorN operator+(double x, const VectorN& N)
{
  VectorN ret(N.dim+1);

  ret.v[0] = x;
  for (int i=1; i<ret.dim; ++i) ret.v[i] = N.v[i-1];

  return ret;
}  

// ベクトルの単項マイナス
VectorN operator-(const VectorN& N)
{
  VectorN ret(N.dim);

  for (int i=0; i<ret.dim; ++i) ret.v[i] = -N.v[i];
  
  return ret;
}

// ベクトル同士の差
VectorN operator-(const VectorN& N1, const VectorN& N2)
{
  if (N1.dim == N2.dim) {
    VectorN ret(N1);
    for (int i=0; i<N1.dim; ++i) ret.v[i] -= N2.v[i];

    return ret;
  }
  else {
    fprintf(stderr, "Dimension error, VectorN -.\n");
    return N1;
  }
}

// ベクトル同士の和
VectorN operator+(const VectorN& N1, const VectorN& N2)
{ 
  if (N1.dim == N2.dim) {
    VectorN ret(N1);
    for (int i=0; i<N1.dim; ++i) ret.v[i] += N2.v[i];
    
    return ret;
  }
  else {
    fprintf(stderr, "Dimension error, VectorN +.\n");
    return N1;
  }
}

// ベクトルを定数倍
VectorN operator*(double r, const VectorN& N)
{
  VectorN ret(N);

  for (int i=0; i<ret.dim; ++i) ret.v[i] *= r;
  
  return ret;
}

// ベクトルを定数倍
VectorN operator*(const VectorN& N, double r)
{ 
  VectorN ret(N);

  for (int i=0; i<ret.dim; ++i) ret.v[i] *= r;
  
  return ret;
}

// ベクトルの外積
VectorN operator%(const VectorN& N1, const VectorN& N2)
{
  if ((N1.dim == N2.dim) && (N1.dim >= 3)) {
    int i;
    VectorN ret(N1.dim);
    
    for (i=0; i<ret.dim-2; ++i)
      ret.v[i] = N1.v[i+1] * N2.v[i+2] - N1.v[i+2] * N2.v[i+1];
    
    ret.v[i] = N1.v[i+1] * N2.v[0] - N1.v[0] * N2.v[i+1];
    ret.v[i+1] = N1.v[0] * N2.v[1] - N1.v[1] * N2.v[0];

    return ret;
  }
  else {
    fprintf(stderr, "Can't calcurate Outer Product.\n");
    return N1;
  }
}

// ベクトルの商
VectorN operator/(const VectorN& N, double r)
{
  if (fabs(r) < eps) {
    fprintf(stderr, "Error : division by zero, VectorN.\n");
    return N;
  } 
  else {
    VectorN ret(N);
    for (int i=0; i<ret.dim; ++i) ret.v[i] /= r;

    return ret;
  }
}

// ベクトル同士の比較
int operator==(const VectorN& N1, const VectorN& N2)
{
  if ((N1.dim == N2.dim) && (fabs(Norm(N1 - N2)) < eps)) return true;
  else return false;
}

// ベクトル同士の比較
int operator!=(const VectorN& N1, const VectorN& N2)
{ 
  if (N1 == N2) return false;
  else return true;
}

// 内積
double operator,(const VectorN& N1, const VectorN& N2)
{
  if (N1.dim == N2.dim) {
    double x = 0.0;
    for(int i=0; i<N1.dim; ++i) x += N1.v[i] * N2.v[i];
    return x;
  }
  else {
    fprintf(stderr, "Dimension error, Inner Product.\n");

    return 0.0;
  }
}

// ベクトルのノルム
double Norm(const VectorN& N)
{
  double x = 0.0;

  for (int i=0; i<N.dim; ++i) x += sqr(N.v[i]);
  
  return sqrt(x);
}

// ベクトルの正規化
VectorN Normalize(const VectorN& N)
{ 
  double d;

  d = Norm(N);

  if (fabs(d) < eps) return N;
  else               return N / d;
}

// 透視作用素
VectorN Proj(const VectorN& N)
{
  if ((N.dim == 3) && (fabs(N.v[N.dim-1]) > eps)) {
    VectorN ret(N);
    
    for (int i=0; i<ret.dim-1; ++i) ret.v[i] /= N.v[N.dim-1];
    ret.v[N.dim-1] = 1.0;
    return ret;
  }
  else {
    fprintf(stderr, "Can't projection.\n");
    return N;
  }
}


