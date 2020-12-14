#include <stdlib.h>
#include <math.h>
#include "MatrixMN.h"
#include "alg.h"

// エディングトンのε
int EddingtonEpsilon(int i, int j, int k)
{
  if((i==0 && j==1 && k==2) || (i==1 && j==2 && k==0) ||
     (i==2 && j==0 && k==1))
    return 1;
  else if((i==0 && j==2 && k==1) || (i==1 && j==0 && k==2) ||
	  (i==2 && j==1 && k==0))
    return -1;
  else return 0;
}

// クロネッカーのデルタ
int KroneckerDelta(int i, int j)
{
  return((i==j)? 1: 0);
}

//
//     Class MatrixMN
//

// constructor
MatrixMN::MatrixMN(int di, int dj)
{
  int i;
  dim_i = di; dim_j = dj;

  m = new double*[dim_i];
  for (i=0; i<dim_i; ++i) m[i] = new double[dim_j];
  
  for (i=0; i<dim_i; ++i)
    for (int j=0; j<dim_j; ++j) m[i][j] = 0.0;
}

MatrixMN::MatrixMN(const VectorN& N1, const VectorN& N2, const VectorN& N3)
{
  if ((N1.Dim() == N2.Dim()) && (N2.Dim() == N3.Dim())) {
    int i;
    dim_i = N1.Dim(); dim_j = 3;

    m = new double*[dim_i];
    for (i=0; i<dim_i; ++i) m[i] = new double[dim_j];

    for (i=0; i<dim_i; ++i) {
      m[i][0] = N1[i]; m[i][1] = N2[i]; m[i][2] = N3[i];
    }
  }
  else {
    dim_i = 0; dim_j = 0;
  }
}

MatrixMN::MatrixMN(double* a, int di, int dj)
{
  int i;
  dim_i = di; dim_j = dj;

  m = new double*[dim_i];
  for (i=0; i<dim_i; ++i) m[i] = new double[dim_j];

  for (i=0; i<dim_i; ++i)
    for (int j=0; j<dim_j; ++j) m[i][j] = a[i * dim_j + j];
}

MatrixMN::MatrixMN(double** a, int di, int dj)
{
  int i;
  dim_i = di; dim_j = dj;

  m = new double*[dim_i];
  for (i=0; i<dim_i; ++i) m[i] = new double[dim_j];

  for (i=0; i<dim_i; ++i)
    for (int j=0; j<dim_j; ++j) m[i][j] = a[i][j];
}

MatrixMN::MatrixMN(const VectorN& N1, const VectorN& N2)
{
  int i;
  dim_i = N1.Dim(); dim_j = N2.Dim();

  m = new double*[dim_i];
  for (i=0; i<dim_i; ++i) m[i] = new double[dim_j];
  
  for (i=0; i<dim_i; ++i)
    for (int j=0; j<dim_j; ++j) m[i][j] = N1[i] * N2[j];
}
/*
MatrixMN::MatrixMN(const Matrix& M, const Vector& N1, const Vector& N2, double x)
{
  dim_i = 4; dim_j = 4;

  m = new double*[dim_i];
  for (int i=0; i<dim_i; ++i) m[i] = new double[dim_j];

  for (i=0; i<3; ++i)
    for (int j=0; j<3; ++j) m[i][j] = M.m[i][j];
  
  for (i=0; i<3; ++i) {
    m[i][3] = N1[i];
    m[3][i] = N2[i];
  }

  m[3][3] = x;
}
*/
// copy constructor
MatrixMN::MatrixMN(const MatrixMN& M)
{
  int i;
  dim_i = M.dim_i; dim_j = M.dim_j;

  m = new double*[dim_i];
  for (i=0; i<dim_i; ++i) m[i] = new double[dim_j];
  
  for(i=0; i<dim_i; ++i)
    for(int j=0; j<dim_j; ++j) m[i][j] = M.m[i][j];
}

// destructor
MatrixMN::~MatrixMN()
{
  for (int i=0; i<dim_i; ++i) delete [] m[i];
  delete [] m;
}

// assignment operator
MatrixMN& MatrixMN::operator=(const MatrixMN& M)
{
  if (this != &M)
    if ((dim_i == M.dim_i) && (dim_j == M.dim_j))
      for (int i=0; i<dim_i; ++i) 
	for (int j=0; j<dim_j; ++j) m[i][j] = M.m[i][j];
    else fprintf(stderr, "Dimension error. MatrixMN =.\n");
    
  return *this;
}

// i行の行ベクトル
double* MatrixMN::operator[](int i)
{
  if ((0 <= i) && (i < dim_i))
    return m[i];
  else{
    fprintf(stderr, "Domain error in MatrixMN[].\n");
    return NULL;
  }
}

const double* MatrixMN::operator[](int i) const
{
  if ((0 <= i) && (i < dim_i))
    return m[i];
  else{
    fprintf(stderr, "Domain error in MatrixMN[].\n");
    return NULL;
  }
}

// j列の列ベクトル
VectorN MatrixMN::operator()(int j) const 
{
  if ((0 <= j) && (j < dim_j)) {
    VectorN ret(dim_i);
    for (int i=0; i<dim_i; ++i) ret[i] = m[i][j];
  
    return ret;
  }
  else{
    fprintf(stderr, "Domain error in MatrixMN().\n");
    return VectorN(0);
  }
}

// MatrixMNを出力
int MatrixMN::Print(FILE *fp, char *str) const 
{ 
  char sp[20];
  int  i, len, cent = (dim_i - 1) / 2;
  
  len = strlen(str);
  for (i=0; i<len; ++i) sp[i] = ' ';
  sp[len] = '\0';
  
  for (i=0; i<dim_i; ++i) {
    fprintf(fp, "%s|", (i == cent) ? str : sp);
    for(int j=0; j<dim_j; ++j) fprintf(fp, "%13.7g ", m[i][j]);
    fprintf(fp, "|\n");
  }

  return true;
}

int MatrixMN::LongPrint(FILE *fp, char *str) const 
{ 
  char sp[20];
  int  i, len, cent = (dim_i - 1) / 2;
  
  len = strlen(str);
  for (i=0; i<len; ++i) sp[i] = ' ';
  sp[len] = '\0';
  
  for (i=0; i<dim_i; ++i) {
    fprintf(fp, "%s|", (i == cent) ? str : sp);
    for(int j=0; j<dim_j; ++j) fprintf(fp, "%22.16g ", m[i][j]);
    fprintf(fp, "|\n");
  }

  return true;
}

// 行列の入力
int MatrixMN::Scan(FILE *fp)
{
  for (int i=0; i<dim_i; ++i) 
    for (int j=0; j<dim_j; ++j)
      if (ScanForScan(fp) == DIGIT) fscanf(fp, "%lf", &m[i][j]);
      else {
	fprintf(stderr, "Scan error, in MatrixMN.\n");
	return false;
      }

  return true;
}

// 対称行列かどうか
int MatrixMN::Symmetric() const
{
  int ret = true;

  if (dim_i != dim_j) {
    fprintf(stderr, "MatrixMN isn't Symmetric.\n");
    ret = false;
  }
  else 
    for (int i=0; i<dim_i-1; ++i)
      for (int j=i+1; j<dim_j; ++j)
	if (fabs(m[i][j] - m[j][i]) > eps) {
	  ret = false;
	  break;
	}
  
  return ret;
}

// 初期化
int MatrixMN::Init(int i, int j)
{
  int k;
  
  if ((i != 0) && (j != 0)) {
    for (k=0; k<dim_i; ++k) delete [] m[k];
    delete [] m;

    dim_i = i; dim_j = j;
    m = new double*[dim_i];
    for (k=0; k<dim_i; ++k) m[k] = new double[dim_j];
  }

  for (k=0; k<dim_i; ++k)
    for (int l=0; l<dim_j; ++l) m[k][l] = 0.0;
  
  return true;
}

// 単位化
int MatrixMN::Unit()
{
  for (int i=0; i<dim_i; ++i)
    for (int j=0; j<dim_j; ++j)
      if (i == j) m[i][j] = 1.0;
      else        m[i][j] = 0.0;

  return true;
}

// 比較
int operator==(const MatrixMN& M1, const MatrixMN& M2)
{
  if ((M1.dim_i == M2.dim_i) && (M1.dim_j == M2.dim_j) &&
      (fabs(Norm(M1 - M2)) <= eps)) return true;
  else                              return false;
}

int operator!=(const MatrixMN& M1, const MatrixMN& M2)
{
  if (M1 == M2) return false;
  else          return true;
}

MatrixMN& MatrixMN::operator+=(const MatrixMN& M)
{
  if ((dim_i == M.dim_i) && (dim_j == M.dim_j))
    for (int i=0; i<dim_i; ++i)
      for (int j=0; j<dim_j; ++j) m[i][j] += M.m[i][j];
  
  return *this;
}

MatrixMN& MatrixMN::operator-=(const MatrixMN& M)
{
  if ((dim_i == M.dim_i) && (dim_j == M.dim_j))
    for (int i=0; i<dim_i; ++i)
      for (int j=0; j<dim_j; ++j) m[i][j] -= M.m[i][j];

  return *this;
}

MatrixMN& MatrixMN::operator*=(double x)
{
  for (int i=0; i<dim_i; ++i)
    for (int j=0; j<dim_j; ++j) m[i][j] *= x;

  return *this;
}

MatrixMN& MatrixMN::operator/=(double x)
{
  if (fabs(x) < eps)
    fprintf(stderr, "Division by ZERO, MatrixMN/=.\n");
  else
    for (int i=0; i<dim_i; ++i)
      for (int j=0; j<dim_j; ++j) m[i][j] /= x;

  return *this;
}

// 次元を1次元増やす
MatrixMN operator+(const MatrixMN& M, double x)
{
  MatrixMN ret(M.dim_i + 1, M.dim_j + 1);
  
  for (int i=0; i<ret.dim_i; ++i)
    for (int j=0; j<ret.dim_j; ++j)
      if ((i == ret.dim_i-1) || (j == ret.dim_j-1)) ret.m[i][j] = 0.0;
      else    	                                    ret.m[i][j] = M.m[i][j];

  ret.m[ret.dim_i - 1][ret.dim_j - 1] = x;

  return ret;
}
  
MatrixMN operator+(double x, const MatrixMN& M)
{
  MatrixMN ret(M.dim_i + 1, M.dim_j + 1);
  
  for (int i=0; i<ret.dim_i; ++i)
    for (int j=0; j<ret.dim_j; ++j)
      if ((i == 0) || (j == 0))	ret.m[i][j] = 0.0;
      else     	                ret.m[i][j] = M.m[i - 1][j - 1];

  ret.m[0][0] = x;

  return ret;
}
  
// 行列とベクトルの積
VectorN operator*(const MatrixMN& M, const VectorN& N)
{
  if (M.dim_j == N.Dim()) {
    VectorN ret(M.dim_i);
  
    for (int i=0; i<ret.Dim(); ++i)
      for (int j=0; j<M.dim_j; ++j) ret[i] += M.m[i][j] * N[j];
 
    return ret;
  }
  else {
    fprintf(stderr, "Can't product, MatrixMN.\n");
    return VectorN(0);
  }
}

// 行列と行列の積
MatrixMN operator*(const MatrixMN& M1, const MatrixMN& M2)
{
  if (M1.dim_j == M2.dim_i) {
    MatrixMN ret(M1.dim_i, M2.dim_j);

    for (int i=0; i<ret.dim_i; ++i)
      for (int j=0; j<ret.dim_j; ++j)
	for (int k=0; k<M1.dim_j; ++k)
	  ret.m[i][j] += M1.m[i][k] * M2.m[k][j];

    return ret;
  }
  else {
    fprintf(stderr, "Can't product, MatrixMN *.\n");
    return M1;
  }
}

// 行列とスカラーの積
MatrixMN operator*(double r, const MatrixMN& M)
{
  MatrixMN ret(M);
  
  for (int i=0; i<ret.dim_i; ++i)
    for (int j=0; j<ret.dim_j; ++j) ret.m[i][j] *= r;
 
  return ret;
}

MatrixMN operator*(const MatrixMN& M, double r)
{
  return r * M;
}

// ベクトルと行列の外積
MatrixMN operator%(const VectorN& N, const MatrixMN& M)
{
  if (N.Dim() == M.dim_i) {
    MatrixMN ret(M.dim_i, M.dim_j);
    VectorN  tmp(M.dim_i);

    for (int j=0; j<M.dim_j; ++j) {
      tmp = N % M(j);
      for (int i=0; i<M.dim_i; ++i) ret.m[i][j] = tmp[i];
    }

    return ret;
  }
  else {
    fprintf(stderr, "Can't product, MatrixMN.\n");
    return MatrixMN(0, 0);
  }
}

// 行列のスカラ分の1
MatrixMN operator/(const MatrixMN& M, double r)
{
  if (fabs(r) > eps) return (1.0 / r) * M;
  else {
    fprintf(stderr, "Error : division by zero in MatrixMN/.\n");
    return M;
  }
}

// 行列の和
MatrixMN operator+(const MatrixMN& M1, const MatrixMN& M2)
{
  if ((M1.dim_i == M2.dim_i) && (M1.dim_j == M2.dim_j)) {
    MatrixMN ret(M1);
  
    for (int i=0; i<ret.dim_i; ++i)
      for (int j=0; j<ret.dim_j; ++j) ret.m[i][j] += M2.m[i][j];
 
    return ret;
  }
  else {
    fprintf(stderr, "Can't addition, MatrixMN +.\n");
    return M1;
  }
}

// 行列の単項マイナス
MatrixMN operator-(const MatrixMN& M)
{
  MatrixMN ret(M.dim_i, M.dim_j);
  
  for (int i=0; i<ret.dim_i; ++i)
    for (int j=0; j<ret.dim_j; ++j) ret.m[i][j] = -M.m[i][j];
 
  return ret;
}

MatrixMN operator-(const MatrixMN& M1, const MatrixMN& M2)
{
  return M1 + (-M2);
}

// 射影行列
MatrixMN Projection(const VectorN& a)
{
  MatrixMN ret(Normalize(a), Normalize(a));

  ret = -ret;
  for (int i=0; i<ret.dim_i; ++i) ret.m[i][i] += 1.0;

  return ret;
}

// 行列の内積
double operator,(const MatrixMN& M1, const MatrixMN& M2)
{
  if ((M1.dim_i == M2.dim_i) && (M1.dim_j == M2.dim_j)) {
    double x = 0.0;

    for (int i=0; i<M1.dim_i; ++i) 
      for (int j=0; j<M1.dim_j; ++j) x += M1.m[i][j] * M2.m[i][j];

    return x;
  }
  else {
    fprintf(stderr, "Dimension error, MatrixMN Inner Product.\n");
    return 0.0;
  }
}

// 行列のノルム
double Norm(const MatrixMN& M)
{
  double x = 0.0;

  for (int i=0; i<M.dim_i; ++i) 
    for (int j=0; j<M.dim_j; ++j) x += sqr(M.m[i][j]);

  return sqrt(x);
}

// 正規化
MatrixMN Normalize(const MatrixMN& M)
{
  return M / Norm(M);
}

// ランク
int Rank(const MatrixMN& M)
{
  MatrixMN X, U;
  VectorN  L;
  int      i, ret = 0;

  if(M.dim_i >= M.dim_j) {
    X.Init(M.dim_j,M.dim_j);
    U.Init(M.dim_j,M.dim_j);
    L.Init(M.dim_j);
    X = Trans(M)*M; }
  else {
    X.Init(M.dim_i,M.dim_i);
    U.Init(M.dim_i,M.dim_i);
    L.Init(M.dim_i);
    X = M*Trans(M); }

  X.Householder(L,U);

  for (i=0; i<L.Dim(); ++i) 
    if (L[i] > eps) ++ret;

  return ret;
}

// 行列の転置
MatrixMN Trans(const MatrixMN& M)
{
  MatrixMN ret(M.dim_j, M.dim_i);
  
  for (int i=0; i<ret.dim_i; ++i) 
    for (int j=0; j<ret.dim_j; ++j) ret.m[i][j] = M.m[j][i];

  return ret;
}

// トレース
double Trace(const MatrixMN& M)
{
  double ret = 0.0;

  if (M.dim_i == M.dim_j) for (int i=0; i<M.dim_i; ++i)  ret += M.m[i][i];

  return ret;
}

// 逆行列
MatrixMN Inverse(const MatrixMN& M)
{
  if (M.dim_i == M.dim_j) {
    int    i, j;
    double det;

    double **mm, *mmp;
    mm = new double*[M.dim_i];
    mm[0] = mmp = new double[M.dim_i*M.dim_j];
    for(i=1; i<M.dim_i; i++) mm[i] = mm[i-1] + M.dim_j;

    double **inv, *invp;
    inv = new double*[M.dim_i];
    inv[0] = invp = new double[M.dim_i*M.dim_j];
    for(i=1; i<M.dim_i; i++) inv[i] = inv[i-1] + M.dim_j;

    MatrixMN  ret(M.dim_i, M.dim_j);

    for(i=0; i<M.dim_i; ++i)
      for(j=0; j<M.dim_i; ++j)
	 mm[i][j] = M.m[i][j];

    det = _alg_matinv(M.dim_i,mm,inv);

    for(i=0; i<ret.dim_i; ++i)
      for(j=0; j<ret.dim_i; ++j)
	 ret[i][j] = inv[i][j];

    delete [] mmp;  delete [] mm;
    delete [] invp; delete [] inv;
    return ret;
  }
  else {
    fprintf(stderr, "MatrixMN isn't regular.\n");
    return MatrixMN(0, 0);
  }
}

// 一般逆行列
MatrixMN GeneralInverse(const MatrixMN& M, int times) 
{
  MatrixMN  V(M.dim_i, M.dim_j), ret(M.dim_i, M.dim_j);
  VectorN   D(M.dim_i);
  int       t;

  if (times == 0) t = Rank(M);
  else            t = times;
  
//  if (fabs(Det(M)) <= eps)
    if (M.Householder(D, V)) {
      for (int i=0; i<t; ++i)
	if (fabs(D[i]) >= eps)
	  ret += 1.0 / D[i] * MatrixMN(V(i), V(i));
      return ret;
    }
    else {
      fprintf(stderr, "MatrixMN is not regular, General Inverse.\n");
      return M;
    }
//  else
//    return Inverse(M);
}

// 行列式
double Det(const MatrixMN& M)
{
  if (M.dim_i == M.dim_j) {
    int    i, j;
    double ret;

    double **mm, *mmp;
    mm = new double*[M.dim_i];
    mm[0] = mmp = new double[M.dim_i*M.dim_j];
    for(i=1; i<M.dim_i; i++) mm[i] = mm[i-1] + M.dim_j;

    for(i=0; i<M.dim_i; ++i)
      for(j=0; j<M.dim_i; ++j)
	 mm[i][j] = M.m[i][j];
    
    ret = _alg_det(M.dim_i,mm);

    delete [] mmp;  delete [] mm;

    return ret;
  }
  else {
    fprintf(stderr, "MatrixMN is not regular, Det.\n");
    return 0.0;
  }
}

// // 余因子行列
// MatrixMN Cofactor(const MatrixMN& M)
// {
//   MatrixMN ret(M.dim_i, M.dim_j), tmp(M.dim_i - 1, M.dim_j - 1);
// 
//   for (int i=0; i<M.dim_i; ++i)
//     for (int j=0; j<M.dim_j; ++j) {
//       for (int k=0, tmpi=0; k<M.dim_i; ++k)
// 	for (int l=0, tmpj=0; l<M.dim_j; ++l)
// 	  if ((k != i) || (l != j)) tmp.m[tmpi++][tmpj++] = M.m[k][l];
//       ret.m[j][i] = pow(-1.0, (double)(i + j)) * Det(tmp);
//     }
// 
//   return ret;
// }

// 固有値・固有ベクトル
int MatrixMN::Householder(VectorN &D, MatrixMN &V) const
{
  int ret = 0;

//  if (Symmetric()) {
    int    i, j;
    double **mm, *mmp, *d;

    mm = new double*[dim_i];
    mm[0] = mmp = new double[dim_i*dim_j];
    for(i=1; i<dim_i; i++) mm[i] = mm[i-1] + dim_j;
    d = new double[dim_i];

    for(i=0; i<dim_i; ++i)
      for(j=0; j<dim_i; ++j)
	 mm[i][j] = m[i][j];

    ret = _alg_eigen(dim_i,mm,d);

    for(i=0; i<dim_i-1; ++i)
      for(j=i+1; j<dim_i; ++j) {
         double wk = mm[i][j];
         mm[i][j]  = mm[j][i];
         mm[j][i]  = wk; }

    D = VectorN(d, dim_i);
    V = MatrixMN(mm, dim_i, dim_j);

    delete [] d; delete [] mmp; delete [] mm;
//  }
  
  if(ret == EXIT_SUCCESS) {
    ret = 1; }
  else {
    ret = 0; }
  return ret;
}

// // 特異値分解 Ｋ＝ＶＬＵt
// int MatrixMN::Svdecomp(MatrixMN& V, MatrixMN& L, MatrixMN& U) const
// {
//   int ret = 0;
//   
//   double *v = new double[dim_i * dim_i],
//          *l = new double[dim_i * dim_j],
//          *u = new double[dim_j * dim_j],
//          *mm  = new double[dim_i * dim_j];
// 
//   for (int i=0; i<dim_i; ++i)
//     for (int j=0; j<dim_j; ++j) mm[i * dim_j + j] = m[i][j];
//   
//   ret = MTXsvdecomp(dim_i, dim_j, mm, v, l, u);
// 
//   V = MatrixMN(v, dim_i, dim_i);
//   L = MatrixMN(l, dim_i, dim_j);
//   U = MatrixMN(u, dim_j, dim_j);
//   
//   delete [] v; delete [] l; delete [] u; delete [] mm;
//   
//   return ret;
// }



