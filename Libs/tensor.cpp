#include "statcomp.h"

#include "tensor.h"
#include "transform.h"

//ijklテンソルクラス

//メモリ確保関数(4次元配列用)
static void memkeep(double****& t, int in, int jn, int kn, int ln)
{
  //メモリ確保
  t = new double***[in];
  for (int i = 0; i < in; i++)
    t[i] = new double**[jn];
  for (int i = 0; i < in; i++)
    for (int j = 0; j < jn; j++)
      t[i][j] = new double*[kn];
  for (int i = 0; i < in; i++)
    for (int j = 0; j < jn; j++)
      for (int k = 0; k < kn; k++)
	t[i][j][k] = new double[ln];

  //全ての要素を0に初期化
  for (int i = 0; i < in; i++)
    for (int j = 0; j < jn; j++)
      for (int k = 0; k < kn; k++)
	for (int l = 0; l < ln; l++)
	  t[i][j][k][l] = 0.0;
}

//1st constructor
Tensor3333::Tensor3333(int i, int j, int k, int l)
{
  in = i; jn = j; kn = k; ln = l;
  memkeep(t, in, jn, kn, ln);
}

//copy constructor
Tensor3333::Tensor3333(const Tensor3333& arg)
{
  in = arg.in; jn = arg.jn; kn = arg.kn; ln = arg.ln;
  memkeep(t, in, jn, kn, ln);
  for (int i = 0; i < in; i++)
    for (int j = 0; j < jn; j++)
      for (int k = 0; k < kn; k++)
	for (int l = 0; l < ln; l++)
	  t[i][j][k][l] = arg.t[i][j][k][l];
}

//destructor
Tensor3333::~Tensor3333()
{
  for (int i = 0; i < in; i++)
    for (int j = 0; j < jn; j++)
      for (int k = 0; k < kn; k++)
	delete [] t[i][j][k];
  for (int i = 0; i < in; i++)
    for (int j = 0; j < jn; j++)
      delete [] t[i][j];
  for (int i = 0; i < in; i++)
    delete [] t[i];
  delete [] t;
}

//要素取り出し
double& Tensor3333::operator()(int i, int j, int k, int l)
{
  if (i >= 0 && i < in && j >= 0 && j < jn &&
      k >= 0 && k < kn && l >= 0 && l < ln)
    return t[i][j][k][l];
  else
    {
      printf("Function operator[] error!\n");
      exit(1);
    }
}

const double& Tensor3333::operator()(int i, int j, int k, int l) const
{
  if (i >= 0 && i < in && j >= 0 && j < jn &&
      k >= 0 && k < kn && l >= 0 && l < ln)
    return t[i][j][k][l];
  else
    {
      printf("Function operator[] error!\n");
      exit(1);
    }
}

//初期化
int Tensor3333::Init()
{
  for (int i = 0; i < in; i++)
    for (int j = 0; j < jn; j++)
      for (int k = 0; k < kn; k++)
        for (int l = 0; l < ln; l++)
          t[i][j][k][l] = 0.0;

  return true;
}

//代入演算
const Tensor3333& Tensor3333::operator=(const Tensor3333& rhs)
{
  if (in == rhs.in && jn == rhs.jn &&
      kn == rhs.kn && ln == rhs.ln)
    {
      for (int i = 0; i < in; i++)
	for (int j = 0; j < jn; j++)
	  for (int k = 0; k < kn; k++)
	    for (int l = 0; l <= 2; l++)
	      t[i][j][k][l] = rhs.t[i][j][k][l];

      return *this;
    }
  else
    {
      printf("Function Tensor3333::operator= error!\n");
      exit(1);
    }
}

//足し算
Tensor3333 Tensor3333::operator+(const Tensor3333& arg) const
{
  if (in == arg.in && jn == arg.jn &&
      kn == arg.kn && ln == arg.ln)
    {
      Tensor3333 answer(in, jn, kn, ln);
      for (int i = 0; i < in; i++)
	for (int j = 0; j < jn; j++)
	  for (int k = 0; k < kn; k++)
	    for (int l = 0; l < ln; l++)
	      answer.t[i][j][k][l]
		= t[i][j][k][l] + arg.t[i][j][k][l];
      return answer;
    }
  else
    {
      printf("Function Tensor3333::operator+or- error!\n");
      exit(1);
    }
}

//引き算
Tensor3333 Tensor3333::operator-(const Tensor3333& arg) const
{
  return *this + (-arg);
}

//かけ算
Tensor3333 Tensor3333::operator*(double r) const
{
  Tensor3333 answer(in, jn, kn, ln);

  for (int i = 0; i < in; i++)
    for (int j = 0; j < jn; j++)
      for (int k = 0; k < kn; k++)
	for (int l = 0; l < ln; l++)
	  answer.t[i][j][k][l] = t[i][j][k][l]*r;

  return answer;
}

Tensor3333 operator*(double r, const Tensor3333& tens)
{
  return tens * r;
}

//かけ算(3333テンソル*33行列用)
Matrix Tensor3333::operator*(const Matrix& arg)
{
  if (in == 3 && jn == 3 && kn == 3 && ln == 3)
    return Type33(Type99(*this)*Type9(arg));
  else
    {
      printf("Function Tensor3333::operator* error!\n");
      exit(1);
    }
}

//割算
Tensor3333 Tensor3333::operator/(double r)
{
  return (*this)*(1.0/r);
}

//符号換え
Tensor3333 operator-(const Tensor3333& arg)
{
  Tensor3333 answer(arg.in, arg.jn, arg.kn, arg.ln);

  for (int i = 0; i < arg.in; i++)
    for (int j = 0; j < arg.jn; j++)
      for (int k = 0; k < arg.kn; k++)
	for (int l = 0; l < arg.ln; l++)
	  answer.t[i][j][k][l] = -arg.t[i][j][k][l];

  return answer;
}

//複製による3333テンソルの要素生成
void Tensor3333::Fill()
{
  for (int i = 0; i <= 2; i++)
    for (int j = 0; j <= 2; j++)
      for (int k = 0; k <= 2; k++)
        for (int l = 0; l <= 2; l++)
          if (3*(i-1)+j > 3*(k-1)+l)
            t[i][j][k][l] = t[k][l][i][j];
}

//複製による(33)(33)テンソルの要素生成
void Tensor3333::FillSym()
{
  int i, j, k, l;

  for (i = 0; i <= 2; i++)
    for (j = 0; j <= 2; j++)
      for (k = 0; k <= 2; k++)
        for (l = 0; l <= 2; l++)
          if (i <= j)
            if (k <= l)
              if (3*(i-1)+j > 3*(k-1)+l)
                t[i][j][k][l] = t[k][l][i][j];

  for (i = 0; i <= 2; i++)
    for (j = 0; j <= 2; j++)
      for (k = 0; k <= 2; k++)
        for (l = 0; l <= 2; l++)
          if (i > j)
            if (k <= l)
              t[i][j][k][l] = t[j][i][k][l];

  for (i = 0; i <= 2; i++)
    for (j = 0; j <= 2; j++)
      for (k = 0; k <= 2; k++)
        for (l = 0; l <= 2; l++)
          if (k > l)
            t[i][j][k][l] = t[i][j][l][k];
}

//(33)(33)テンソルの固有値と固有行列を求める
int Tensor3333::EigenSym(VectorN& val, Matrix* mat)
{
  MatrixMN M = Type66(*this);
  MatrixMN eigenvectors(6, 6);  //固有ベクトル

  if (M.Householder(val, eigenvectors) == 0)
    {
      printf("Function Tensor3333::EigenSym error!\n");
      return 0;
    }
  else
    {
      for (int i = 0; i <= 5; i++)
	mat[i] = Type33sym(eigenvectors(i));

      return 1;
    }
}

//3333テンソルの固有値と固有行列を求める
int Tensor3333::Eigen(VectorN& val, Matrix* mat)
{
  MatrixMN M = Type99(*this);
  MatrixMN eigenvectors(9, 9);  //固有ベクトル

  if (M.Householder(val, eigenvectors) == 0)
    {
      printf("Function Tensor3333::Eigen error!\n");
      return 0;
    }
  else
    {
      for (int i = 0; i <= 8; i++)
	mat[i] = Type33(eigenvectors(i));

      return 1;
    }
}


//要素を出力する
void Tensor3333::Print()
{
  for (int i = 0; i < in; i++)
    for (int j = 0; j < jn; j++)
      for (int k = 0; k < kn; k++)
	{
	  for (int l = 0; l < ln; l++)
	    printf("t[%d][%d][%d][%d] = %g  ",
		   i, j, k, l, t[i][j][k][l]);
	  printf("\n");
	}
}

void Tensor3333::Fprint(FILE* fp)
{
  for (int i = 0; i < in; i++)
    for (int j = 0; j < jn; j++)
      for (int k = 0; k < kn; k++)
	{
	  for (int l = 0; l < ln; l++)
	    fprintf(fp, "t[%d][%d][%d][%d] = %g  ",
		    i, j, k, l, t[i][j][k][l]);

	  fprintf(fp, "\n");
	}
}

//(33)(33)テンソルの一般逆を求める
Tensor3333 GeneralInverseSym(const Tensor3333& T, int times)
{
  if (T.in == 3 && T.jn == 3 && T.kn == 3 && T.ln == 3)
    return Type3333sym(GeneralInverse(Type66(T), times));
  else
    {
      printf("Function Tensor3333::GeneralInverseSym error!\n");
      exit(1);
    }
}

//3333テンソルの一般逆を求める
Tensor3333 GeneralInverse(const Tensor3333& T, int times)
{
  if (T.in == 3 && T.jn == 3 && T.kn == 3 && T.ln == 3)
    return Type3333(GeneralInverse(Type99(T), times));
  else
    {
      printf("Function Tensor3333::GeneralInverse error!\n");
      exit(1);
    }
}

//テンソル積(引数の２つの行列は３×３に限る)
Tensor3333 TensorProduct(const Matrix& m1, const Matrix& m2)
{
  Tensor3333 ret;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
	for (int l = 0; l < 3; l++)
	  ret(i, j, k, l) = m1[i][j]*m2[k][l];
  return ret;
}
