#include "statcomp.h"
#include "tensor.h"

//3333テンソル->99行列
MatrixMN Type99(const Tensor3333& tens)
{
  if ((tens.DimI() == 3) && (tens.DimJ() == 3) &&
      (tens.DimK() == 3) && (tens.DimL() == 3))
    {
      MatrixMN ret(9, 9);
      
      for (int i = 0; i <= 2; i++)
	for (int j = 0; j <= 2; j++)
	  for (int k = 0; k <= 2; k++)
	    for (int l = 0; l <= 2; l++)
	      ret[3*i+j][3*k+l] = tens(i, j, k, l);
      return ret;
    }
    else
      {
	printf("Function Type99 error!\n");
	exit(1);
      }
}


//(33)(33)テンソル->66行列
MatrixMN Type66(const Tensor3333& tens)
{
  if ((tens.DimI() == 3) && (tens.DimJ() == 3) &&
      (tens.DimK() == 3) && (tens.DimL() == 3))
    {
      MatrixMN ret(6, 6);
      double root2 = sqrt(2.0);
      
      ret[0][0] = tens(0, 0, 0, 0);
      ret[1][0] = tens(1, 1, 0, 0);
      ret[2][0] = tens(2, 2, 0, 0);
      ret[3][0] = root2*tens(1, 2, 0, 0);
      ret[4][0] = root2*tens(2, 0, 0, 0);
      ret[5][0] = root2*tens(0, 1, 0, 0);
	
      ret[0][1] = tens(0, 0, 1, 1);
      ret[1][1] = tens(1, 1, 1, 1);
      ret[2][1] = tens(2, 2, 1, 1);
      ret[3][1] = root2*tens(1, 2, 1, 1);
      ret[4][1] = root2*tens(2, 0, 1, 1);
      ret[5][1] = root2*tens(0, 1, 1, 1);
	
      ret[0][2] = tens(0, 0, 2, 2);
      ret[1][2] = tens(1, 1, 2, 2);
      ret[2][2] = tens(2, 2, 2, 2);
      ret[3][2] = root2*tens(1, 2, 2, 2);
      ret[4][2] = root2*tens(2, 0, 2, 2);
      ret[5][2] = root2*tens(0, 1, 2, 2);
	
      ret[0][3] = root2*tens(0, 0, 1, 2);
      ret[1][3] = root2*tens(1, 1, 1, 2);
      ret[2][3] = root2*tens(2, 2, 1, 2);
      ret[3][3] = 2.0*tens(1, 2, 1, 2);
      ret[4][3] = 2.0*tens(2, 0, 1, 2);
      ret[5][3] = 2.0*tens(0, 1, 1, 2);
	
      ret[0][4] = root2*tens(0, 0, 2, 0);
      ret[1][4] = root2*tens(1, 1, 2, 0);
      ret[2][4] = root2*tens(2, 2, 2, 0);
      ret[3][4] = 2.0*tens(1, 2, 2, 0);
      ret[4][4] = 2.0*tens(2, 0, 2, 0);
      ret[5][4] = 2.0*tens(0, 1, 2, 0);
	
      ret[0][5] = root2*tens(0, 0, 0, 1);
      ret[1][5] = root2*tens(1, 1, 0, 1);
      ret[2][5] = root2*tens(2, 2, 0, 1);
      ret[3][5] = 2.0*tens(1, 2, 0, 1);
      ret[4][5] = 2.0*tens(2, 0, 0, 1);
      ret[5][5] = 2.0*tens(0, 1, 0, 1);
	
      return ret;
    }
  else
    {
      printf("Function Type66 error!\n");
      exit(1);
    }
}


//6ベクトル->(33)行列
Matrix Type33sym(const VectorN& vect)
{
  if (vect.Dim() == 6)
    {
      double root2 = sqrt(2.0);
      
      return Matrix(root2*vect[0], vect[5], vect[4],
		    vect[5], root2*vect[1], vect[3],
		    vect[4], vect[3], root2*vect[2])/root2;
    }
  else
      {
	printf("Function Type33sym error!\n");
	exit(1);
      }
}


//66行列->(33)(33)テンソル
Tensor3333 Type3333sym(const MatrixMN& M)
{
  if ((M.DimI() == 6) && (M.DimJ() == 6))
    {
      Tensor3333 T;
      double root2 = sqrt(2.0);
      
      T(0, 0, 0, 0) = M[0][0];
      T(1, 1, 0, 0) = M[1][0];
      T(2, 2, 0, 0) = M[2][0];
      T(1, 2, 0, 0) = M[3][0]/root2;
      T(2, 0, 0, 0) = M[4][0]/root2;
      T(0, 1, 0, 0) = M[5][0]/root2;

      T(0, 0, 1, 1) = M[0][1];
      T(1, 1, 1, 1) = M[1][1];
      T(2, 2, 1, 1) = M[2][1];
      T(1, 2, 1, 1) = M[3][1]/root2;
      T(2, 0, 1, 1) = M[4][1]/root2;
      T(0, 1, 1, 1) = M[5][1]/root2;

      T(0, 0, 2, 2) = M[0][2];
      T(1, 1, 2, 2) = M[1][2];
      T(2, 2, 2, 2) = M[2][2];
      T(1, 2, 2, 2) = M[3][2]/root2;
      T(2, 0, 2, 2) = M[4][2]/root2;
      T(0, 1, 2, 2) = M[5][2]/root2;

      T(0, 0, 1, 2) = M[0][3]/root2;
      T(1, 1, 1, 2) = M[1][3]/root2;
      T(2, 2, 1, 2) = M[2][3]/root2;
      T(1, 2, 1, 2) = M[3][3]/2.0;
      T(2, 0, 1, 2) = M[4][3]/2.0;
      T(0, 1, 1, 2) = M[5][3]/2.0;

      T(0, 0, 2, 0) = M[0][4]/root2;
      T(1, 1, 2, 0) = M[1][4]/root2;
      T(2, 2, 2, 0) = M[2][4]/root2;
      T(1, 2, 2, 0) = M[3][4]/2.0;
      T(2, 0, 2, 0) = M[4][4]/2.0;
      T(0, 1, 2, 0) = M[5][4]/2.0;

      T(0, 0, 0, 1) = M[0][5]/root2;
      T(1, 1, 0, 1) = M[1][5]/root2;
      T(2, 2, 0, 1) = M[2][5]/root2;
      T(1, 2, 0, 1) = M[3][5]/2.0;
      T(2, 0, 0, 1) = M[4][5]/2.0;
      T(0, 1, 0, 1) = M[5][5]/2.0;

      for (int i = 0; i <= 2; i++)
	for (int j = 0; j <= 2; j++)
	  for (int k = 0; k <= 2; k++)
	    for (int l = 0; l <= 2; l++)
	      if (T(i, j, k, l) != 0.0)
		{
		  T(i, j, l, k) = T(i, j, k, l);
		  T(j, i, k, l) = T(i, j, k, l);
		  T(j, i, l, k) = T(i, j, k, l);
		}
      
      return T;
    }
  else
    {
      printf("Function Type3333sym error!\n");
      exit(1);
    }
}


//9ベクトル->33行列
Matrix Type33(const VectorN& vec)
{
  if (vec.Dim() == 9)
    {
      Matrix ret;
      for (int i = 0; i <= 2; i++)
	for (int j = 0; j <= 2; j++)
	  ret[i][j] = vec[3*i + j];

	return ret;
    }
  else
    {
      printf("Function Type33 error!\n");
      exit(1);
    }
}


//33行列->9ベクトル
VectorN Type9(const Matrix& mat)
{
  if (mat.DimI() == 3 && mat.DimJ() == 3)
    {
      VectorN ret(9);
      for (int i = 0; i <= 2; i++)
	for (int j = 0; j <= 2; j++)
	  ret[3*i+j] = mat[i][j];

      return ret;
    }
  else
    {
      printf("Function Type9 error!\n");
      exit(1);
    }
}


//99行列->3333テンソル
Tensor3333 Type3333(const MatrixMN& M)
{
  Tensor3333 T;
  int kap, lam;

  for (int i = 1; i <= 3; i++)
    for (int j = 1; j <= 3; j++)
      for (int k = 1; k <= 3; k++)
	for (int l = 1; l <= 3; l++) 
	  {
	    kap = 3*(i-1)+j;
	    lam = 3*(k-1)+l;
	    T(i-1, j-1, k-1, l-1) = M[kap-1][lam-1];
	  }
  
  return T;
}
