#include <statcomp.h>
#include <tool.h>

//3333テンソルのトレースを計算する
double Trace(Tensor3333 T)
{
  double answer = 0.0;
  if ((T.DimI() == 3) && (T.DimJ() == 3) && 
      (T.DimK() == 3) && (T.DimL() == 3))
    {
      for (int m = 0; m < 3; m++)
	for (int n = 0; n < 3; n++)
	  answer += T(m, n, m ,n);
    }
  else
    {
      fprintf(stderr, "Tensor3333 dimension error in Trace.\n");
      exit(1);
    }
  
  return answer;
}


//射影テンソルを計算する
Tensor3333 ProjTensor(const Matrix& M)
{
  Tensor3333 P;

  if (M.DimI() == 3 && M.DimJ() == 3)
    {
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++)
              {
                P(i, j, k, l)
                  = KroneckerDelta(i, k)*KroneckerDelta(j, l) 
                  - M[i][j]*M[k][l];
              }
    }
  else
    {
      fprintf(stderr, "Matrix dimension error in ProjTensor.\n");
      exit(1);
    }

  return P;
}


//疑似射影テンソルを計算する
Tensor3333 Pseudo(Matrix M)
{
  Tensor3333 P;
  double Mnorm;
  
  if (M.DimI() == 3 && M.DimJ() == 3)
    {
      Mnorm = pow(Norm(M), 2);
      for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
	  for (int k = 0; k < 3; k++)
	    for (int l = 0; l < 3; l++)
	      {
		P(i, j, k, l)
		  = KroneckerDelta(i, k)*KroneckerDelta(j, l) 
		  - M[i][j]*M[k][l]/Mnorm;
	      } 
    }
  else
    {
      fprintf(stderr, "Matrix dimension error in Pseudo.\n");
      exit(1);
    }

  return P;      
}
