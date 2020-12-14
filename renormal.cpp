//----------------------------------------------------------------------//
//  fundamental: calculate fundamental matrix by renormalization        //
//                                                                      //
//  1999/2/12                                                           //
//  Programmed by Hitoshi Mishima                                       //
//  Department of Computer Science,                                     //
//  Gunma University Faculty of Engineering                             //
//----------------------------------------------------------------------//

#include <statcomp.h>
#include <tool.h>

int Renormal(ImagePosition* P1, ImagePosition* P2, int point_num,
	     VectorN& eigen_val, Matrix* eigen_mat, double& J,
	     int dim = 1, double min = 1.0e-8, int itmax = 50)
//  P1        : array of image points(camera1).         input
//  P2        : array of image points(camera2).         input
//  point_num : number of image points.                 input
//  eigen_val : eigenvalues of moment tensor            output
//  eigen_mat : eigenmatrices of moment tensor          output
//  J         : residual of renormalization             output
//  dim       : dimension of renormalization            input
//               (1 or 2, default: 1)
//  min       : constant for convergence test           input
//               (default: 1.0e-8)
//  itmax     : maximum iteration number                input
//               (>=1, default: 50)
//
//Return value
//  iteration number
{
  int i, j, k, l, a, iter = 0;
  double lambda9, W = 1.0, c = 0.0, J2 = 0.0;
  Vector x1, x2;
  Matrix V1, V2, F9;


  /******テンソルと残差Jを初期化する******/
  Tensor3333 M_hat, M, N1, N2;
  J = 1.0e+10;

  /******各テンソルの特定の要素を計算する******/
  for (a = 0; a < point_num; a++)
    {
      //画像位置からベクトルと共分散行列を取得
      x1 = P1[a].Pos();      V1 = P1[a].Cov();
      x2 = P2[a].Pos();      V2 = P2[a].Cov();

      //各テンソルの要素計算
      for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	  for (k = 0; k < 3; k++)
	    for (l = 0; l < 3; l++)
	      {
		if (3*(i-1)+j <= 3*(k-1)+l)
		  {
		    M(i, j, k, l) = M(i, j, k, l) + x1[i]*x2[j]*x1[k]*x2[l];

		    N1(i, j, k, l) = N1(i, j, k, l) +
		      V1[i][k]*x2[j]*x2[l] + V2[j][l]*x1[i]*x1[k];

		    if (dim == 2)
		      N2(i, j, k, l) = N2(i, j, k, l) + V1[i][k]*V2[j][l];
		  }
	      }
    }


  /******各テンソルの残りの要素をすべて埋める******/
  M.Fill();
  N1.Fill();
  if (dim == 2) N2.Fill();


  /******各テンソルを更新する******/
  M = M/(double)point_num;
  N1 = N1/(double)point_num;
  if (dim == 2)
    N2 = N2/(double)point_num;
  M_hat = M;


  /******M_hatの9個の固有値・固有行列を求める******/
  if (M_hat.Eigen(eigen_val, eigen_mat) == 0)
    return 0;
  else
    {
      lambda9 = eigen_val[8]; //最小固有値
      if (lambda9 < 0.0) lambda9 *= -1.0;
      F9 = eigen_mat[8];      //対応するノルム1の固有行列
//       for (int num1 = 0; num1 < 9; num1++)
// 	printf("\t固有値 %d : %22.16f\n", num1+1, eigen_val[num1]);
//       printf("Eigenmatrix_9 = \n");
//       eigen_mat[8].LongPrint();
    }


  /******くりこみ法の反復******/
  do
    {
      iter++; //反復回数+1

      //反復回数がitmaxを越えたら終了(くりこみ法失敗)
      if (iter > itmax)
        {
//           printf("Renormalization exceed %d iteration.\n", itmax);
          return 0;
        }

      //定数cの更新
      double inp1 = (F9, N1*F9); //F9とN1*F9の内積
      double inp2 = (F9, N2*F9); //F9とN2*F9の内積

      if (dim == 2)
	{
	  double D = pow(inp1 - 2.0*c*inp2, 2) - 4.0*lambda9*inp2;
	  if (D >= 0.0)
	    c += (inp1 - 2.0*c*inp2 - sqrt(D))/(2.0*inp2);
	  else
	    c += lambda9/inp1;
	}
      else
	c = c + lambda9/inp1;

      //各テンソルの初期化
      M.Init();
      N1.Init();
      if (dim == 2)
	N2.Init();

      //各テンソルの特定の要素を計算
      for (a = 0; a < point_num; a++)
	{
	  //画像位置からベクトルと共分散行列を取得
	  x1 = P1[a].Pos();      V1 = P1[a].Cov();
	  x2 = P2[a].Pos();      V2 = P2[a].Cov();

	  //重みWの更新
	  Matrix TF9 = Trans(F9);
	  if (dim == 1)
	    W = 1.0/((x2, TF9*V1*F9*x2) + (x1, F9*V2*TF9*x1));
	  else
	    W = 1.0/((x2, TF9*V1*F9*x2) + (x1, F9*V2*TF9*x1)
		     + c*(V1*F9, F9*V2));

	  //各テンソルの要素計算
	  for (i = 0; i < 3; i++)
	    for (j = 0; j < 3; j++)
	      for (k = 0; k < 3; k++)
		for (l = 0; l < 3; l++)
		  {
		    if (3*(i-1)+j <= 3*(k-1)+l)
		      {
			M(i, j, k, l) += W*x1[i]*x2[j]*x1[k]*x2[l];

			N1(i, j, k, l) = N1(i, j, k, l) +
			  W*(V1[i][k]*x2[j]*x2[l]+ V2[j][l]*x1[i]*x1[k]);

			if (dim == 2)
			  {
			    N2(i, j, k, l)
			      += W*V1[i][k]*V2[j][l];
			  }
		      }
		  }
	}

      //各テンソルの残りの要素をすべて埋める
      M.Fill();
      N1.Fill();
      if (dim == 2) N2.Fill();

      //各テンソルの更新
      M = M/(double)(point_num);
      N1 = N1/(double)(point_num);
      M_hat = M - c*N1;
      if (dim == 2)
	{
	  N2 = N2/(double)(point_num);
	  M_hat = M_hat + c*c*N2;
	}

      //残差Jを更新
      J2 = J;
      J = (F9, M*F9);

      if (J2 < J) J = J2;
      else
	{
	  //M_hatの9個の固有値・固有行列を求める
	  if (M_hat.Eigen(eigen_val, eigen_mat) == 0)
	    return 0;
	  else
	    {
	      lambda9 = eigen_val[8]; //最小固有値
	      F9 = eigen_mat[8];      //対応するノルム1の固有行列
// 	      printf("\t反復%i回目 : 最小固有値 = %22.16f\n", iter, lambda9);
	    }
	}
    } while ((J2 > J) && ((fabs(lambda9)) > min));

  return iter;
}
