#include <statcomp.h>
#include <tool.h>
#include "correct.h"
#include "extra.h"

/*
 *3333テンソルM_hatの9個の固有値VectorN(9)と対応する固有行列の
 *正規直交系Matrix*, 及び二乗ノイズレベルの推定値e_hat2を与えて、
 *補正した基礎行列F, その標準偏位F_plus, F_minusを返す。
 *  ※(PDP = primary deviation pair)
 */

void PDP(int point_num, VectorN eigen_vec, Matrix* eigen_mat,
	 double e_hat2, Matrix& F, Matrix& F_plus, Matrix& F_minus,
	 int crr_flag = 1)
{
  /****************************************************************
   *引数
   *  int       point_num   データ点数                input
   *  VectorN   eigen_vec   M_hatの9個の固有値        input
   *                        9ベクトルの各要素
   *  Matrix    *eigen_mat  M_hatの9個の固有行列      input
   *                        eigen_vecの各要素に対応
   *  double    e_hat2      二乗ノイズレベルの推定値  input
   *  Matrix&   F           補正した基礎行列          output
   *  Matrix&   F_plus      Fの標準偏位               output
   *  Matrix&   F_minus     Fの   〃                  output
   *  int       crr_flag    最適補正実行フラグ        input
   *                        1のときFの最適補正実行
   *                        デフォルト値: 1
   *
   ****************************************************************/

  int i, j, k, l, m, n, p, q;
  Tensor3333 M_hat, V_0F, W, Pi, P;


  /******M_hatの最小固有値に対応する固有行列をFとする******/
  F = eigen_mat[8];


  /******Fの正規化共分散テンソルV_0Fを計算する******/
  for (int a = 0; a < 8; a++)
  /******Eq (22) ******/
    V_0F = V_0F + TensorProduct(eigen_mat[a], eigen_mat[a])/eigen_vec[a];
  V_0F = V_0F/(double)point_num;


  /******FをdetF=0となるように補正する******/
  if (crr_flag == 1)
    F = Correct(F, V_0F, 1.0e-12);


  /******モーメントテンソルM_hatを計算し直す******/
  for (int a = 0; a < 8; a++)
    M_hat = M_hat + eigen_vec[a]*TensorProduct(eigen_mat[a], eigen_mat[a]);

  /******テンソルWを計算する******/
  /******Eq. 10 ******/
  Pi = Pseudo(Cofactor(F));
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
        for (l = 0; l < 3; l++)
          for (m = 0; m < 3; m++)
            for (n = 0; n < 3; n++)
              for (p = 0; p < 3; p++)
                for (q = 0; q < 3; q++)
                  {
                    W(i, j, k, l) = W(i, j, k, l) +
                      Pi(i, j, m, n)*Pi(k, l, p, q)*M_hat(m, n, p, q);
                  }


  /******テンソルWの9個の固有値・固有行列を求める******/
  VectorN eigen_mu(9);
  Matrix eigen_G[9];
  if (W.Eigen(eigen_mu, eigen_G))
    {
      //標準偏位F_plus, F_minusを計算する
      double mu = sqrt(e_hat2/(eigen_mu[6]*point_num));
      F_plus  = Normalize(F + mu*eigen_G[6]);
      F_minus = Normalize(F - mu*eigen_G[6]);
    }
  else
    {
      printf("Cannot calculate eigenvalues of Tensor 'W'\n");
      exit(1);
    }
}
