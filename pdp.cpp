#include <statcomp.h>
#include <tool.h>
#include "correct.h"
#include "extra.h"

/*
 *3333�ƥ󥽥�M_hat��9�Ĥθ�ͭ��VectorN(9)���б������ͭ�����
 *����ľ���Matrix*, �ڤ����Υ�����٥�ο�����e_hat2��Ϳ���ơ�
 *�����������ù���F, ����ɸ���а�F_plus, F_minus���֤���
 *  ��(PDP = primary deviation pair)
 */

void PDP(int point_num, VectorN eigen_vec, Matrix* eigen_mat,
	 double e_hat2, Matrix& F, Matrix& F_plus, Matrix& F_minus,
	 int crr_flag = 1)
{
  /****************************************************************
   *����
   *  int       point_num   �ǡ�������                input
   *  VectorN   eigen_vec   M_hat��9�Ĥθ�ͭ��        input
   *                        9�٥��ȥ�γ�����
   *  Matrix    *eigen_mat  M_hat��9�Ĥθ�ͭ����      input
   *                        eigen_vec�γ����Ǥ��б�
   *  double    e_hat2      ���Υ�����٥�ο�����  input
   *  Matrix&   F           �����������ù���          output
   *  Matrix&   F_plus      F��ɸ���а�               output
   *  Matrix&   F_minus     F��   ��                  output
   *  int       crr_flag    ��Ŭ�����¹ԥե饰        input
   *                        1�ΤȤ�F�κ�Ŭ�����¹�
   *                        �ǥե������: 1
   *
   ****************************************************************/

  int i, j, k, l, m, n, p, q;
  Tensor3333 M_hat, V_0F, W, Pi, P;


  /******M_hat�κǾ���ͭ�ͤ��б������ͭ�����F�Ȥ���******/
  F = eigen_mat[8];


  /******F����������ʬ���ƥ󥽥�V_0F��׻�����******/
  for (int a = 0; a < 8; a++)
  /******Eq (22) ******/
    V_0F = V_0F + TensorProduct(eigen_mat[a], eigen_mat[a])/eigen_vec[a];
  V_0F = V_0F/(double)point_num;


  /******F��detF=0�Ȥʤ�褦����������******/
  if (crr_flag == 1)
    F = Correct(F, V_0F, 1.0e-12);


  /******�⡼���ȥƥ󥽥�M_hat��׻���ľ��******/
  for (int a = 0; a < 8; a++)
    M_hat = M_hat + eigen_vec[a]*TensorProduct(eigen_mat[a], eigen_mat[a]);

  /******�ƥ󥽥�W��׻�����******/
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


  /******�ƥ󥽥�W��9�Ĥθ�ͭ�͡���ͭ��������******/
  VectorN eigen_mu(9);
  Matrix eigen_G[9];
  if (W.Eigen(eigen_mu, eigen_G))
    {
      //ɸ���а�F_plus, F_minus��׻�����
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
