#include <statcomp.h>
#include <tool.h>
#include "extra.h"

/*
 *行列Fとその正規化共分散テンソルV_0Fを与えて、
 *FがdetF=0となるように最適に補正する
 */
Matrix Correct(Matrix F, Tensor3333 V_0F, double eps_d = 1.0e-12)
{
  /*******************************************************
   *引数
   *  Matrix      F      (基礎)行列               input
   *  Tensor3333  V_0F   Fの正規化共分散テンソル  input
   *  double      eps_d  反復終了条件             input
   *                     デフォルト値: 1.0e-12
   *
   *返り値
   *  補正した行列
   *
   ******************************************************/

  int i, j, k, l, m, n, p, q, iter = 0;
  Matrix F_trs, T_prd, co_F = F;
  Tensor3333 P, v0fp, v0f = V_0F;

  while (fabs(Det(co_F)) > eps_d)
    {
      iter++;
      if (iter > 20)
	{
	  fprintf(stderr, "Iteration exceeded 20 times.\n");
	  exit(1);
	}

      //co_Fを更新する
      F_trs = Trans(Cofactor(co_F));
      T_prd = v0f*F_trs;
      /******Eq (23) ******/
      co_F = Normalize(co_F - Det(co_F)*T_prd/(F_trs, T_prd));

      //co_Fの射影テンソルPを計算する
      /******Eq (24) ******/
      P = ProjTensor(co_F);

//      //動作確認用の出力
//      printf("PF = \n");
//      (P*co_F).Print();
//      printf("\n");
      //0になることを確認済み
      v0fp = v0f;

      //co_Fの正規化共分散テンソルv0fを更新する
      for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	  for (k = 0; k < 3; k++)
	    for (l = 0; l < 3; l++)
	      for (m = 0; m < 3; m++)
		for (n = 0; n < 3; n++)
		  for (p = 0; p < 3; p++)
		    for (q = 0; q < 3; q++)
		      {
		          /******Eq (25) ******/
			v0f(i, j, k, l)
			  += P(i, j, m, n)*P(k, l, p, q)*v0fp(m, n, p, q);
		      }

      //動作確認用の出力
      //printf("V[F]F = \n");
      //(v0f*co_F).Print(); //v0f*co_Fは3333テンソル×33行列(=33行列)
    }

  return co_F;
}
