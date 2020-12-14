#include <statcomp.h>
#include <tool.h>
#include "extra.h"

/*
 *����F�Ȥ�����������ʬ���ƥ󥽥�V_0F��Ϳ���ơ�
 *F��detF=0�Ȥʤ�褦�˺�Ŭ����������
 */
Matrix Correct(Matrix F, Tensor3333 V_0F, double eps_d = 1.0e-12)
{
  /*******************************************************
   *����
   *  Matrix      F      (����)����               input
   *  Tensor3333  V_0F   F����������ʬ���ƥ󥽥�  input
   *  double      eps_d  ȿ����λ���             input
   *                     �ǥե������: 1.0e-12
   *
   *�֤���
   *  ������������
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

      //co_F�򹹿�����
      F_trs = Trans(Cofactor(co_F));
      T_prd = v0f*F_trs;
      /******Eq (23) ******/
      co_F = Normalize(co_F - Det(co_F)*T_prd/(F_trs, T_prd));

      //co_F�μͱƥƥ󥽥�P��׻�����
      /******Eq (24) ******/
      P = ProjTensor(co_F);

//      //ư���ǧ�Ѥν���
//      printf("PF = \n");
//      (P*co_F).Print();
//      printf("\n");
      //0�ˤʤ뤳�Ȥ��ǧ�Ѥ�
      v0fp = v0f;

      //co_F����������ʬ���ƥ󥽥�v0f�򹹿�����
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

      //ư���ǧ�Ѥν���
      //printf("V[F]F = \n");
      //(v0f*co_F).Print(); //v0f*co_F��3333�ƥ󥽥��33����(=33����)
    }

  return co_F;
}
