//
// ���ꤳ��ˡ�ˤ����ù���Ȥ���ɸ���а̤κ�Ŭ�׻�
//
// Usage: fundamental file1 file2
// 
//      flie1  : ��ư����pixel��ɸ        ascii-file
//      file2  : ��ư���pixel��ɸ        ascii-file

#define MAIN
#include <statcomp.h>
#include <tool.h>
#include "renormal.h"
#include "pdp.h"
#include "correct.h"
#include "extra.h"
#include "filedeal.h"

int main(int argc, char *argv[])
{
  char *usg = "Usage: %s file1 file2\n";
  char *msg = "%s: Data number disagreement!\n";

  if (argc != 3)
    {
      fprintf(stderr, usg, argv[0]); 
      exit(1);
    }


  int a, data_num, data_num2, dim, itmax;
  double **d1, **d2;
  VectorN eigen_val(9);
  Matrix eigen_mat[9], F_hat, F_p, F_m, c_F;
  double J, e_hat2, min, f0;
  double xm1, ym1, Lx1, Ly1, xm2, ym2, Lx2, Ly2;

  data_num = data_num2 = 0;
  xm1 = ym1 = Lx1 = Ly1 = xm2 = ym2 = Lx2 = Ly2 = 0.0;
  f0 = 600.0; dim = 1; itmax = 50; min = 1.0e-10;


  //������ɸ�ե������������ɤ߹��ߤʤ�����������Ѵ��ν���
  Text2(argv[1], d1, data_num,  xm1, ym1, Lx1, Ly1); //��ư��
  Text2(argv[2], d2, data_num2, xm2, ym2, Lx2, Ly2); //��ư��

  if (data_num != data_num2)
    {
      fprintf(stderr, msg, argv[0]);
      exit(1);
    }

  
  //��������󥰰���(�ǥե���Ⱦ�����Υ)�γ���
  f0 = Lx1;
  if (f0 < Ly1)
    f0 = Ly1;
  if (f0 < Lx2)
    f0 = Lx2;
  if (f0 < Ly2)
    f0 = Ly2;
  f0 *= 2.0;


  //�������֤�����
  ImagePosition P1[data_num], P2[data_num];
  for (a = 0; a < data_num; a++)
    {
       P1[a] = ImagePosition(ym1-d1[a][1], d1[a][0]-xm1, f0);
       P2[a] = ImagePosition(ym2-d2[a][1], d2[a][0]-xm2, f0);
    }

  for (a = 0; a < data_num; ++a) delete [] d1[a];
  delete [] d1;
  for (a = 0; a < data_num; ++a) delete [] d2[a];
  delete [] d2;


  //���ϥե����ޥå��Ѵ��ѹ���A1, A2�η׻�
  Matrix A1, A2;
  A1 = Matrix(0,    -1/f0,  ym1/f0,
  	      1/f0,  0,    -xm1/f0,
  	      0,     0,       1);
  A2 = Matrix(0,    -1/f0,  ym2/f0,
 	      1/f0,  0,    -xm2/f0,
 	      0,     0,       1);


  //���ꤳ��ˡ�ˤ����ù�������
  if (Renormal(P1, P2, data_num, eigen_val, eigen_mat, J, dim, min, itmax))
    {
      //���ù�����ಽȽ��
      J = fabs(J);
      e_hat2 = J/(1.0 - 8/(double)data_num);
      if (e_hat2 > eigen_val[7]*data_num)
        {
          printf("Warning: Degeneracy!\n");
          exit(1);
        }
      else
	{
	  //��Ŭ�����ڤ�ɸ���а̷׻���Ԥ�
	  PDP(data_num, eigen_val, eigen_mat, e_hat2, c_F, F_p, F_m);
	  c_F = Normalize(Trans(A1)*c_F*A2);
	  c_F.Print(stdout, "F  = ");  putchar('\n');
	  F_p = Normalize(Trans(A1)*F_p*A2);
	  F_p.Print(stdout, "F+ = ");  putchar('\n');
	  F_m = Normalize(Trans(A1)*F_m*A2);
	  F_m.Print(stdout, "F- = ");  putchar('\n');
	}
      //printf("Calculation succeeded.\n");
    }
  else
    printf("Calculation failure.\n");
  

  exit(0);
}

