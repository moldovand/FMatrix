#ifndef _Correct_h_
#define _Correct_h_

/*
 *����F�Ȥ�����������ʬ���ƥ󥽥�V_0F��Ϳ���ơ�
 *F��detF=0�Ȥʤ�褦�˺�Ŭ����������
 */
Matrix Correct(Matrix F, Tensor3333 V_0F, double eps_d = 1.0e-12);

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


#endif //_Correct_h_
