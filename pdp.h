#ifndef _pdp_h_
#define _pdp_h_

/*
 *3333�ƥ󥽥�M_hat��9�Ĥθ�ͭ��VectorN(9)���б������ͭ�����
 *����ľ���Matrix*, �ڤ����Υ�����٥�ο�����e_hat2��Ϳ���ơ�
 *�����������ù���F, ����ɸ���а�F_plus, F_minus���֤���
 *  ��(PDP = primary deviation pair)
 */

void PDP(int point_num, VectorN eigen_vec, Matrix* eigen_mat, 
	 double e_hat2, Matrix& F, Matrix& F_plus, Matrix& F_minus, 
	 int crr_flag = 1);

/****************************************************************
 *����
 *  int       point_num   �ǡ�������                input
 *  VectorN   eigen_vec   M_hat��9�Ĥθ�ͭ��        input
 *                        9�٥��ȥ�γ�����
 *  Matrix    *eigen_mat  M_hat��9�Ĥθ�ͭ����      input
 *                        eigen_vec�γ����Ǥ��б�
 *  double    e_hat2      ���Υ�����٥�ο�����  input
 *  Matrix&   F           �����������ù���          output
 *  Matrix&   F_plus      F��ɸ���а�(+)            output
 *  Matrix&   F_minus     F��   ��   (-)            output
 *  int       crr_flag    ��Ŭ�����¹ԥե饰        input
 *                        1�ΤȤ�F�κ�Ŭ�����¹�
 *                        �ǥե������: 1
 *
 ****************************************************************/


#endif //_pdp_h_
