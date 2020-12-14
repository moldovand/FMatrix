#ifndef _pdp_h_
#define _pdp_h_

/*
 *3333テンソルM_hatの9個の固有値VectorN(9)と対応する固有行列の
 *正規直交系Matrix*, 及び二乗ノイズレベルの推定値e_hat2を与えて、
 *補正した基礎行列F, その標準偏位F_plus, F_minusを返す。
 *  ※(PDP = primary deviation pair)
 */

void PDP(int point_num, VectorN eigen_vec, Matrix* eigen_mat, 
	 double e_hat2, Matrix& F, Matrix& F_plus, Matrix& F_minus, 
	 int crr_flag = 1);

/****************************************************************
 *引数
 *  int       point_num   データ点数                input
 *  VectorN   eigen_vec   M_hatの9個の固有値        input
 *                        9ベクトルの各要素
 *  Matrix    *eigen_mat  M_hatの9個の固有行列      input
 *                        eigen_vecの各要素に対応
 *  double    e_hat2      二乗ノイズレベルの推定値  input
 *  Matrix&   F           補正した基礎行列          output
 *  Matrix&   F_plus      Fの標準偏位(+)            output
 *  Matrix&   F_minus     Fの   〃   (-)            output
 *  int       crr_flag    最適補正実行フラグ        input
 *                        1のときFの最適補正実行
 *                        デフォルト値: 1
 *
 ****************************************************************/


#endif //_pdp_h_
