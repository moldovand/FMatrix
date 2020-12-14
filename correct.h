#ifndef _Correct_h_
#define _Correct_h_

/*
 *行列Fとその正規化共分散テンソルV_0Fを与えて、
 *FがdetF=0となるように最適に補正する
 */
Matrix Correct(Matrix F, Tensor3333 V_0F, double eps_d = 1.0e-12);

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


#endif //_Correct_h_
