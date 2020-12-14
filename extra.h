#ifndef _extra_h_
#define _extra_h_

//3333テンソルのトレースを計算する
double Trace(Tensor3333 T);

//射影テンソルを計算する
Tensor3333 ProjTensor(const Matrix& M);

//疑似射影テンソルを計算する
Tensor3333 Pseudo(Matrix M);

#endif //_extra_h_
