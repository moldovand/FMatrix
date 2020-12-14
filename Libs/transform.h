#ifndef _TRANSFORM_
#define _TRANSFORM_

//3333テンソル->99行列
MatrixMN Type99(const Tensor3333&);

//(33)(33)テンソル->66行列
MatrixMN Type66(const Tensor3333&);

//6ベクトル->(33)行列
Matrix Type33sym(const VectorN&);

//66行列->(33)(33)テンソル
Tensor3333 Type3333sym(const MatrixMN&);

//9ベクトル->33行列
Matrix Type33(const VectorN&);

//33行列->9ベクトル
VectorN Type9(const Matrix&);

//99行列->3333テンソル
Tensor3333 Type3333(const MatrixMN&);

#endif //_TRANSFORM_
