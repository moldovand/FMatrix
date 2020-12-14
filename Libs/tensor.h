#ifndef _TENSOR_
#define _TENSOR_

//ijklテンソルクラス
class Tensor3333
{
  int in, jn, kn, ln;       //それぞれの要素数
  double**** t;             //要素

public:
  //1st constructor
  Tensor3333(int i = 3, int j = 3, int k = 3, int l = 3);

  //copy constructor
  Tensor3333(const Tensor3333&);

  //destructor
  ~Tensor3333();

  //要素取り出し
  double&       operator()(int, int, int, int);
  const double& operator()(int, int, int, int) const;


  //代入演算  
  const Tensor3333& operator=(const Tensor3333&);

  //算術演算
  Tensor3333 operator+(const Tensor3333&) const;  //足し算
  Tensor3333 operator-(const Tensor3333&) const;  //引き算
  Tensor3333 operator*(double) const;             //かけ算
  friend Tensor3333 operator*(double, const Tensor3333&);
  //3333テンソル*33行列
  Matrix operator*(const Matrix&);
  Tensor3333 operator/(double);  //割算  
  //符号換え
  friend Tensor3333 operator-(const Tensor3333&);
  //テンソル積(引数の２つの行列は３×３に限る)
  friend Tensor3333 TensorProduct(const Matrix&, const Matrix&);

  //それぞれの次元
  int DimI() const { return in; }
  int DimJ() const { return jn; }
  int DimK() const { return kn; }
  int DimL() const { return ln; }

  //初期化
  int Init();

  //複製による3333テンソルの要素生成
  void Fill();

  //複製による(33)(33)テンソルの要素生成
  void FillSym();

  //(33)(33)テンソルの全ての固有値と全ての固有行列を求める
  int EigenSym(VectorN&, Matrix*);

  //3333テンソルの全ての固有値と全ての固有行列を求める
  int Eigen(VectorN&, Matrix*);
                                                 
  //(33)(33)テンソルの一般逆を求める
  friend Tensor3333 GeneralInverseSym(const Tensor3333&, int times = 0);

  //3333テンソルの一般逆を求める
  friend Tensor3333 GeneralInverse(const Tensor3333&, int times = 0);
  
  //要素出力
  void Print();
  void Fprint(FILE*);
};

#endif //_TENSOR_
