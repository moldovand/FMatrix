#ifndef _ImagePosition_h_
#define _ImagePosition_h_

//
//     Class  ImagePosition (≤Ë¡¸∞Ã√÷)
//

class ImagePosition
{
private:
  Vector _x;  // ≤Ë¡¸∞Ã√÷ ( Z-Vector )
  Matrix _V0; // ∂¶ ¨ª∂π‘ŒÛ
  double _f;  // æ«≈¿µ˜Œ•

public:
  // constructor
  ImagePosition();
  ImagePosition(double x, double y, double f);
  ImagePosition(const Vector& x, double f);
  ImagePosition(const Vector& x, const Matrix& V, double f);

  // copy constructor
  ImagePosition(const ImagePosition&);

  // destructor
  ~ImagePosition(){};

  // assingment operator
  ImagePosition& operator=(const ImagePosition&);

  const Vector&  Pos()  const  { return _x;  }
  const Matrix&  Cov()  const  { return _V0; }
  const double&  Foc()  const  { return _f;  }

  void SetPos(double x, double y);
  void SetPos(const Vector& x);
  void SetCov(double dx, double dy, double dxy);
  void SetCov(const MatrixMN& S);
  void SetCov(const Matrix&   V)  { _V0 = V; }
  void SetFoc(double f) { _f = f; }

  void ChangeFoc(double f);

  int  Print(FILE *fp = stdout, char *str = "#ImagePosition") const;
  int  LongPrint(FILE *fp = stdout, char *str = "#ImagePosition") const;
  int  Scan(FILE *fp = stdin);
};

#endif  // _ImagePosition_h_

