#ifndef _Matrix_h_
#define _Matrix_h_

#include "MatrixMN.h"

//
//     Class  Matrix (3-Dimensional Square Matrix)
//

class Matrix : public MatrixMN
{ 
private:

public:
  // constructor
  Matrix();
  Matrix(const Vector&, const Vector&, const Vector&);
  Matrix(const Vector&, double);
  Matrix(double*);
  Matrix(double**);
  Matrix(double, double, double,
	 double, double, double,
	 double, double, double);
  Matrix(const Vector&, const Vector&);
  Matrix(double, double, double);

  // copy constructor
  Matrix(const Matrix&);

  // destructor
  ~Matrix();

  // assignment operator
  Matrix& operator=(const Matrix&);
  
  double* operator[](int i);                       // ¬Ëiπ‘
  const double* operator[](int i) const;
  Vector  operator()(int j) const;                 // ¬ËjŒÛ

  int     Print(FILE *fp = stdout, char *str = "") const;
  int     LongPrint(FILE *fp = stdout, char *str = "") const;
  int     Scan(FILE *fp = stdin);

  // ¬–æŒπ‘ŒÛ?
  int Symmetric() const;
  // ΩÈ¥¸≤Ω
  int     Init();
  // √±∞Ã≤Ω
  int     Unit();

  // »Ê≥”
  friend int operator==(const Matrix&, const Matrix&);
  friend int operator!=(const Matrix&, const Matrix&);

  Matrix& operator+=(const Matrix&);
  Matrix& operator-=(const Matrix&);
  Matrix& operator*=(double);
  Matrix& operator/=(double);

  // º°∏µ§Ú4, 4 º°∏µ§À¡˝§‰§π
  friend MatrixMN operator+(const Matrix&, double);
  friend MatrixMN operator+(double, const Matrix&);
  
  // ªªΩ—±Èªª
  friend Vector  operator*(const Matrix&, const Vector&);
  friend Matrix  operator*(const Matrix&, const Matrix&);
  friend Matrix  operator*(double, const Matrix&);
  friend Matrix  operator*(const Matrix&, double);
  friend Matrix  operator%(const Vector&, const Matrix&);  // ≥∞¿—
  friend Matrix  operator%(const Matrix&, const Vector&);  // ≥∞¿—
  friend Matrix  operator%(const Matrix&, const Matrix&);  // ≥∞¿—
  friend Matrix  operator/(const Matrix&, double);
  friend Matrix  operator+(const Matrix&, const Matrix&);
  friend Matrix  operator-(const Matrix&);         
  friend Matrix  operator-(const Matrix&, const Matrix&);

  // ºÕ±∆π‘ŒÛ
  friend Matrix  Projection(const Vector&);
  // ∆‚¿—
  friend double  operator,(const Matrix&, const Matrix&);      
  // •Œ•Î•‡
  friend double  Norm(const Matrix&);
  // ¿µµ¨≤Ω
  friend Matrix  Normalize(const Matrix&);
  // •È•Û•Ø
  friend int     Rank(const Matrix&);
  // ≈æ√÷
  friend Matrix  Trans(const Matrix&);
  // •»•Ï°º•π
  friend double  Trace(const Matrix&);
  // µ’π‘ŒÛ
  friend Matrix  Inverse(const Matrix&);
  // ∞Ï»Ãµ’π‘ŒÛ
  friend Matrix  GeneralInverse(const Matrix&, int times = 0);
  // π‘ŒÛº∞
  friend double  Det(const Matrix&);
  // Õæ∞¯ª“π‘ŒÛ
  friend Matrix  Cofactor(const Matrix&);
  // ∏«Õ≠√Õ°¶∏«Õ≠•Ÿ•Ø•»•Î
  int    Householder(Vector&, Matrix&) const;
  // ∆√∞€√Õ ¨≤Ú £À°·£÷£Ã£’t
//  int    Svdecomp(Matrix&, Matrix&, Matrix&) const;
  // ≤Û≈æº¥§»≤Û≈æ≥—(≈Ÿ)§À ¨≤Ú§π§Î
  int    Rtrans(Vector&, double&) const;
};

extern const Matrix Im;
extern const Matrix Zm;
extern const Matrix Pk;

#endif  // _Matrix_h_
