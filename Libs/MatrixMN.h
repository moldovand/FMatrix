#ifndef _MatrixMN_h_
#define _MatrixMN_h_

#include "VectorN.h"
#include "Vector.h"
//#include "Matrix.h"

//
//     Class  MatrixMN (MN-Matrix)
//

class MatrixMN     
{ 
private:
protected:
  double **m;
  int      dim_i, dim_j;

public:
  // constructor
  MatrixMN(int di = 3, int dj = 3);
  MatrixMN(const VectorN&, const VectorN&, const VectorN&);
  MatrixMN(double*, int, int);
  MatrixMN(double**, int, int);
  MatrixMN(const VectorN&, const VectorN&);
//  MatrixMN(const Matrix&, const Vector&, const Vector&, double);

  // copy constructor
  MatrixMN(const MatrixMN&);

  // destructor
  ~MatrixMN();
  
  // assignment operator
  MatrixMN& operator=(const MatrixMN&);
  
  double* operator[](int i);                       // ¬Ëiπ‘
  const double* operator[](int i) const;
  VectorN operator()(int j) const;                 // ¬ËjŒÛ
  int     DimI() const { return dim_i; }
  int     DimJ() const { return dim_j; }
  
  int     Print(FILE *fp = stdout, char *str = "") const;
  int     LongPrint(FILE *fp = stdout, char *str = "") const;
  int     Scan(FILE *fp = stdin);

  // ¬–æŒπ‘ŒÛ?
  int Symmetric() const;
  // ΩÈ¥¸≤Ω
  int     Init(int i = 0, int j = 0);
  // √±∞Ã≤Ω
  int     Unit();
  
  // »Ê≥”
  friend int operator==(const MatrixMN&, const MatrixMN&);
  friend int operator!=(const MatrixMN&, const MatrixMN&);
  
  MatrixMN& operator+=(const MatrixMN&);
  MatrixMN& operator-=(const MatrixMN&);
  MatrixMN& operator*=(double);
  MatrixMN& operator/=(double);
  
  // º°∏µ§Údim_i + 1, dim_j + 1 º°∏µ§À¡˝§‰§π
  friend MatrixMN  operator+(const MatrixMN&, double);
  friend MatrixMN  operator+(double, const MatrixMN&);
  
  // ªªΩ—±Èªª
  friend VectorN   operator*(const MatrixMN&, const VectorN&);
  friend MatrixMN  operator*(const MatrixMN&, const MatrixMN&);
  friend MatrixMN  operator*(double, const MatrixMN&);
  friend MatrixMN  operator*(const MatrixMN&, double);
  friend MatrixMN  operator%(const VectorN&, const MatrixMN&);
  friend MatrixMN  operator/(const MatrixMN&, double);
  friend MatrixMN  operator+(const MatrixMN&, const MatrixMN&);
  friend MatrixMN  operator-(const MatrixMN&);         
  friend MatrixMN  operator-(const MatrixMN&, const MatrixMN&);

  // ºÕ±∆π‘ŒÛ
  friend MatrixMN  Projection(const VectorN&);
  // ∆‚¿—
  friend double    operator,(const MatrixMN&, const MatrixMN&);      
  // •Œ•Î•‡
  friend double    Norm(const MatrixMN&);
  // ¿µµ¨≤Ω
  friend MatrixMN  Normalize(const MatrixMN&);
  // •È•Û•Ø
  friend int       Rank(const MatrixMN&);
  // ≈æ√÷
  friend MatrixMN  Trans(const MatrixMN&);
  // •»•Ï°º•π
  friend double    Trace(const MatrixMN&);
  // µ’π‘ŒÛ
  friend MatrixMN  Inverse(const MatrixMN&);
  // ∞Ï»Ãµ’π‘ŒÛ
  friend MatrixMN  GeneralInverse(const MatrixMN&, int times = 0);
  // π‘ŒÛº∞
  friend double    Det(const MatrixMN&);
  // Õæ∞¯ª“π‘ŒÛ
//  friend MatrixMN  Cofactor(const MatrixMN&);
  // ∏«Õ≠√Õ°¶∏«Õ≠•Ÿ•Ø•»•Î
  int    Householder(VectorN&, MatrixMN&) const;
  // ∆√∞€√Õ ¨≤Ú £À°·£÷£Ã£’t
//  int    Svdecomp(MatrixMN&, MatrixMN&, MatrixMN&) const;
};

int EddingtonEpsilon(int, int, int);
int KroneckerDelta(int, int);
     
#endif  // _MatrixMN_h_


