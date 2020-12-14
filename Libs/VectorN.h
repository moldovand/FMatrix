#ifndef _VectorN_h_
#define _VectorN_h_

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include "FileScan.h"
#include "defVM.h"

//
//     Class  VectorN (N-Dimension Vector)
//

class VectorN     
{ 
private:
protected:
  double *v;
  int     dim;

public:
  // constructor
  VectorN(int d = 3);                                        
  VectorN(double*, int d);

  // copy constructor
  VectorN(const VectorN&);

  // destructor
  ~VectorN();

  // assingment operator
  VectorN& operator=(const VectorN&);
  
  double&       operator[](int i)       { return v[i]; }
  const double& operator[](int i) const { return v[i]; }
  int     Dim() const { return dim; }
  
  int     Print(FILE *fp = stdout, char *str = "") const;
  int     LongPrint(FILE *fp = stdout, char *str = "") const;
  int     Scan(FILE *fp = stdin);

  // ½é´ü²½
  int     Init(int d = 0);
  
  // Èæ³Ó
  friend int operator==(const VectorN&, const VectorN&);     
  friend int operator!=(const VectorN&, const VectorN&);

  VectorN& operator+=(const VectorN&);
  VectorN& operator-=(const VectorN&);
  VectorN& operator*=(double);
  VectorN& operator/=(double);

  // ¼¡¸µ¤ò1¼¡¸µÁý¤ä¤¹
  friend VectorN  operator+(const VectorN&, double);
  friend VectorN  operator+(double, const VectorN&);
  
  // »»½Ñ±é»»
  friend VectorN  operator-(const VectorN&);
  friend VectorN  operator-(const VectorN&, const VectorN&);
  friend VectorN  operator+(const VectorN&, const VectorN&);
  friend VectorN  operator*(double, const VectorN&);
  friend VectorN  operator*(const VectorN&, double);
  friend VectorN  operator%(const VectorN&, const VectorN&); // ³°ÀÑ
  friend VectorN  operator/(const VectorN&, double);
  // ÆâÀÑ
  friend double   operator,(const VectorN&, const VectorN&);
  // ¥Î¥ë¥à
  friend double   Norm(const VectorN&);
  // Àµµ¬²½
  friend VectorN  Normalize(const VectorN&);                       
  // Æ©»ëºîÍÑÁÇ
  friend VectorN  Proj(const VectorN&);
};

double sqr(double);
int    sgn(double);
double cbrt(double);

#endif  // _VectorN_h_


