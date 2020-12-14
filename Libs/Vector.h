#ifndef _Vector_h_
#define _Vector_h_

#include "VectorN.h"

//
//     Class  Vector (3-Dimension Vector)
//

class Vector : public VectorN 
{ 
private:
  
public:
  // constructor
  Vector();                                        
  Vector(double, double, double);
  Vector(double*);

  // copy constructor
  Vector(const Vector&);

  // destructor
  ~Vector();

  // assingment operator
  Vector& operator=(const Vector&);
  
  double&       operator[](int i)       { return v[i]; }
  const double& operator[](int i) const { return v[i]; }

  int     Print(FILE *fp = stdout, char *str = "") const;
  int     LongPrint(FILE *fp = stdout, char *str = "") const;
  int     Scan(FILE *fp = stdin);

  // Èæ³Ó
  friend int operator==(const Vector&, const Vector&);     
  friend int operator!=(const Vector&, const Vector&);

  Vector& operator+=(const Vector&);
  Vector& operator-=(const Vector&);
  Vector& operator*=(const double);
  Vector& operator/=(const double);

  // ¼¡¸µ¤ò1¼¡¸µÁý¤ä¤¹
  friend VectorN operator+(const Vector&, double);
  friend VectorN operator+(double, const Vector&);

  // »»½Ñ±é»»
  friend Vector  operator-(const Vector&);
  friend Vector  operator-(const Vector&, const Vector&);
  friend Vector  operator+(const Vector&, const Vector&);
  friend Vector  operator*(double, const Vector&);
  friend Vector  operator*(const Vector&, double);
  friend Vector  operator%(const Vector&, const Vector&); // ³°ÀÑ
  friend Vector  operator/(const Vector&, double);
  // ÆâÀÑ
  friend double  operator,(const Vector&, const Vector&);
  // ¥¹¥«¥é£³½ÅÀÑ
  friend double  ScalarTP(const Vector&, const Vector&, const Vector&);
  // ¥Î¥ë¥à
  friend double  Norm(const Vector&);
  // Àµµ¬²½
  friend Vector  Normalize(const Vector&);                       
  // Æ©»ëºîÍÑÁÇ
  friend Vector  Proj(const Vector&);
};

extern const Vector Zv;
extern const Vector Iv;
extern const Vector Jv;
extern const Vector Kv;

#endif  // _Vector_h_

