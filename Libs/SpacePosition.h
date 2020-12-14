#ifndef _SpacePosition_h_
#define _SpacePosition_h_

//
//     Class  SpacePosition
//


class SpacePosition
{ 
private:
  Vector _pos;
  Matrix _cov;

public:
  // constructor
  SpacePosition();
  SpacePosition(double x, double y, double z);
  SpacePosition(double x, double y, double z, 
		double sx,  double sy,  double sz);
  SpacePosition(double x, double y, double z,
		double sx,  double sy,  double sz,
		double gxy, double gyz, double gzx);

  // copy constructor
  SpacePosition(const SpacePosition&);

  // destructor
  ~SpacePosition(){};

  // assingment operator
  SpacePosition& operator=(const SpacePosition&);

  void SetPos(double x, double y, double z);
  void SetCov() { _cov = Im; }
  void SetCov(double sx, double sy, double sz);
  void SetCov(double  sx, double  sy, double  sz,
	      double gxy, double gyz, double gzx);

  const Vector&  Pos()  const { return _pos; }
  const Matrix&  Cov()  const { return _cov; }

  int     Print(FILE *fp = stdout, char *str = "#SpacePosition") const;
  int     LongPrint(FILE *fp = stdout, char *str = "#SpacePosition") const;
  int     Scan(FILE *fp = stdin);


  // friend function, friend class
  friend int _SPosSaveGI(char *filename, SpacePosition *spos,
			 char *type, int n);
  friend SpacePosition* _SPosLoadGI(char *filename, int n);
};

#endif  // _SpacePosition_h_

