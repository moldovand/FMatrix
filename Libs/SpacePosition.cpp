#include "statcomp.h"

//
//     Class SpacePosition
//

// constructor
SpacePosition::SpacePosition()
{
  _pos = Zv;
  _cov = Im;
}

SpacePosition::SpacePosition(double x, double y, double z)
{
  _pos = Vector( x, y, z);
  _cov = Im;
}

SpacePosition::SpacePosition(double  x, double  y, double z,
			     double sx, double sy, double sz)
{
  _pos = Vector( x, y, z);
  _cov = Matrix( sqr(sx), 0.0,     0.0,
		 0.0,     sqr(sy), 0.0,
		 0.0,     0.0,     sqr(sz) );
}

SpacePosition::SpacePosition(double x,   double y,   double z,
			     double sx,  double sy,  double sz,
			     double gxy, double gyz, double gzx)
{
  _pos = Vector(x, y, z);
  _cov = Matrix( sqr(sx), gxy,     gzx,
		 gxy,     sqr(sy), gyz,
		 gzx,     gyz,     sqr(sz) );
}

// copy constructor
SpacePosition::SpacePosition(const SpacePosition& p)
{
  _pos = p._pos;
  _cov = p._cov;
}

// assignment operator
SpacePosition& SpacePosition::operator=(const SpacePosition &p)
{
  if (this != &p) {
    _pos = p._pos;
    _cov = p._cov;
  }

  return *this;
}

void SpacePosition::SetPos(double x, double y, double z)
{
  _pos = Vector(x, y, z);
}

void SpacePosition::SetCov(double sx, double sy, double sz)
{
  _cov = Matrix( sqr(sx), 0.0,     0.0,
		 0.0,     sqr(sy), 0.0,
		 0.0,     0.0,     sqr(sz) );
}

void SpacePosition::SetCov(double  sx, double  sy, double  sz,
			   double gxy, double gyz, double gzx)
{
  _cov = Matrix( sqr(sx), gxy,     gzx,
		 gxy,     sqr(sy), gyz,
		 gzx,     gyz,     sqr(sz) );
}

// 画像点の出力
int SpacePosition::Print(FILE* fp, char *str) const
{
  fprintf(fp, "%s\n", str);
  _pos.Print(fp, " Pos = ");
  _cov.Print(fp, " Cov = ");
  fprintf(fp, "%%\n");

  return true;
}

int SpacePosition::LongPrint(FILE* fp, char *str) const
{
  fprintf(fp, "%s\n", str);
  _pos.LongPrint(fp, " Pos = ");
  _cov.LongPrint(fp, " Cov = ");
  fprintf(fp, "%%\n");

  return true;
}

// 画像点の入力
int SpacePosition::Scan(FILE* fp)
{
  char info[15];

  if (ScanForScan(fp) != IDETA) {
    fprintf(stderr, "Format error, SpacePosition.\n");
    return false;
  }
  fscanf(fp, " %s ", info);
  if (strcmp(info, "#SpacePosition") != 0) {
    fprintf(stderr, "Format error, SpacePosition.\n");
    return false;
  }
  
  _pos.Scan(fp);
  _cov.Scan(fp); 

  return true;
}
