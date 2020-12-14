#include "statcomp.h"

//
//     Class ImagePosition
//

// constructor
ImagePosition::ImagePosition()
{
  _x  = Vector( 0.0, 0.0, 1.0 );
  _V0 = Pk;
  _f  = 600.0;
}

ImagePosition::ImagePosition(double x, double y, double f=600.0)
{
  _x  = Vector( x/f, y/f, 1.0 );
  _V0 = Pk;
  _f  = f;
}

ImagePosition::ImagePosition(const Vector& x, double f=600.0)
{
  _x  = x;
  _V0 = Pk;
  _f  = f;
}

ImagePosition::ImagePosition(const Vector& x, const Matrix& V, double f=600.0)
{
  _x  = x;
  _V0 = V;
  _f  = f;

}

// copy constructor
ImagePosition::ImagePosition(const ImagePosition& ipos)
{
  _x  = ipos._x;
  _V0 = ipos._V0;
  _f  = ipos._f;
}

// assignment operator
ImagePosition& ImagePosition::operator=(const ImagePosition &ipos)
{
  if (this != &ipos) {
    _x  = ipos._x;
    _V0 = ipos._V0;
    _f  = ipos._f;
  }

  return *this;
}

void ImagePosition::SetPos(double x, double y)
{
  _x  = Vector( x/_f, y/_f, 1.0 );
}

void ImagePosition::SetPos(const Vector& x)
{
  if ( x[2] != 1.0 ) {
    fprintf(stderr, "Argument isn't Z-Vector. (Error in ImagePosition) \n");
    return;
  }
  else
    _x = x;
}

void ImagePosition::SetCov(double dx, double dy = -1.0, double gxy=0.0)
{
  if ( dy == -1.0 )  dy = dx;
  _V0 = Matrix( sqr(dx),  gxy,      0.0,
		gxy,      sqr(dy),  0.0,
		0.0,      0.0,      0.0  );
  _V0 /= sqr(_f);
}

void ImagePosition::SetCov(const MatrixMN& S)
{
  if ( ( S.DimI() != 2 ) || ( S.DimJ() != 2 ) ) {
    fprintf(stderr, "Argument isn't 22-Matrix. (Error in ImagePosition) \n");
    return;
  }
  else {
    _V0 = Matrix( S[0][0],  S[0][1],  0.0,
		  S[1][0],  S[1][1],  0.0,
		  0.0,      0.0,      0.0  );
    _V0 /= sqr(_f);
  }
}

void ImagePosition::ChangeFoc(double f)
{
  double tmp = _f / f;

  _x  = Matrix( tmp, tmp, 1.0) * _x;
  _V0 = sqr(tmp) * _V0;
  _f  = f;
}

// 画像点の出力
int ImagePosition::Print(FILE* fp, char *str) const
{
  fprintf(fp, "%s\n", str);
  _x.Print(fp, " Pos = ");
  _V0.Print(fp, " Cov = ");
  fprintf(fp, " Foc = %13.6f \n", _f);
  fprintf(fp, "%%\n");

  return true;
}

int ImagePosition::LongPrint(FILE* fp, char *str) const
{
  fprintf(fp, "%s\n", str);
  _x.LongPrint(fp, " Pos = ");
  _V0.LongPrint(fp, " Cov = ");
  fprintf(fp, " Foc = %22.16f \n", _f);
  fprintf(fp, "%%\n");

  return true;
}

// 画像点の入力
int ImagePosition::Scan(FILE* fp)
{
  char info[15];

  if (ScanForScan(fp) != IDETA) {
    fprintf(stderr, "Format error, ImagePosition.\n");
    return false;
  }
  fscanf(fp, " %s ", info);
  if (strcmp(info, "#ImagePosition") != 0) {
    fprintf(stderr, "Format error, ImagePosition.\n");
    return false;
  }

  _x.Scan(fp);
  _V0.Scan(fp);

  fscanf(fp, " %s ", info);
  if (strcmp(info, "Foc") != 0) {
    fprintf(stderr, "Format error, ImagePosition.\n");
    return false;
  }
  if (ScanForScan(fp) == DIGIT) fscanf(fp, "%lf", &_f);
  else {
    fprintf(stderr, "Scan error, in ImagePosition.\n");
    return false;
  }

  return true;
}
