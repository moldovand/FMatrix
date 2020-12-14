#include <statcomp.h>
#include <ctype.h>


//ファイルポインタを次の行の先頭に置く
void NextLine(FILE* fp)
{
  while (!feof(fp))
    if (fgetc(fp) == '\n') break;
}


//ファイルポインタを次の数の先頭に置く
void NextNumber(FILE* fp)
{
  char c = '\0';
  while (!feof(fp))
    {
      c = fgetc(fp);
      if (isdigit(c) || c == '-') break;
    }
  
  ungetc(c, fp);
}


//ポイントデータを配列に格納
void Text2(char* fname, double**& vec, const int& pnum, 
	   double& xm, double& ym, double& Lx, double& Ly)
{
  int i, j;
  FILE* fp = fopen(fname, "r");
  NextNumber(fp);
    
  //ポイント数を読み込む
  fscanf(fp, "%d", &pnum);
  vec = new double*[pnum];
  for (i = 0; i < pnum; ++i)
    vec[i] = new double[2];
  for (i = 0; i < pnum; ++i)
    for (j = 0; j < 2; ++j) vec[i][j] = 0.0;
  NextNumber(fp);
    
  //読み込んだ座標値を配列に格納しながら座標値の最大値・最小値を選択
  double x = 0.0, y = 0.0;
  double sumx = 0.0, sumy = 0.0;
  double xmax = -1.0e+8, xmin = 1.0e+8;
  double ymax = -1.0e+8, ymin = 1.0e+8;
  for (i = 0; i < pnum; i++)
    {
      fscanf(fp, "%lf", &x);
      NextNumber(fp);
      fscanf(fp, "%lf", &y);

      vec[i][0] = x;
      vec[i][1] = y;
      sumx += x; sumy += y;
      if (xmax < x) xmax = x;
      if (xmin > x) xmin = x;
      if (ymax < y) ymax = y;
      if (ymin > y) ymin = y;

      NextNumber(fp);
    }

  //座標値の平均値などを算出
  xm = sumx/(double)pnum;
  ym = sumy/(double)pnum;
  Lx = xmax - xmin;
  Ly = ymax - ymin;

  fclose(fp);
}
