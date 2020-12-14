#ifndef _filedeal_h_
#define _filedeal_h_

//ファイルポインタを次の行の先頭に置く
void NextLine(FILE* fp);

//ファイルポインタを次の数の先頭に置く
void NextNumber(FILE* fp);

//ポイントデータをベクトルに格納
void Text2(char* fname, double**& vec, const int& pnum, 
	   double& xm, double& ym, double& Lx, double& Ly);

#endif //_filedeal_h_
