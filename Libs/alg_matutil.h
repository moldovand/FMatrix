/***********************************************************
	alg_matutil.h -- 行列 (ヘッダファイル)
***********************************************************/
#ifndef MATUTIL
#define MATUTIL

#define SCALAR double
typedef SCALAR *vector, **matrix;

void   _alg_error(char *message);
vector _alg_newvec(int n);
matrix _alg_newmat(int nrow, int ncol);
vector _alg_new_vector(int n);
matrix _alg_new_matrix(int nrow, int ncol);
void   _alg_free_vector(vector v);
void   _alg_free_matrix(matrix a);
double _alg_innerproduct(int n, vector u, vector v);

#endif  /* MATUTIL */
