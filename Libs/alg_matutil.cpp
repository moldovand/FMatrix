/***********************************************************
	alg_matutil.c -- 行列
***********************************************************/
/* 行列操作の小道具集 */
#include <stdio.h>
#include <stdlib.h>
#include "alg_matutil.h"

void _alg_error(char *message)
{
  fprintf(stderr, "\n%s\n", message);  exit(EXIT_FAILURE);
}

vector _alg_newvec(int n)
{
  return (SCALAR *)malloc(sizeof(SCALAR) * n);
}

matrix _alg_newmat(int nrow, int ncol)
{
  int i;
  matrix a;
  a = (SCALAR **)malloc((nrow + 1) * sizeof(void *));
  if (a == NULL) return NULL;  /* 記憶領域不足 */
  for (i = 0; i < nrow; i++) {
    a[i] = (SCALAR *)malloc(sizeof(SCALAR) * ncol);
    if (a[i] == NULL) {
      while (--i >= 0) free(a[i]);
      free(a);  return NULL;  /* 記憶領域不足 */
    }
  }
  a[nrow] = NULL;  /* 行の数を自動判断するための工夫 */
  return a;
}

vector _alg_new_vector(int n)
{
  vector v;
  v = _alg_newvec(n);
  if (v == NULL) _alg_error("記憶領域不足.");
  return v;
}

matrix _alg_new_matrix(int nrow, int ncol)
{
  matrix a;
  a = _alg_newmat(nrow, ncol);
  if (a == NULL) _alg_error("記憶領域不足.");
  return a;
}

void _alg_free_vector(vector v)
{
  free(v);
}

void _alg_free_matrix(matrix a)
{
  matrix b;
  b = a;
  while (*b != NULL) free(*b++);
  free(a);
}

double _alg_innerproduct(int n, vector u, vector v)
{
  int i, n5;
  double s;
  s = 0;  n5 = n % 5;
  for (i = 0; i < n5; i++) s += u[i]*v[i];
  for (i = n5; i < n; i += 5)
    s += u[i]*v[i] + u[i+1]*v[i+1] + u[i+2]*v[i+2]
      + u[i+3]*v[i+3] + u[i+4]*v[i+4];
  return s;
}
