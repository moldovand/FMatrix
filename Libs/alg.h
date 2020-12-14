#ifndef _ALG_H_
#define _ALG_H_

/***********************************************************/
/* alg.h -- Header file for functions                      */
/*          from the book "Algorithm Encyclopedia"         */
/***********************************************************/
#define SCALAR double

double _alg_gauss(int n, SCALAR **a, SCALAR *b, SCALAR *x); /* Not used */
double _alg_matinv(int n, SCALAR **a, SCALAR **a_inv);
double _alg_det(int n, SCALAR **a);
int    _alg_eigen(int n, SCALAR **a, SCALAR *d);

#undef SCALAR

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS
#endif

#ifndef EXIT_FAILURE
#define EXIT_FAILURE
#endif


#endif  /* _ALG_H */
