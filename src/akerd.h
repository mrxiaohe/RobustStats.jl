#ifndef __AKERD__
#define __AKERD__
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
void akerd_C(double *pts, int npts, double *x, int n, double *alam, double hval, double *output) ;
void alam_C(double *fhat, int n, double gm, double aval, double *output) ;
#endif