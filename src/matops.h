#ifndef MATOPS_H
#define MATOPS_H

void inverse(double *a, double *b); // a = b^-1
void transpose(double *xc, double *x); // helper function for LAPACK row_major -> column major

#endif // !MATOPS_H
