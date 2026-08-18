#ifndef APELIB_H
#define APELIB_H

void ape_smear(double *V, double *U, double alpha, long int sdir[4]);
void ape_smear_fast(double *V, double *U, double alpha, long int sdir[4]);
void gram_schmidt(double *U);
void gram_schmidt_per(double *link);
void newton_reunit(double *m);

#endif // !APELIB_H
