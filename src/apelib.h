#ifndef APELIB_H
#define APELIB_H

#include <stdbool.h>
// #include <stdint.h>
// #include "sizes.h"

void newton_reunit(double *m);
void ape_smear(double *restrict V, double *restrict U, const double alpha, const bool sdir[D]);
void gram_schmidt(double *U);
void gram_schmidt_per(double *link);

#endif // !APELIB_H
