#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <lapacke.h>
#include "sizes.h"

void inverse(double *a, double *b) {

    memcpy(a,b,X*sizeof(double));

    int *IPIV = malloc(NC*sizeof(int));

    // lapack_int info_rf, info_ri;
    //
    // info_rf = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, NC, NC, (lapack_complex_double*)a, NC, IPIV);
    // info_ri = LAPACKE_zgetri(LAPACK_ROW_MAJOR, NC, (lapack_complex_double*)a, NC, IPIV);
    LAPACKE_zgetrf(LAPACK_ROW_MAJOR, NC, NC, (lapack_complex_double*)a, NC, IPIV);
    LAPACKE_zgetri(LAPACK_ROW_MAJOR, NC, (lapack_complex_double*)a, NC, IPIV);

}
