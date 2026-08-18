#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#define LAPACK_ILP64
#include <openblas64/lapack.h>
// #include <lapacke.h>
#include "sizes.h"
#include "matops.h"

void inverse(double *a, double *b) {

    // Better to use LAPACK directly. Mainly for portability

    transpose(a,b);
    // memcpy(a,b,X*sizeof(double));
    
    lapack_int *IPIV = malloc(NC*sizeof(int32_t));
    lapack_int info_rf, info_ri;
    lapack_complex_double work[16];
    lapack_int lwork = 100*NC;
    //
    // info_rf = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, NC, NC, (lapack_complex_double*)a, NC, IPIV);
    // info_ri = LAPACKE_zgetri(LAPACK_ROW_MAJOR, NC, (lapack_complex_double*)a, NC, IPIV);
    // LAPACKE_zgetrf(LAPACK_ROW_MAJOR, NC, NC, (lapack_complex_double*)a, NC, IPIV);
    // LAPACKE_zgetri(LAPACK_ROW_MAJOR, NC, (lapack_complex_double*)a, NC, IPIV);

    LAPACK_zgetrf(&NC, &NC, (lapack_complex_double*)a, &NC, IPIV, &info_rf);
    LAPACK_zgetri(&NC, (lapack_complex_double*)a, &NC, IPIV, work, &lwork, &info_ri);

    transpose(a,a);

    free(IPIV);

}

void transpose(double *xc, double *x) {
    double y[X];
    for (int64_t i = 0; i < X; i++) {
        y[i] = x[i];
    }
    xc[0] = y[0];       
    xc[1] = y[1];      
    xc[2] = y[6];       
    xc[3] = y[7];      
    xc[4] = y[12];      
    xc[5] = y[13];     
    xc[6] = y[2];       
    xc[7] = y[3];      
    xc[8] = y[8];       
    xc[9] = y[9];
    xc[10] = y[14];
    xc[11] = y[15];
    xc[12] = y[4];
    xc[13] = y[5];
    xc[14] = y[10];
    xc[15] = y[11];
    xc[16] = y[16];
    xc[17] = y[17];
}
