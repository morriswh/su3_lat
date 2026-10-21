#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <sys/resource.h>
#include <inttypes.h>
// #include <lapacke.h>
#include <complex.h>
// #include <gmp.h>
// #include <mpfr.h>
#include <math.h>
#include <time.h>
#include "sizes.h"
#include "matops.h"
#include "su3lib.h"
#include "apelib.h"

static inline void inv_adjoint(double *restrict z, const double *restrict x) {
    // x^{-dagger} = conj(cof(x)) / conj(det x), no LAPACK needed for 3x3
    const double complex *a = (const double complex *)x;
    double complex c[9];
    // cofactor rows: r0 = r1 x r2, r1 = r2 x r0, r2 = r0 x r1
    c[0] = a[4]*a[8] - a[5]*a[7];
    c[1] = a[5]*a[6] - a[3]*a[8];
    c[2] = a[3]*a[7] - a[4]*a[6];
    c[3] = a[7]*a[2] - a[8]*a[1];
    c[4] = a[8]*a[0] - a[6]*a[2];
    c[5] = a[6]*a[1] - a[7]*a[0];
    c[6] = a[1]*a[5] - a[2]*a[4];
    c[7] = a[2]*a[3] - a[0]*a[5];
    c[8] = a[0]*a[4] - a[1]*a[3];
    const double complex d = a[0]*c[0] + a[1]*c[1] + a[2]*c[2];
    const double complex s = conj(1.0/d);
    double complex *o = (double complex *)z;
    for (int k = 0; k < 9; k++) {
        o[k] = conj(c[k])*s;
    }
}

void newton_reunit(double *m) {
    double detm[2];
    det(detm, m);
    const double complex s = 1.0/cpow(detm[0] + detm[1]*I, 1.0/3.0);

    double x[X], xai[X];
    double complex *xc = (double complex *)x;
    const double complex *mc = (const double complex *)m;
    for (int32_t i = 0; i < 9; i++) {
        xc[i] = s*mc[i];
    }

    for (int32_t it = 0; it < 100; it++) {
        inv_adjoint(xai, x);
        double delta = 0.0;

        for (int32_t i = 0; i < X; i++) {
            const double xn = 0.5*(x[i] + xai[i]);
            delta += fabs(xn - x[i]);
            x[i] = xn;
        }

        if (delta < 1e-9) {
            break;
        }
    }
    memcpy(m, x, X*sizeof(double));
}

void ape_smear(double *restrict V, double *restrict U, const double alpha, const bool sdir[D]) {
    // make coefficients
    int32_t n_mu = 0;
    for (int32_t mu = 0; mu < D; mu++) if(sdir[mu]) n_mu++;
    const double alpha_norm = alpha / (2.0 * (n_mu -1));
    const double one_minus_alpha = 1 - alpha;

    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < LVOLUME; i++) {
        for (int32_t mu = 0; mu < D; mu++) {
            if ( !sdir[mu] ) {
                memcpy(&V[X*(D*i + mu)], &U[X*(D*i + mu)], sizeof(double)*X);
                continue;
            }
            double staple[X], staple_sum[X];
            memset(staple_sum, 0, sizeof(staple_sum));
            for (int32_t nu = 0; nu < D; nu++) {
                if ( sdir[nu] && nu != mu ) {
                    staples(staple, U, i, mu, nu);
                    for (int32_t j = 0; j < X; j++) {
                        staple_sum[j] += staple[j];
                    }
                }
            }
            for (int32_t j = 0; j < X; j++) {
                V[X*(D*i + mu)+j] = one_minus_alpha * U[X*(D*i + mu)+j] + alpha_norm * staple_sum[j];
            }
            newton_reunit(&V[X*(D*i + mu)]);
        }
    }
}

// 
void gram_schmidt(double *U) {
    int64_t i;
    #pragma omp parallel for private(i)
    for (i = 0; i < LLENGTH; i++) {

        double temp[X];
        double dot_real, dot_imag;
        double norm_u2, norm_u;
        double norm_v2, norm_v;

        memcpy(temp,&U[X*i],sizeof(double)*X);

        int64_t p = X*i; // position in big array for config
        // get real and imaginary parts of dot product
        dot_real = temp[0]*temp[6] + temp[1]*temp[7]
                 + temp[2]*temp[8] + temp[3]*temp[9]
                 + temp[4]*temp[10] + temp[5]*temp[11];
        dot_imag = temp[0]*temp[7] - temp[1]*temp[6]
                 + temp[2]*temp[9] - temp[3]*temp[8]
                 + temp[4]*temp[11] - temp[5]*temp[10];

        // norm of u vector
        norm_u2 = temp[0]*temp[0] + temp[1]*temp[1]
                + temp[2]*temp[2] + temp[3]*temp[3]
                + temp[4]*temp[4] + temp[5]*temp[5];
        norm_u = sqrt(norm_u2);

        // construct orthogonal vector for second row
        temp[6] = temp[6] - (temp[0]*dot_real - temp[1]*dot_imag)/norm_u2; 
        temp[7] = temp[7] - (temp[0]*dot_imag + temp[1]*dot_real)/norm_u2; 
        temp[8] = temp[8] - (temp[2]*dot_real - temp[3]*dot_imag)/norm_u2; 
        temp[9] = temp[9] - (temp[2]*dot_imag + temp[3]*dot_real)/norm_u2; 
        temp[10] = temp[10] - (temp[4]*dot_real - temp[5]*dot_imag)/norm_u2; 
        temp[11] = temp[11] - (temp[4]*dot_imag + temp[5]*dot_real)/norm_u2; 

        // norm of v vector
        norm_v2 = temp[6]*temp[6] + temp[7]*temp[7]
                + temp[8]*temp[8] + temp[9]*temp[9]
                + temp[10]*temp[10] + temp[11]*temp[11];
        norm_v = sqrt(norm_v2);

        U[p] = temp[0]/norm_u; // u1.real
        U[p+1] = temp[1]/norm_u; // u1.imag
        U[p+2] = temp[2]/norm_u; // u2.real
        U[p+3] = temp[3]/norm_u; // u2.imag
        U[p+4] = temp[4]/norm_u; // u3.real
        U[p+5] = temp[5]/norm_u; // u3.imag
        U[p+6] = temp[6]/norm_v; // v1.real
        U[p+7] = temp[7]/norm_v; // v1.imag
        U[p+8] = temp[8]/norm_v;  // v2.real
        U[p+9] = temp[9]/norm_v;  // v2.imag
        U[p+10] = temp[10]/norm_v; // v3.real
        U[p+11] = temp[11]/norm_v; // v3.imag
        U[p+12] =   (U[p+2]*U[p+10] - U[p+3]*U[p+11])
                  - (U[p+4]*U[p+8] - U[p+5]*U[p+9]); // re(u2.v3-u3.v2)
        U[p+13] = - (U[p+2]*U[p+11] + U[p+3]*U[p+10])
                  + (U[p+4]*U[p+9] + U[p+5]*U[p+8]); // im(u2.v3-u3.v2)
        U[p+14] =   (U[p+4]*U[p+6] - U[p+5]*U[p+7])
                  - (U[p]*U[p+10] - U[p+1]*U[p+11]); // re(u3.v1-u1.v3)
        U[p+15] = - (U[p+4]*U[p+7] + U[p+5]*U[p+6])
                  + (U[p]*U[p+11] + U[p+1]*U[p+10]); // im(u3.v1-u1.v3)
        U[p+16] =   (U[p]*U[p+8] - U[p+1]*U[p+9])
                  - (U[p+2]*U[p+6] - U[p+3]*U[p+7]); // re(u1.v2-u2.v1)
        U[p+17] = - (U[p]*U[p+9] + U[p+1]*U[p+8])
                  + (U[p+2]*U[p+7] + U[p+3]*U[p+6]); // im(u1.v2-u2.v1)
    }
}

void gram_schmidt_per(double *link) {
    double temp[X];
    double dot_real, dot_imag;
    double norm_u2, norm_u;
    double norm_v2, norm_v;

    memcpy(temp,link,sizeof(double)*X);

    // get real and imaginary parts of dot product
    dot_real = temp[0]*temp[6] + temp[1]*temp[7]
             + temp[2]*temp[8] + temp[3]*temp[9]
             + temp[4]*temp[10] + temp[5]*temp[11];
    dot_imag = temp[0]*temp[7] - temp[1]*temp[6]
             + temp[2]*temp[9] - temp[3]*temp[8]
             + temp[4]*temp[11] - temp[5]*temp[10];

    // norm of u vector
    norm_u2 = temp[0]*temp[0] + temp[1]*temp[1]
            + temp[2]*temp[2] + temp[3]*temp[3]
            + temp[4]*temp[4] + temp[5]*temp[5];
    norm_u = sqrt(norm_u2);

    // construct orthogonal vector for second row
    temp[6] = temp[6] - (temp[0]*dot_real - temp[1]*dot_imag)/norm_u2; 
    temp[7] = temp[7] - (temp[0]*dot_imag + temp[1]*dot_real)/norm_u2; 
    temp[8] = temp[8] - (temp[2]*dot_real - temp[3]*dot_imag)/norm_u2; 
    temp[9] = temp[9] - (temp[2]*dot_imag + temp[3]*dot_real)/norm_u2; 
    temp[10] = temp[10] - (temp[4]*dot_real - temp[5]*dot_imag)/norm_u2; 
    temp[11] = temp[11] - (temp[4]*dot_imag + temp[5]*dot_real)/norm_u2; 

    // norm of v vector
    norm_v2 = temp[6]*temp[6] + temp[7]*temp[7]
            + temp[8]*temp[8] + temp[9]*temp[9]
            + temp[10]*temp[10] + temp[11]*temp[11];
    norm_v = sqrt(norm_v2);

    link[0] = temp[0]/norm_u; // u1.real
    link[1] = temp[1]/norm_u; // u1.imag
    link[2] = temp[2]/norm_u; // u2.real
    link[3] = temp[3]/norm_u; // u2.imag
    link[4] = temp[4]/norm_u; // u3.real
    link[5] = temp[5]/norm_u; // u3.imag
    link[6] = temp[6]/norm_v; // v1.real
    link[7] = temp[7]/norm_v; // v1.imag
    link[8] = temp[8]/norm_v;  // v2.real
    link[9] = temp[9]/norm_v;  // v2.imag
    link[10] = temp[10]/norm_v; // v3.real
    link[11] = temp[11]/norm_v; // v3.imag
    link[12] =   (link[2]*link[10] - link[3]*link[11])
              - (link[4]*link[8] - link[5]*link[9]); // re(u2.v3-u3.v2)
    link[13] = - (link[2]*link[11] + link[3]*link[10])
              + (link[4]*link[9] + link[5]*link[8]); // im(u2.v3-u3.v2)
    link[14] =   (link[4]*link[6] - link[5]*link[7])
              - (link[0]*link[10] - link[1]*link[11]); // re(u3.v1-u1.v3)
    link[15] = - (link[4]*link[7] + link[5]*link[6])
              + (link[0]*link[11] + link[1]*link[10]); // im(u3.v1-u1.v3)
    link[16] =   (link[0]*link[8] - link[1]*link[9])
              - (link[2]*link[6] - link[3]*link[7]); // re(u1.v2-u2.v1)
    link[17] = - (link[0]*link[9] + link[1]*link[8])
              + (link[2]*link[7] + link[3]*link[6]); // im(u1.v2-u2.v1)
}
