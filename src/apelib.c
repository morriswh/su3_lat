#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <inttypes.h>
#include <lapacke.h>
#include <complex.h>
// #include <gmp.h>
// #include <mpfr.h>
#include <math.h>
#include <time.h>
#include "sizes.h"
#include "matops.h"
#include "su3lib.h"
#include "apelib.h"

enum dim{mu_x, mu_y, mu_z, mu_t};

void ape_smear_yz(double *V, double *U, double alpha) {

    // make coefficients
    double one_minus_alpha = 1 - alpha;
    double alpha_div_six = alpha/6.0;

    int i, mu, nu;
    #pragma omp parallel for collapse(2) private(i, mu, nu)
    for (i = 0; i < LVOLUME; i++) {
        for (mu = 0; mu < D; mu++) {
            if ( mu == mu_y || mu == mu_z ) {
                double base[X], staple[X], staple_sum[X];
                bzero(staple_sum, sizeof(staple_sum));
                memcpy(base, &U[X*(D*i + mu)], sizeof(double)*X);
                for (nu = 1; nu < D-1; nu++) { // what am I supposed to include here?
                    if (nu == mu) {            // defer to what GLU does for now
                        continue;
                    }
                    else {
                        staples(staple, U, i, mu, nu);
                        for (int j = 0; j < X; j++) {
                            staple_sum[j] += staple[j];
                        }
                    }
                }
                cmmul(base, one_minus_alpha, base);
                cmmul(staple_sum, alpha_div_six, staple_sum);
                for (int j = 0; j < X; j++) {
                    base[j] += staple_sum[j];
                }
                memcpy(&V[X*(D*i + mu)], base, sizeof(double)*X);
                newton_reunit(&V[X*(D*i + mu)]);
            }
            else {
                memcpy(&V[X*(D*i + mu)], &U[X*(D*i + mu)], sizeof(double)*X);
            }
        }
    }
}

void gram_schmidt(double *U) {

    int i;
    #pragma omp parallel for private(i)
    for (i = 0; i < LLENGTH; i++) {

        double temp[X];
        double dot_real, dot_imag;
        double norm_u2, norm_u;
        double norm_v2, norm_v;

        memcpy(temp,&U[X*i],sizeof(double)*X);

        int p = X*i; // position in big array for config
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

void newton_reunit(double *m) {

    // Either rewrite this without using complex.h, or rewrite everything
    // to use complex.h. Might be better to just do that

    double complex mc[9];
    for (int i = 0; i < 9; i++) {
        mc[i] = m[2*i] + m[2*i+1]*I;
    }

    double detm[2];
    det(detm, m);

    double complex detmc = detm[0] + detm[1]*I;

    double complex one_div_det_cube_root = 1.0/cpowf(detmc, 1.0/3.0);

    for (int i = 0; i < 9; i++) {
        mc[i] = one_div_det_cube_root*mc[i];
    }

    for (int i = 0; i < 9; i++) {
        m[2*i] = creal(mc[i]);
        m[2*i+1] = cimag(mc[i]);
    }

    double delta = 1.0;

    double x[X], xai[X], xn[X];
    memcpy(x, mc, X*sizeof(double));

    while (delta > 10E-15) {

        delta = 0.0;
        double xa[X], xi[X];
        adjoint(xa, x);
        inverse(xai, xa);

        for (int i = 0; i < X; i++) {
            xn[i] = 0.5*(x[i] + xai[i]);
        }

        adjoint(xa,xn);
        matmul(xi,xa,xn);
        for (int i = 0; i < X; i++) {
            delta += fabs(xi[i]);
        }
        delta -= 3.0;
        // printf("%.15lf\n", delta);

        memcpy(x, xn, X*sizeof(double));
    }

    memcpy(m, x, X*sizeof(double));

}











