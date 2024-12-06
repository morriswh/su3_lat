#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <inttypes.h>
// #include <gmp.h>
// #include <mpfr.h>
#include <math.h>
#include <time.h>
#include "sizes.h"
#include "su3lib.h"
#include "apelib.h"

enum dim{mu_x, mu_y, mu_z, mu_t};

void ape_smear_yz(double *V, double *U, int steps, double alpha) {

    // make coefficients
    double one_minus_alpha_matrix[X], alpha_matrix_div_six[X];

    memset(one_minus_alpha_matrix, 0, sizeof(one_minus_alpha_matrix));
    one_minus_alpha_matrix[0] = 1 - alpha;
    one_minus_alpha_matrix[8] = 1 - alpha;
    one_minus_alpha_matrix[16] = 1 - alpha;

    memset(alpha_matrix, 0, sizeof(alpha_matrix));
    alpha_matrix_div_six[0] = alpha/6.0;
    alpha_matrix_div_six[8] = alpha/6.0;
    alpha_matrix_div_six[16] = alpha/6.0;

    double base[X], staple_sum[X];
    memset(staple_sum, 0, sizeof(staple_sum));

    // int mu = mu_y;
    int i, mu, nu, j;
    for (i = 0; i < LVOLUME; i++) {
        for (mu = mu_y; mu <= mu_z; mu++) {
            memcpy(base, &U[X*(D*i + mu)]);
            for (nu = 0; nu < D; nu++) {
                if (nu == mu) {
                    continue;
                }
                staples(staple, U, i, mu, nu);
                for (j = 0; j < X; j++) {
                    staple_sum[j] += staple[j];
                }
            }
            matmul(base, one_minus_alpha_matrix, base);
            matmul(staple_sum, alpha_matrix_div_six, staple_sum);
            for (j = 0; j < X; j++) {
                base[j] += staple_sum[j];
            }
            memcpy(&V[X*(D*i + mu)], base);
        }
    }
    gram_schmidt(V);
}

void gram_schmidt(double *U) {

    double temp[X];
    double dot_real, dot_imag;
    double norm_u2, norm_u;
    double norm_v2, norm_v;
    int i, j, p;
    for (i = 0; i < LLENGTH; i++) {
        memcpy(temp,&U[X*i],sizeof(double)*X)
        // for (j = 0; j < X; j++) {
        //     temp[j] = U[];
        // }
        p = X*i; // position in big array for config
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

int ipow(int base, int exp)
{
    // borrowed this from a stackexchange post
    int result = 1;
    for (;;)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        if (!exp)
            break;
        base *= base;
    }

    return result;
}

int mod(int value, int period) {
    // avoid taking unnecessary operations
    if (value >= 0 && value < period) {
        return value;
    }
    else {
        return (value%period + period)%period;
    }
    // return (value%period + period)%period;
}

