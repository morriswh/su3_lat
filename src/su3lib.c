#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <inttypes.h>
// #include <gmp.h>
// #include <mpfr.h>
#include <math.h>
#include <time.h>
#include "sizes.h"
#include "su3lib.h"

#define LENGTH(X) (sizeof X / sizeof X[0])

enum dim{mu_x, mu_y, mu_z, mu_t};

void unit_config(double *U) {
    // make a unit gauge configuration
    for (int64_t i = 0; i < LLENGTH; ++i) {
        int64_t p = X*i; // position in big array for config
        // printf("%ld\n", p);

        U[p]    = 1.0;
        U[p+1]  = 0.0;
        U[p+2]  = 0.0;
        U[p+3]  = 0.0;
        U[p+4]  = 0.0;
        U[p+5]  = 0.0;

        U[p+6]  = 0.0;
        U[p+7]  = 0.0;
        U[p+8]  = 1.0;
        U[p+9]  = 0.0;
        U[p+10] = 0.0;
        U[p+11] = 0.0;

        U[p+12] = 0.0;
        U[p+13] = 0.0;
        U[p+14] = 0.0;
        U[p+15] = 0.0;
        U[p+16] = 1.0;
        U[p+17] = 0.0;
    }
}

void rnd_config(double *U) {
    // make a random gauge configuration
    // double U[CLENGTH];
    double temp[S_X];
    double dot_real, dot_imag;
    double norm_u2, norm_u;
    double norm_v2, norm_v;
    int64_t i, j, p;
    for (i = 0; i < LLENGTH; i++) {
        for (j = 0; j < S_X; j++) {
            temp[j] = 2*(double)rand()/(double)RAND_MAX - 1;
        }
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

void rnd_gauge_tr(double *Omega) {
    // make a random gauge transformation matrix
    // double Omega[OLENGTH];
    double temp[S_X];
    double dot_real, dot_imag;
    double norm_u2, norm_u;
    double norm_v2, norm_v;
    int64_t i, j, p;
    for (i = 0; i < LVOLUME; i++) {
        for (j = 0; j < S_X; j++) {
            temp[j] = 2*(double)rand()/(double)RAND_MAX - 1;
        }
        p = X*i; // position in big array
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

        Omega[p] = temp[0]/norm_u; // u1.real
        Omega[p+1] = temp[1]/norm_u; // u1.imag
        Omega[p+2] = temp[2]/norm_u; // u2.real
        Omega[p+3] = temp[3]/norm_u; // u2.imag
        Omega[p+4] = temp[4]/norm_u; // u3.real
        Omega[p+5] = temp[5]/norm_u; // u3.imag
        Omega[p+6] = temp[6]/norm_v; // v1.real
        Omega[p+7] = temp[7]/norm_v; // v1.imag
        Omega[p+8] = temp[8]/norm_v;  // v2.real
        Omega[p+9] = temp[9]/norm_v;  // v2.imag
        Omega[p+10] = temp[10]/norm_v; // v3.real
        Omega[p+11] = temp[11]/norm_v; // v3.imag
        Omega[p+12] =   (Omega[p+2]*Omega[p+10] - Omega[p+3]*Omega[p+11])
                      - (Omega[p+4]*Omega[p+8] - Omega[p+5]*Omega[p+9]); // re(u2.v3-u3.v2)
        Omega[p+13] = - (Omega[p+2]*Omega[p+11] + Omega[p+3]*Omega[p+10])
                      + (Omega[p+4]*Omega[p+9] + Omega[p+5]*Omega[p+8]); // im(u2.v3-u3.v2)
        Omega[p+14] =   (Omega[p+4]*Omega[p+6] - Omega[p+5]*Omega[p+7])
                      - (Omega[p]*Omega[p+10] - Omega[p+1]*Omega[p+11]); // re(u3.v1-u1.v3)
        Omega[p+15] = - (Omega[p+4]*Omega[p+7] + Omega[p+5]*Omega[p+6])
                      + (Omega[p]*Omega[p+11] + Omega[p+1]*Omega[p+10]); // im(u3.v1-u1.v3)
        Omega[p+16] =   (Omega[p]*Omega[p+8] - Omega[p+1]*Omega[p+9])
                      - (Omega[p+2]*Omega[p+6] - Omega[p+3]*Omega[p+7]); // re(u1.v2-u2.v1)
        Omega[p+17] = - (Omega[p]*Omega[p+9] + Omega[p+1]*Omega[p+8])
                      + (Omega[p+2]*Omega[p+7] + Omega[p+3]*Omega[p+6]); // im(u1.v2-u2.v1)
    }
}

void get_index(int64_t i) {
    // int64_t x = i%X;
    // int64_t mu = (i/X)%D;
    // int64_t l = ((i/X)/D)%L;
    // int64_t t = (((i/X)/D)/L)%T;
    // printf("%ld , %ld , %ld , %ld \n", t, l, mu, x);
    // int64_t mu = i%D;
    int64_t x = i%L;
    int64_t y = (i/L)%L;
    int64_t z = (i/LSQUARE)%L;
    int64_t t = (i/LCUBE)%T;
    // printf("%ld , %ld , %ld , %ld , %ld \n", x, y, z, t, mu);
    printf("%ld , %ld , %ld , %ld \n", x, y, z, t);
}

void get_link(double *link, double *U, int64_t pos) {
    pos = mod(pos, LLENGTH);
    pos *= X;
    memcpy(link,&U[pos],sizeof(double)*X);
    // for (int64_t i = 0; i < X; i++) {
    //     link[i] = U[pos+i];
    // }
}

void get_site(double *site, double *S, int64_t spos) {
    spos = mod(spos, LVOLUME);
    spos *= X;
    memcpy(site,&S[spos],sizeof(double)*X);
    // for (int64_t i = 0; i < X; i++) {
    //     site[i] = S[spos+i];
    // }
}

int64_t ipow(int64_t base, int64_t exp)
{
    // borrowed this from a stackexchange post
    int64_t result = 1;
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

int64_t mod(int64_t value, int64_t period) {
    // avoid taking unnecessary operations
    if (value >= 0 && value < period) {
        return value;
    }
    else {
        return (value%period + period)%period;
    }
    // return (value%period + period)%period;
}

void matmul(double *z, double *x, double *y) {
    double u[X];
    double v[X];
    memcpy(u,x,X*sizeof(double));
    memcpy(v,y,X*sizeof(double));
    // for (int64_t i = 0; i < X; i++) {
    //     u[i] = x[i];
    //     v[i] = y[i];
    // }

    z[0]  = u[0]*v[0]   - u[1]*v[1]   + u[2]*v[6]          // m[0] + i m[1]     m[2] + i m[3]     m[4] + i m[5]
          - u[3]*v[7]   + u[4]*v[12]  - u[5]*v[13];        // 
    z[1]  = u[0]*v[1]   + u[1]*v[0]   + u[2]*v[7]          // 
          + u[3]*v[6]   + u[4]*v[13]  + u[5]*v[12];        // 
    z[2]  = u[0]*v[2]   - u[1]*v[3]   + u[2]*v[8]          // m[6] + i m[7]     m[8] + i m[9]     m[10] + i m[11]
          - u[3]*v[9]   + u[4]*v[14]  - u[5]*v[15];        // 
    z[3]  = u[0]*v[3]   + u[1]*v[2]   + u[2]*v[9]          // 
          + u[3]*v[8]   + u[4]*v[15]  + u[5]*v[14];        // 
    z[4]  = u[0]*v[4]   - u[1]*v[5]   + u[2]*v[10]         // m[12] + i m[13]   m[14] + i m[15]   m[16] + i m[17]
          - u[3]*v[11]  + u[4]*v[16]  - u[5]*v[17];
    z[5]  = u[0]*v[5]   + u[1]*v[4]   + u[2]*v[11]
          + u[3]*v[10]  + u[4]*v[17]  + u[5]*v[16];
    z[6]  = u[6]*v[0]   - u[7]*v[1]   + u[8]*v[6]
          - u[9]*v[7]   + u[10]*v[12] - u[11]*v[13];
    z[7]  = u[6]*v[1]   + u[7]*v[0]   + u[8]*v[7]
          + u[9]*v[6]   + u[10]*v[13] + u[11]*v[12];
    z[8]  = u[6]*v[2]   - u[7]*v[3]   + u[8]*v[8]
          - u[9]*v[9]   + u[10]*v[14] - u[11]*v[15];
    z[9]  = u[6]*v[3]   + u[7]*v[2]   + u[8]*v[9]
          + u[9]*v[8]   + u[10]*v[15] + u[11]*v[14];
    z[10] = u[6]*v[4]   - u[7]*v[5]   + u[8]*v[10]
          - u[9]*v[11]  + u[10]*v[16] - u[11]*v[17];
    z[11] = u[6]*v[5]   + u[7]*v[4]   + u[8]*v[11]
          + u[9]*v[10]  + u[10]*v[17] + u[11]*v[16];
    z[12] = u[12]*v[0]  - u[13]*v[1]  + u[14]*v[6]
          - u[15]*v[7]  + u[16]*v[12] - u[17]*v[13];
    z[13] = u[12]*v[1]  + u[13]*v[0]  + u[14]*v[7]
          + u[15]*v[6]  + u[16]*v[13] + u[17]*v[12];
    z[14] = u[12]*v[2]  - u[13]*v[3]  + u[14]*v[8]
          - u[15]*v[9]  + u[16]*v[14] - u[17]*v[15];
    z[15] = u[12]*v[3]  + u[13]*v[2]  + u[14]*v[9]
          + u[15]*v[8]  + u[16]*v[15] + u[17]*v[14];
    z[16] = u[12]*v[4]  - u[13]*v[5]  + u[14]*v[10]
          - u[15]*v[11] + u[16]*v[16] - u[17]*v[17];
    z[17] = u[12]*v[5]  + u[13]*v[4]  + u[14]*v[11]
          + u[15]*v[10] + u[16]*v[17] + u[17]*v[16];
}

void cmmul(double *z, double c, double *x) {
    for (int32_t i = 0; i < X; ++i) {
        z[i] = c*x[i];
    }
}

void immul(double *z, double c, double *x) {
    // multiply an imaginary number by a complex matrix
    for (int32_t i = 0; i < X/2; ++i) {
        int32_t j = 2*i;
        z[j]   = -c*x[j+1]; // real 
        z[j+1] = c*x[j]; // imag
    }
}

void matadd(double *z, double *x, double *y) {
    for (int32_t i = 0; i < X; ++i) {
        z[i] = x[i] + y[i];
    }
}

void adjoint(double *xc, double *x) {
    double y[X];
    for (int64_t i = 0; i < X; i++) {
        y[i] = x[i];
    }
    xc[0] = y[0];       // m[0] + i m[1]         m[2] + i m[3]         m[4] + i m[5]
    xc[1] = -y[1];      // 
    xc[2] = y[6];       // 
    xc[3] = -y[7];      // 
    xc[4] = y[12];      // m[6] + i m[7]         m[8] + i m[9]         m[10] + i m[11]
    xc[5] = -y[13];     // 
    xc[6] = y[2];       // 
    xc[7] = -y[3];      // 
    xc[8] = y[8];       // m[12] + i m[13]       m[14] + i m[15]       m[16] + i m[17]
    xc[9] = -y[9];
    xc[10] = y[14];
    xc[11] = -y[15];
    xc[12] = y[4];
    xc[13] = -y[5];
    xc[14] = y[10];
    xc[15] = -y[11];
    xc[16] = y[16];
    xc[17] = -y[17];
}

void print_link(double *link) {
    // printf("[(%lf)+1j*(%lf),(%lf)+1j*(%lf),(%lf)+1j*(%lf)]\n",
    //        link[0], link[1], link[2], link[3], link[4], link[5]);
    // printf("[(%lf)+1j*(%lf),(%lf)+1j*(%lf),(%lf)+1j*(%lf)]\n",
    //        link[6], link[7], link[8], link[9], link[10], link[11]);
    // printf("[(%lf)+1j*(%lf),(%lf)+1j*(%lf),(%lf)+1j*(%lf)]\n",
    //        link[12], link[13], link[14], link[15], link[16], link[17]);
    // printf("np.array([[(%.15.15lf)+1j*(%.15lf),(%.15lf)+1j*(%.15lf),(%.15lf)+1j*(%.15lf)], [(%.15lf)+1j*(%.15lf),(%.15lf)+1j*(%.15lf),(%.15lf)+1j*(%.15lf)],[(%.15lf)+1j*(%.15lf),(%.15lf)+1j*(%.15lf),(%.15lf)+1j*(%.15lf)]])\n",
    printf("[[(%.15lf)+1j*(%.15lf),(%.15lf)+1j*(%.15lf),(%.15lf)+1j*(%.15lf)], [(%.15lf)+1j*(%.15lf),(%.15lf)+1j*(%.15lf),(%.15lf)+1j*(%.15lf)],[(%.15lf)+1j*(%.15lf),(%.15lf)+1j*(%.15lf),(%.15lf)+1j*(%.15lf)]]\n",
           link[0], link[1], link[2], link[3], link[4], link[5],
           link[6], link[7], link[8], link[9], link[10], link[11],
           link[12], link[13], link[14], link[15], link[16], link[17]);
}

void trace(double *result, double *matrix) {
    result[0] = matrix[0] + matrix[8] + matrix[16];
    result[1] = matrix[1] + matrix[9] + matrix[17];
}

void det(double *result, double *m) {
    double min1_real, min2_real, min3_real;
    double min1_imag, min2_imag, min3_imag;

    min1_real = -m[10]*m[14] + m[11]*m[15] + m[8]*m[16] - m[9]*m[17];
    min1_imag = -m[11]*m[14] - m[10]*m[15] + m[9]*m[16] + m[8]*m[17];

    min2_real = -m[10]*m[12] + m[11]*m[13] + m[6]*m[16] - m[7]*m[17];
    min2_imag = -m[11]*m[12] - m[10]*m[13] + m[7]*m[16] + m[6]*m[17];

    min3_real = -m[8]*m[12] + m[9]*m[13] + m[6]*m[14] - m[7]*m[15];
    min3_imag = -m[9]*m[12] - m[8]*m[13] + m[7]*m[14] + m[6]*m[15];

    result[0] = min1_real*m[0] - min1_imag*m[1] - min2_real*m[2] + min2_imag*m[3] + min3_real*m[4] - min3_imag*m[5];
    result[1] = min1_imag*m[0] + min1_real*m[1] - min2_imag*m[2] - min2_real*m[3] + min3_imag*m[4] + min3_real*m[5];
}

void plaquette(double *plaq, double *U, int64_t spos, int64_t mu, int64_t nu) {
    // spos is in LVOLUME now
    int64_t step_mu, step_nu, stride_mu, stride_nu;

    if (mu == 0) {
        step_mu = 1;
        stride_mu = L;
    }
    else if (mu == 1) {
        step_mu = L;
        stride_mu = LSQUARE;
    }
    else if (mu == 2) {
        step_mu = LSQUARE;
        stride_mu = LCUBE;
    }
    else {
        step_mu = LCUBE;
        stride_mu = LVOLUME;
    }

    if (nu == 0) {
        step_nu = 1;
        stride_nu = L;
    }
    else if (nu == 1) {
        step_nu = L;
        stride_nu = LSQUARE;
    }
    else if (nu == 2) {
        step_nu = LSQUARE;
        stride_nu = LCUBE;
    }
    else {
        step_nu = LCUBE;
        stride_nu = LVOLUME;
    }

    int64_t pos = D*spos + mu;
    double temp[X];
    get_link(plaq, U, mod(pos, LLENGTH));
    // matmul(plaq,plaq,temp);

    pos = D*(spos - mod(spos, stride_mu) + mod(spos + step_mu, stride_mu)) + nu;
    get_link(temp, U, mod(pos, LLENGTH));
    matmul(plaq,plaq,temp);

    // pos = D*(spos+step_nu) + mu;
    pos = D*(spos - mod(spos, stride_nu) + mod(spos + step_nu, stride_nu)) + mu;
    get_link(temp, U, mod(pos, LLENGTH));
    adjoint(temp, temp);
    matmul(plaq,plaq,temp);

    pos = D*spos + nu;
    get_link(temp, U, mod(pos, LLENGTH));
    adjoint(temp,temp);
    matmul(plaq,plaq,temp);
}

void all_plaquettes(double *U, double *a_plaq, double *s_plaq, double *t_plaq) {
    double splaq = 0.0, tplaq = 0.0;
    int64_t i;
    #pragma omp parallel for private(i) reduction(+:splaq) reduction(+:tplaq)
    for (i = 0; i < LVOLUME; i++) {
        double tr[2], plaq[X], help = 0.0;
        int64_t mu, nu;
        for (mu = 0; mu < D-1; mu++) {
            for (nu = mu+1; nu < D-1; nu++) {
                plaquette(plaq, U, i, mu, nu);
                trace(tr, plaq);
                help += tr[0];
            }
        }
        splaq = splaq + help;
        help = 0.0;
        for (mu = 0; mu < nu; mu++) {
            plaquette(plaq, U, i, mu, nu);
            trace(tr, plaq);
            help += tr[0];
        }
        tplaq = tplaq + help;
    }
    *s_plaq = 2.0 * splaq / (double)( ( D-1 ) * ( D-2 ) * LVOLUME * NC );
    *t_plaq = 2.0 * tplaq / (double)( ( D-1 ) * ( D-2 ) * LVOLUME * NC );
    *a_plaq = 0.5*( *s_plaq + *t_plaq );
}

void cloverleaf(double *clover, double *U, int64_t spos, int64_t mu, int64_t nu) {
    memset(clover, 0, X*sizeof(double));

    double temp[X], leaf[X];
    // double clover_l1[X], clover_l2[X], clover_l3[X], clover_l4[X];
    int64_t pos;
    int64_t step_mu, step_nu, stride_mu, stride_nu;

    if (mu == 0) {
        step_mu = 1;
        stride_mu = L;
    }
    else if (mu == 1) {
        step_mu = L;
        stride_mu = LSQUARE;
    }
    else if (mu == 2) {
        step_mu = LSQUARE;
        stride_mu = LCUBE;
    }
    else {
        step_mu = LCUBE;
        stride_mu = LVOLUME;
    }

    if (nu == 0) {
        step_nu = 1;
        stride_nu = L;
    }
    else if (nu == 1) {
        step_nu = L;
        stride_nu = LSQUARE;
    }
    else if (nu == 2) {
        step_nu = LSQUARE;
        stride_nu = LCUBE;
    }
    else {
        step_nu = LCUBE;
        stride_nu = LVOLUME;
    }

            // Top right leaf

            pos = D*spos + mu;
            get_link(leaf, U, mod(pos, LLENGTH));

            pos = D*(spos - mod(spos, stride_mu) + mod(spos + step_mu, stride_mu)) + nu;
            get_link(temp, U, mod(pos, LLENGTH));
            matmul(leaf, leaf, temp);

            pos = D*(spos - mod(spos, stride_nu) + mod(spos + step_nu, stride_nu)) + mu;
            get_link(temp, U, mod(pos, LLENGTH));
            adjoint(temp, temp);
            matmul(leaf, leaf, temp);

            pos = D*spos + nu;
            get_link(temp, U, mod(pos, LLENGTH));
            adjoint(temp, temp);
            matmul(leaf, leaf, temp);

            for (int64_t i = 0; i < X; i++) {
                clover[i] += leaf[i];
            }

            // Top left leaf

            pos = D*spos + nu;
            get_link(leaf, U, mod(pos, LLENGTH));

            pos = D*(spos - mod(spos, stride_mu) - mod(spos, stride_nu)
                + mod(spos - step_mu, stride_mu) + mod(spos + step_nu, stride_nu)) + mu;
            get_link(temp, U, mod(pos, LLENGTH));
            adjoint(temp, temp);
            matmul(leaf, leaf, temp);

            pos = D*(spos - mod(spos, stride_mu) + mod(spos - step_mu, stride_mu)) + nu;
            get_link(temp, U, mod(pos, LLENGTH));
            adjoint(temp, temp);
            matmul(leaf, leaf, temp);

            pos = D*(spos - mod(spos, stride_mu) + mod(spos - step_mu, stride_mu)) + mu;
            get_link(temp, U, mod(pos, LLENGTH));
            matmul(leaf, leaf, temp);

            for (int64_t i = 0; i < X; i++) {
                clover[i] += leaf[i];
            }

            // Bottom left leaf

            pos = D*(spos - mod(spos, stride_mu) + mod(spos - step_mu, stride_mu)) + mu;
            get_link(leaf, U, mod(pos, LLENGTH));
            adjoint(leaf, leaf);

            pos = D*(spos - mod(spos, stride_mu) - mod(spos, stride_nu)
                + mod(spos - step_mu, stride_mu) + mod(spos - step_nu, stride_nu)) + nu;
            get_link(temp, U, mod(pos, LLENGTH));
            adjoint(temp, temp);
            matmul(leaf, leaf, temp);

            pos = D*(spos - mod(spos, stride_mu) - mod(spos, stride_nu)
                + mod(spos - step_mu, stride_mu) + mod(spos - step_nu, stride_nu)) + mu;
            get_link(temp, U, mod(pos, LLENGTH));
            matmul(leaf, leaf, temp);

            pos = D*(spos - mod(spos, stride_nu) + mod(spos - step_nu, stride_nu)) + nu;
            get_link(temp, U, mod(pos, LLENGTH));
            matmul(leaf, leaf, temp);

            for (int64_t i = 0; i < X; i++) {
                clover[i] += leaf[i];
            }

            // Bottom right leaf

            pos = D*(spos - mod(spos, stride_nu) + mod(spos - step_nu, stride_nu)) + nu;
            get_link(leaf, U, mod(pos, LLENGTH));
            adjoint(leaf, leaf);

            pos = D*(spos - mod(spos, stride_nu) + mod(spos - step_nu, stride_nu)) + mu;
            get_link(temp, U, mod(pos, LLENGTH));
            matmul(leaf, leaf, temp);

            pos = D*(spos - mod(spos, stride_mu) - mod(spos, stride_nu)
                + mod(spos + step_mu, stride_mu) + mod(spos - step_nu, stride_nu)) + nu;
            get_link(temp, U, mod(pos, LLENGTH));
            matmul(leaf, leaf, temp);

            pos = D*spos + mu;
            get_link(temp, U, mod(pos, LLENGTH));
            adjoint(temp, temp);
            matmul(leaf, leaf, temp);

            for (int64_t i = 0; i < X; i++) {
                clover[i] += leaf[i];
            }

    adjoint(temp, clover);
    for (int64_t i = 0; i < X; i++) {
        clover[i] -= temp[i];
        clover[i] *= 1.0/8.0;
    }
    // print_link(clover);
}

void top_charge_density(double *topo_charge_density, double *U, int64_t spos) {
    double tcharge[2] = {0};
    // memcpy(tcharge, topo_charge, 2*sizeof(double));
    // memset(tcharge, 0, 2*sizeof(double));

    int64_t i;
    double C1[X], C2[X], C3[X], C4[X], C5[X], C6[X];
    double O1[X], O2[X], O3[X], O4[X], O5[X], O6[X];
    double PT[X];
    double AT[2];

    cloverleaf(C1, U, spos, 0, 1);
    cloverleaf(C2, U, spos, 0, 2);
    cloverleaf(C3, U, spos, 0, 3);
    cloverleaf(C4, U, spos, 2, 1);
    cloverleaf(C5, U, spos, 1, 3);
    cloverleaf(C6, U, spos, 3, 2);

    matmul(O1, C1, C6);
    matmul(O2, C2, C5);
    matmul(O3, C3, C4);
    matmul(O4, C4, C3);
    matmul(O5, C5, C2);
    matmul(O6, C6, C1);

    for (i = 0; i < X; i++) {
        PT[i] = O1[i] + O2[i] + O3[i] + O4[i] + O5[i] + O6[i];
    }

    trace(AT, PT);
    tcharge[0] += AT[0];
    tcharge[1] += AT[1];

    topo_charge_density[0] = tcharge[0]/(8.0*pow(PI, 2));
    topo_charge_density[1] = tcharge[1]/(8.0*pow(PI, 2));
}

void top_charge(double *topo_charge, double *U) {
    double tcharge[2] = {0};

    #pragma omp parallel for reduction(+: tcharge)
    for (int64_t spos = 0; spos < LVOLUME; spos++) {

        double tcharged[2] = {0};

        top_charge_density(tcharged, U, spos);

        tcharge[0] += tcharged[0];
        tcharge[1] += tcharged[1];
    }

    topo_charge[0] = tcharge[0];
    topo_charge[1] = tcharge[1];
}

void directional_sum(double *lattice, int64_t *directions, 
                     double *out_array, int64_t out_array_size) {

    // First, we determine the size of the output array, and the size of
    // elements to be summed over. We also create two shape arrays to help
    // with the logic for the directional sums.

    int64_t sum_shape[] = {1, 1, 1, 1}; 
    int64_t out_shape[] = {1, 1, 1, 1}; 
    int64_t out_size = 1;
    int64_t sum_size = 1;

    for (int64_t i = 0; i < D-1; i ++) {
        if (directions[i] == 0) {
            out_size *= L;
            out_shape[i] = L;
        } else {
            sum_size *= L;
            sum_shape[i] = L;
        }
    }
    if (directions[3] == 0) {
        out_size *= T;
        out_shape[3] = T;
    } else {
        sum_size *= T;
        sum_shape[3] = T;
    }

    if (out_size != out_array_size) {
        printf("output array size mismatch!\n");
        exit(1);
    }

    for (int64_t i = 0; i < out_size; i++) {
        for (int64_t j = 0; j < sum_size; j++) {
            int64_t vali = i;
            int64_t valj = j;
            int64_t pos = 0;
            int64_t coef = 1;
            for (int64_t k = 0; k < D; k++) {
                pos += coef*((1-directions[k])*mod(vali, out_shape[k]) + 
                             directions[k]*mod(valj, sum_shape[k]));
                vali /= out_shape[k];
                valj /= sum_shape[k];
                coef *= L;
            }
            out_array[i] += lattice[pos];
        }
    }
}

void staples(double *staple, double *U, int64_t spos, int64_t mu, int64_t nu) {
    // spos is in LVOLUME now
    int64_t step_mu, step_nu, stride_mu, stride_nu;

    if (mu == 0) {
        step_mu = 1;
        stride_mu = L;
    }
    else if (mu == 1) {
        step_mu = L;
        stride_mu = LSQUARE;
    }
    else if (mu == 2) {
        step_mu = LSQUARE;
        stride_mu = LCUBE;
    }
    else {
        step_mu = LCUBE;
        stride_mu = LVOLUME;
    }

    if (nu == 0) {
        step_nu = 1;
        stride_nu = L;
    }
    else if (nu == 1) {
        step_nu = L;
        stride_nu = LSQUARE;
    }
    else if (nu == 2) {
        step_nu = LSQUARE;
        stride_nu = LCUBE;
    }
    else {
        step_nu = LCUBE;
        stride_nu = LVOLUME;
    }

    double temp[X], staple_h[X];

    // first staple
    int64_t pos = D*spos + nu;
    get_link(staple, U, mod(pos, LLENGTH));

    pos = D*(spos - mod(spos, stride_nu) + mod(spos + step_nu, stride_nu)) + mu;
    get_link(temp, U, mod(pos, LLENGTH));
    matmul(staple, staple, temp);

    pos = D*(spos - mod(spos, stride_mu) + mod(spos + step_mu, stride_mu)) + nu;
    get_link(temp, U, mod(pos, LLENGTH));
    adjoint(temp, temp);
    matmul(staple, staple, temp);

    // second staple
    pos = D*(spos - mod(spos, stride_nu) + mod(spos - step_nu, stride_nu)) + nu;
    get_link(staple_h, U, mod(pos, LLENGTH));
    adjoint(staple_h, staple_h);

    pos = D*(spos - mod(spos, stride_nu) + mod(spos - step_nu, stride_nu)) + mu;
    get_link(temp, U, mod(pos, LLENGTH));
    matmul(staple_h, staple_h, temp);
    
    pos = D*(spos - mod(spos, stride_nu) + mod(spos - step_nu, stride_nu) - mod(spos, stride_mu) + mod(spos + step_mu, stride_mu)) + nu;
    get_link(temp, U, mod(pos, LLENGTH));
    matmul(staple_h, staple_h, temp);

    for (int64_t i = 0; i < X; i++) {
        staple[i] += staple_h[i];
    }
}

void link_trace(double *ltr, double *U) {
    double res = 0;
    for (int64_t i = 0; i < LLENGTH; i++) {
        double tr[2];
        trace(tr, &U[X*i]);
        res += tr[0];
    }
    *ltr = res / ( LLENGTH * NC );
}

void transform(double *Up, double *U, double *Om) {
    int64_t pos, opos, fopos;
    int64_t offs, period, h;
    double link[X];
    double omegaL[X];
    double omegaR[X];
    for (int64_t i = 0; i < LLENGTH; i++) {
        pos = X*i;
        opos = X*(i/D);
        // offs = X*ipow(L,i%D); 
        h = i%D;
        
        if (h == 0) {
            offs = X;
            period = X*L;
        }
        else if (h == 1) {
            offs = X*L;
            period = X*LSQUARE;
        }
        else if (h == 2) {
            offs = X*LSQUARE;
            period = X*LCUBE;
        }
        else {
            offs = X*LCUBE;
            period = X*LVOLUME;
        }

        fopos = opos - mod(opos, period) + mod(opos + offs, period);
        for (int64_t j = 0; j < X; j++) {
            link[j] = U[pos+j];
            omegaL[j] = Om[opos+j];
            omegaR[j] = Om[mod(fopos+j, OLENGTH)];
        }
        adjoint(omegaR,omegaR);

        matmul(link,omegaL,link);
        matmul(link,link,omegaR);

        for (int64_t k = 0; k < X; k++) {
            Up[pos+k] = link[k];
        }
    }
}

void gauge_action(double *action, double *U) {
    double plaq[X];
    double sum = 0.0;
    double tr[2];
    int64_t i = 0, j = 0, k = 0;
    // I don't understand 'reduction' yet, but apparently it is needed
    #pragma omp parallel for reduction(+:sum) private(tr, plaq, i, j, k)
    for (i = 0; i < LVOLUME; i++) {
        double help = 0.0;
        for (j = 0; j < D; j++) {
            for (k = j+1; k < D; k++) {
                plaquette(plaq, U, i, j, k);
                trace(tr, plaq);
                help += tr[0];
            }
        }
        sum += help;
    }
    *action = 2*beta*(3*T*L - sum);
}

