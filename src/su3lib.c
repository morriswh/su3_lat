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

void rnd_config(double *U) {

    // make a random gauge configuration
    // double U[CLENGTH];
    double temp[S_X];
    double dot_real, dot_imag;
    double norm_u2, norm_u;
    double norm_v2, norm_v;
    int i, j, p;
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
    int i, j, p;
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

void get_index(int i) {
    // int x = i%X;
    // int mu = (i/X)%D;
    // int l = ((i/X)/D)%L;
    // int t = (((i/X)/D)/L)%T;
    // printf("%d, %d, %d, %d\n", t, l, mu, x);
    // int mu = i%D;
    int x = i%L;
    int y = (i/L)%L;
    int z = (i/LSQUARE)%L;
    int t = (i/LCUBE)%T;
    // printf("%d, %d, %d, %d, %d\n", x, y, z, t, mu);
    printf("%d, %d, %d, %d\n", x, y, z, t);
}

void get_link(double *link, double *U, int pos) {
    pos = mod(pos, LLENGTH);
    pos *= X;
    memcpy(link,&U[pos],sizeof(double)*X);
    // for (int i = 0; i < X; i++) {
    //     link[i] = U[pos+i];
    // }
}

void get_site(double *site, double *S, int spos) {
    spos = mod(spos, LVOLUME);
    spos *= X;
    memcpy(site,&S[spos],sizeof(double)*X);
    // for (int i = 0; i < X; i++) {
    //     site[i] = S[spos+i];
    // }
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

void matmul(double *z, double *x, double *y) {
    double u[X];
    double v[X];
    memcpy(u,x,X*sizeof(double));
    memcpy(v,y,X*sizeof(double));
    // for (int i = 0; i < X; i++) {
    //     u[i] = x[i];
    //     v[i] = y[i];
    // }

    z[0]  = u[0]*v[0]   - u[1]*v[1]   + u[2]*v[6]          // m[0] + i m[1]         m[2] + i m[3]         m[4] + i m[5]
          - u[3]*v[7]   + u[4]*v[12]  - u[5]*v[13];        // 
    z[1]  = u[0]*v[1]   + u[1]*v[0]   + u[2]*v[7]          // 
          + u[3]*v[6]   + u[4]*v[13]  + u[5]*v[12];        // 
    z[2]  = u[0]*v[2]   - u[1]*v[3]   + u[2]*v[8]          // m[6] + i m[7]         m[8] + i m[9]         m[10] + i m[11]
          - u[3]*v[9]   + u[4]*v[14]  - u[5]*v[15];        // 
    z[3]  = u[0]*v[3]   + u[1]*v[2]   + u[2]*v[9]          // 
          + u[3]*v[8]   + u[4]*v[15]  + u[5]*v[14];        // 
    z[4]  = u[0]*v[4]   - u[1]*v[5]   + u[2]*v[10]         // m[12] + i m[13]       m[14] + i m[15]       m[16] + i m[17]
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
    z[0]  = c*x[0];
    z[1]  = c*x[1];
    z[2]  = c*x[2];
    z[3]  = c*x[3];
    z[4]  = c*x[4];
    z[5]  = c*x[5];
    z[6]  = c*x[6];
    z[7]  = c*x[7];
    z[8]  = c*x[8];
    z[9]  = c*x[9];
    z[10] = c*x[10];
    z[11] = c*x[11];
    z[12] = c*x[12];
    z[13] = c*x[13];
    z[14] = c*x[14];
    z[15] = c*x[15];
    z[16] = c*x[16];
    z[17] = c*x[17];
}

void adjoint(double *xc, double *x) {
    double y[X];
    for (int i = 0; i < X; i++) {
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

void plaquette(double *plaq, double *U, int spos, int mu, int nu) {
    // spos is in LVOLUME now
    int step_mu, step_nu, stride_mu, stride_nu;

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

    int pos = D*spos + mu;
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
    int i;
    #pragma omp parallel for private(i) reduction(+:splaq) reduction(+:tplaq)
    for (i = 0; i < LVOLUME; i++) {
        double tr[2], plaq[X], help = 0.0;
        int mu, nu;
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

void cloverleaf(double *clover, double *U, int spos, int mu, int nu) {

    memset(clover, 0, X*sizeof(double));

    double temp[X], leaf[X];
    // double clover_l1[X], clover_l2[X], clover_l3[X], clover_l4[X];
    int pos;
    int step_mu, step_nu, stride_mu, stride_nu;

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

            for (int i = 0; i < X; i++) {
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

            for (int i = 0; i < X; i++) {
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

            for (int i = 0; i < X; i++) {
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

            for (int i = 0; i < X; i++) {
                clover[i] += leaf[i];
            }

    adjoint(temp, clover);
    for (int i = 0; i < X; i++) {
        clover[i] -= temp[i];
        clover[i] *= 1.0/8.0;
    }
    // print_link(clover);
}

void top_charge_density(double *topo_charge_density, double *U, int spos) {

    double tcharge[2] = {0};
    // memcpy(tcharge, topo_charge, 2*sizeof(double));
    // memset(tcharge, 0, 2*sizeof(double));

    int i;
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
    for (int spos = 0; spos < LVOLUME; spos++) {

        double tcharged[2] = {0};

        top_charge_density(tcharged, U, spos);

        tcharge[0] += tcharged[0];
        tcharge[1] += tcharged[1];
    }

    topo_charge[0] = tcharge[0];
    topo_charge[1] = tcharge[1];
}

void staples(double *staple, double *U, int spos, int mu, int nu) {
    // spos is in LVOLUME now
    int step_mu, step_nu, stride_mu, stride_nu;

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
    int pos = D*spos + nu;
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

    for (int i = 0; i < X; i++) {
        staple[i] += staple_h[i];
    }

}

void link_trace(double *ltr, double *U) {
    double res = 0;
    for (int i = 0; i < LLENGTH; i++) {
        double tr[2];
        trace(tr, &U[X*i]);
        res += tr[0];
    }
    *ltr = res / ( LLENGTH * NC );
}

void transform(double *Up, double *U, double *Om) {
    int pos, opos, fopos;
    int offs, period, h;
    double link[X];
    double omegaL[X];
    double omegaR[X];
    for (int i = 0; i < LLENGTH; i++) {
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
        for (int j = 0; j < X; j++) {
            link[j] = U[pos+j];
            omegaL[j] = Om[opos+j];
            omegaR[j] = Om[mod(fopos+j, OLENGTH)];
        }
        adjoint(omegaR,omegaR);

        matmul(link,omegaL,link);
        matmul(link,link,omegaR);

        for (int k = 0; k < X; k++) {
            Up[pos+k] = link[k];
        }
    }
}

void gauge_action(double *action, double *U) {
    double plaq[X];
    double sum = 0.0;
    double tr[2];
    int i = 0, j = 0, k = 0;
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

void read_nersc(double *U, char *config, char *header, int *header_length) {
    // need to figure out the issue with fread and doubles for the NERSC config
    // for now the union method works I guess
    union {
        long int retriever;
        double element;
    } package; 

    FILE *config_file = fopen(config, "rb");
    int safety = 2000;
    char check_end[] = "END_HEADER\n";
    *header_length = 0;
    for (int i = 0; i < safety; i++) {
        header[i] = fgetc(config_file);
        *header_length += 1;
        if ( i > 10) {
            if (strncmp(check_end, &header[i-10], sizeof(check_end)) == 0) {
                break;
            }
        }
    }

    // read in configuration data
    for (int i = 0; i != CLENGTH; i++) {
        fread(&package.retriever, sizeof(package.retriever), 1, config_file);
        package.retriever = be64toh(package.retriever);
        U[i] = package.element;
    }
    fclose(config_file);
}

void write_nersc(double *U, char *out_file, char *header, int header_length) {
    union {
        long int sender;
        double element;
    } package; 

    FILE *config_file = fopen(out_file, "wb");
    fwrite(header, sizeof(char), header_length, config_file);
    // fwrite(U, sizeof(double), CLENGTH, config_file);

    for (int i = 0; i != CLENGTH; i++) {
        package.element = U[i];
        package.sender = htobe64(package.sender);
        fwrite(&package.sender, sizeof(package.sender), 1, config_file);
    }
    fclose(config_file);
}

void find_value(char *value, char *label, int label_length, char *header, int header_length) {
    for (int i = 0; i < header_length; i++) {
        if (strncmp(&header[i], label, label_length*sizeof(char)) == 0) {
            i = i + label_length + 3; // the 3 is for ' = '
            int j = 0;
            while (header[i] != '\n') {
                value[j] = header[i];
                i++;
                j++;
            }
            value[j] = '\0';
            break;
        }
    }
}

void make_header(char *header, int *header_length, double *U, char* label) {
    // First gather all the information
    // Ultimately this should take a structure as an argument with
    // all the relavant information

    uint32_t cks = checksum(U);
    char cksm[20];
    snprintf(cksm, 20, "%x", cks);

    double tr;
    link_trace(&tr, U);
    char ltr[30];
    snprintf(ltr, 30*sizeof(char), "%.14lf", tr);

    double tplaq, splaq, aplaq;
    all_plaquettes(U, &aplaq, &splaq, &tplaq);
    char plaq[30];
    snprintf(plaq, 30*sizeof(char), "%.14lf", aplaq);

    /*char label[30];*/
    /*snprintf(label, 30*sizeof(char), "T%dL%d_ape", T, L);*/

    // char *user = getlogin();
    char user[10];
    snprintf(user, 10*sizeof(char), "wayne");

    time_t t = time(NULL);
    struct tm tms = *localtime(&t);
    char today[80];
    snprintf(today, 80*sizeof(char), 
             "%d-%02d-%02d %02d:%02d:%02d", 
             tms.tm_year + 1900, tms.tm_mon + 1, 
             tms.tm_mday, tms.tm_hour, tms.tm_min, tms.tm_sec);

    char *inputs[] = {
        cksm,
        ltr,
        plaq,
        label,
        user,
        today
    };

    // Now make the header
    char *new_header = malloc(sizeof(char)*1000);
    // char line[100];
    int old_index = 0, new_index = 0;
    // char line_end[] = "END_HEADER\n";
    char *lines[] = {
        "CHECKSUM\0",
        "LINK_TRACE\0",
        "PLAQUETTE\0",
        "ENSEMBLE_LABEL\0",
        "CREATOR_MACHINE\0",
        "CREATION_DATE\0"
    };
    while( old_index < *header_length ){
        new_header[new_index] = header[old_index];
        for (int j = 0; j < 6; j++ ) {
            int prop_size = 0;
            while (lines[j][prop_size] != '\0') {
                prop_size++;
            };
            if (strncmp(&header[old_index], lines[j], prop_size*sizeof(char)) == 0) {
                memcpy(&new_header[new_index], lines[j], prop_size*sizeof(char));
                new_index += prop_size;
                memcpy(&new_header[new_index], " = ", 3*sizeof(char));
                new_index += 3;
                int i = 0;
                while (inputs[j][i] != '\0') {
                    new_header[new_index] = inputs[j][i];
                    i++;
                    new_index++;
                };
                new_header[new_index] = '\n';
                while (header[old_index] != '\n') {
                    old_index++;
                };
            };
        };
        new_index++;
        old_index++;
    };
    new_header[new_index] = '\0';
    // // bzero(header, sizeof(char)*(*header_length));
    *header_length = new_index;
    memcpy(header, new_header, sizeof(char)*(*header_length+1));
    free(new_header);
}

uint32_t checksum(double *U) {
    // crc checksum adapted from Hudspith's GLU code

    uint32_t res = 0;
    for (int i = 0; i < CLENGTH; i++) {
        uint32_t *buf = (uint32_t*)(U + i);
        res += *buf + *(buf+1);
    }

    return res;
}

