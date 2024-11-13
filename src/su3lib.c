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

// #define LENGTH(X) (sizeof X / sizeof X[0])

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
    for (int i = 0; i < X; i++) {
        u[i] = x[i];
        v[i] = y[i];
    }

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
    printf("np.array([[(%lf)+1j*(%lf),(%lf)+1j*(%lf),(%lf)+1j*(%lf)], [(%lf)+1j*(%lf),(%lf)+1j*(%lf),(%lf)+1j*(%lf)],[(%lf)+1j*(%lf),(%lf)+1j*(%lf),(%lf)+1j*(%lf)]]\n",
           link[0], link[1], link[2], link[3], link[4], link[5],
           link[6], link[7], link[8], link[9], link[10], link[11],
           link[12], link[13], link[14], link[15], link[16], link[17]);
}

void trace(double *result, double *matrix) {
    result[0] = matrix[0] + matrix[8] + matrix[16];
    result[1] = matrix[1] + matrix[9] + matrix[17];
}

void plaquette(double *plaq, double *U, int spos, int mu, int nu) {
    // spos is in LVOLUME now
    int offs_br, offs_tl, period_br, period_tl;

    if (mu == 0) {
        offs_br = 1;
        period_br = L;
    }
    else if (mu == 1) {
        offs_br = L;
        period_br = LSQUARE;
    }
    else if (mu == 2) {
        offs_br = LSQUARE;
        period_br = LCUBE;
    }
    else {
        offs_br = LCUBE;
        period_br = LVOLUME;
    }

    if (nu == 0) {
        offs_tl = 1;
        period_tl = L;
    }
    else if (nu == 1) {
        offs_tl = L;
        period_tl = LSQUARE;
    }
    else if (nu == 2) {
        offs_tl = LSQUARE;
        period_tl = LCUBE;
    }
    else {
        offs_tl = LCUBE;
        period_tl = LVOLUME;
    }

    int pos = D*spos + mu;
    double temp[X];
    get_link(plaq, U, mod(pos, LLENGTH));
    // matmul(plaq,plaq,temp);

    pos = D*(spos - mod(spos, period_br) + mod(spos + offs_br, period_br)) + nu;
    get_link(temp, U, mod(pos, LLENGTH));
    matmul(plaq,plaq,temp);

    // pos = D*(spos+offs_tl) + mu;
    pos = D*(spos - mod(spos, period_tl) + mod(spos + offs_tl, period_tl)) + mu;
    get_link(temp, U, mod(pos, LLENGTH));
    adjoint(temp, temp);
    matmul(plaq,plaq,temp);

    pos = D*spos + nu;
    get_link(temp, U, mod(pos, LLENGTH));
    adjoint(temp,temp);
    matmul(plaq,plaq,temp);
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
    double sum = 0;
    double tr[2];
    int i = 0, j = 0, k = 0;
    // I don't understand 'reduction' yet, but apparently it is needed
    #pragma omp parallel for reduction(+:sum) private(tr, plaq, j, k)
    for (i = 0; i < LVOLUME; i++) {
        for (j = 0; j < D; j++) {
            for (k = j+1; k < D; k++) {
                plaquette(plaq, U, i, j, k);
                trace(tr, plaq);
                sum += tr[0];
            }
        }
    }
    *action = 2*beta*(3*T*L - sum);
}

void read_nersc(double *U, char *config) {
    // need to figure out the issue with fread and doubles for the NERSC config
    // for now the union method works I guess
    union {
        long int retriever;
        double element;
    } package; 
    FILE *config_file = fopen(config, "rb");
    // read header, should makethis more general maybe 
    char header[50];
    for (int i = 0; i < 22; i++) {
        fgets(header, 50, config_file);
        printf("%s", header);
    }
    // read in configuration data
    for (int i = 0; i != CLENGTH; i++) {
        fread(&package.retriever, sizeof(package.retriever), 1, config_file);
        package.retriever = be64toh(package.retriever);
        U[i] = package.element;
    }
    fclose(config_file);
}


// void read_nersc(double *U) {
//     // need to figure out the issue with fread and doubles for the NERSC config
//     // for now the union method works I guess
//     union {
//         long int retriever[CLENGTH];
//         double U[CLENGTH];
//     } package; 
//     FILE *config = fopen("thyp2.1_gauge_nohead", "rb");
//     // const size_t ret_code = 
//     fread(package.retriever, sizeof(package.retriever[0]), CLENGTH, config);
//     for (int i = 0; i != CLENGTH; i++) {
//         package.retriever[i] = be64toh(package.retriever[i]);
//         U[i] = package.U[i];
//     }
// }
