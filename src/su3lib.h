#ifndef SU3LIB_H
#define SU3LIB_H

// declaration of functions
void get_index(int64_t i); // returns index from position in config array
void get_site(double *site, double *S, int64_t spos); // 12 reals for su3
int64_t ipow(int64_t base, int64_t exp); // exponentiation by squaring
void immul(double *z, double c, double *x); // imaginary scalar times su(3) matrix
void matadd(double *z, double *x, double *y); // sum of two su(3) matrices
void print_link(double *link); // print64_t gauge link, needs cosmetic attention
void det(double *result, double *matrix); // determinant of su(3) matrix
void transform(double *Up, double *U, double *Omega); // gauge transformation
void plaquette(double *plaq, double *U, int64_t spos, int32_t mu, int32_t nu); // plaquette
void all_plaquettes(double *U, double *a_plaq, double *s_plaq, double *t_plaq);
void cloverleaf(double *clover, double *U, int64_t spos, int64_t mu, int64_t nu); // clover
void top_charge_density(double *topo_charge_density, double *U, int64_t spos); // at one site
void top_charge(double *topo_charge, double *U);
void directional_sum(double *lattice, int64_t *directions, 
                     double *out_array, int64_t out_array_size);
void staples(double *staple, double *U, int64_t spos, int64_t mu, int64_t nu); // staple
void link_trace(double *ltr, double *U);
void gauge_action(double *action, double *U);
void unit_config(double *U);
void rnd_config(double *U);
void rnd_gauge_tr(double *Omega);

/* frequently used functions
 * made the below functions static inline functions because it apparently leads to increased efficiency
 */
// void matmul(double *z, double *x, double *y); // product of two su(3) matrices
// void adjoint(double *xc, double *x); // get conjugate transpose
// void trace(double *result, double *matrix); // trace of su(3) matrix
// void cmmul(double *z, double c, double *x); // scalar times su(3) matrix
// int64_t mod(int64_t value, int64_t period); // modulo function
// void get_link(double *link, double *U, int64_t pos); // 12 reals for su3

static inline void matmul(double *z, const double *x, const double *y) {
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

static inline void matmul_na(double * z, const double * x, const double * y) {
    double u[X];
    double v[X];
    memcpy(u,x,X*sizeof(double));
    memcpy(v,y,X*sizeof(double));

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            // z_ij = sum_k u_ik * conj(v_jk)
            double re = 0.0, im = 0.0;
            for (int k = 0; k < 3; k++) {
                const double ur = u[6*i + 2*k], ui = u[6*i + 2*k + 1];
                const double vr = v[6*j + 2*k], vi = v[6*j + 2*k + 1];
                // terms added one at a time, in the same order as matmul,
                // so results match adjoint + matmul euactlv
                re += ur*vr;
                re += ui*vi;
                im -= ur*vi;
                im += ui*vr;
            }
            z[6*i + 2*j]     = re;
            z[6*i + 2*j + 1] = im;
        }
    }
}

static inline void matmul_na_nc(double *restrict z, const double *restrict x, const double *restrict y) {
    //no copy
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            // z_ij = sum_k x_ik * conj(y_jk)
            double re = 0.0, im = 0.0;
            for (int k = 0; k < 3; k++) {
                const double xr = x[6*i + 2*k], xi = x[6*i + 2*k + 1];
                const double yr = y[6*j + 2*k], yi = y[6*j + 2*k + 1];
                // terms added one at a time, in the same order as matmul,
                // so results match adjoint + matmul exactly
                re += xr*yr;
                re += xi*yi;
                im -= xr*yi;
                im += xi*yr;
            }
            z[6*i + 2*j]     = re;
            z[6*i + 2*j + 1] = im;
        }
    }
}

static inline void matmul_nc(double *restrict z, const double *restrict x, const double *restrict y) {
    // no copy

    z[0]  = x[0]*y[0]   - x[1]*y[1]   + x[2]*y[6]          // m[0] + i m[1]     m[2] + i m[3]     m[4] + i m[5]
          - x[3]*y[7]   + x[4]*y[12]  - x[5]*y[13];        // 
    z[1]  = x[0]*y[1]   + x[1]*y[0]   + x[2]*y[7]          // 
          + x[3]*y[6]   + x[4]*y[13]  + x[5]*y[12];        // 
    z[2]  = x[0]*y[2]   - x[1]*y[3]   + x[2]*y[8]          // m[6] + i m[7]     m[8] + i m[9]     m[10] + i m[11]
          - x[3]*y[9]   + x[4]*y[14]  - x[5]*y[15];        // 
    z[3]  = x[0]*y[3]   + x[1]*y[2]   + x[2]*y[9]          // 
          + x[3]*y[8]   + x[4]*y[15]  + x[5]*y[14];        // 
    z[4]  = x[0]*y[4]   - x[1]*y[5]   + x[2]*y[10]         // m[12] + i m[13]   m[14] + i m[15]   m[16] + i m[17]
          - x[3]*y[11]  + x[4]*y[16]  - x[5]*y[17];
    z[5]  = x[0]*y[5]   + x[1]*y[4]   + x[2]*y[11]
          + x[3]*y[10]  + x[4]*y[17]  + x[5]*y[16];
    z[6]  = x[6]*y[0]   - x[7]*y[1]   + x[8]*y[6]
          - x[9]*y[7]   + x[10]*y[12] - x[11]*y[13];
    z[7]  = x[6]*y[1]   + x[7]*y[0]   + x[8]*y[7]
          + x[9]*y[6]   + x[10]*y[13] + x[11]*y[12];
    z[8]  = x[6]*y[2]   - x[7]*y[3]   + x[8]*y[8]
          - x[9]*y[9]   + x[10]*y[14] - x[11]*y[15];
    z[9]  = x[6]*y[3]   + x[7]*y[2]   + x[8]*y[9]
          + x[9]*y[8]   + x[10]*y[15] + x[11]*y[14];
    z[10] = x[6]*y[4]   - x[7]*y[5]   + x[8]*y[10]
          - x[9]*y[11]  + x[10]*y[16] - x[11]*y[17];
    z[11] = x[6]*y[5]   + x[7]*y[4]   + x[8]*y[11]
          + x[9]*y[10]  + x[10]*y[17] + x[11]*y[16];
    z[12] = x[12]*y[0]  - x[13]*y[1]  + x[14]*y[6]
          - x[15]*y[7]  + x[16]*y[12] - x[17]*y[13];
    z[13] = x[12]*y[1]  + x[13]*y[0]  + x[14]*y[7]
          + x[15]*y[6]  + x[16]*y[13] + x[17]*y[12];
    z[14] = x[12]*y[2]  - x[13]*y[3]  + x[14]*y[8]
          - x[15]*y[9]  + x[16]*y[14] - x[17]*y[15];
    z[15] = x[12]*y[3]  + x[13]*y[2]  + x[14]*y[9]
          + x[15]*y[8]  + x[16]*y[15] + x[17]*y[14];
    z[16] = x[12]*y[4]  - x[13]*y[5]  + x[14]*y[10]
          - x[15]*y[11] + x[16]*y[16] - x[17]*y[17];
    z[17] = x[12]*y[5]  + x[13]*y[4]  + x[14]*y[11]
          + x[15]*y[10] + x[16]*y[17] + x[17]*y[16];
}

static inline void adjoint(double *xc, double *x) {
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

static inline void trace(double *result, double *matrix) {
    result[0] = matrix[0] + matrix[8] + matrix[16];
    result[1] = matrix[1] + matrix[9] + matrix[17];
}

static inline void cmmul(double *z, double c, double *x) {
    for (int32_t i = 0; i < X; ++i) {
        z[i] = c*x[i];
    }
}

static inline int64_t mod(int64_t value, int64_t period) {
    // avoid taking unnecessary operations
    if (value >= 0 && value < period) {
        return value;
    }
    else {
        return (value%period + period)%period;
    }
    // return (value%period + period)%period;
}

static inline void get_link(double *link, double *U, int64_t pos) {
    // pos = mod(pos, LLENGTH);
    pos *= X;
    memcpy(link,&U[pos],sizeof(double)*X);
}

// static inline const double *link_ptr(const double *U, int64_t pos) {
//     return U[pos];
// }

#endif // !SU3LIB_H
