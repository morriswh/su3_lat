#ifndef SU3LIB_H
#define SU3LIB_H

// declaration of functions
void get_index(int64_t i); // returns index from position in config array
void get_link(double *link, double *U, int64_t pos); // 12 reals for su3
void get_site(double *site, double *S, int64_t spos); // 12 reals for su3
int64_t ipow(int64_t base, int64_t exp); // exponentiation by squaring
int64_t mod(int64_t value, int64_t period); // modulo function
void matmul(double *z, double *x, double *y); // product of two su(3) matrices
void cmmul(double *z, double c, double *x); // scalar times su(3) matrix
void immul(double *z, double c, double *x); // imaginary scalar times su(3) matrix
void matadd(double *z, double *x, double *y); // sum of two su(3) matrices
void adjoint(double *xc, double *x); // get conjugate transpose
void print_link(double *link); // print64_t gauge link, needs cosmetic attention
void trace(double *result, double *matrix); // trace of su(3) matrix
void det(double *result, double *matrix); // determinant of su(3) matrix
void transform(double *Up, double *U, double *Omega); // gauge transformation
void plaquette(double *plaq, double *U, int64_t spos, int64_t mu, int64_t nu); // plaquette
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

#endif // !SU3LIB_H
