#ifndef SU3LIB_H
#define SU3LIB_H

// declaration of functions
void get_index(long int i); // returns index from position in config array
void get_link(double *link, double *U, long int pos); // 12 reals for su3
void get_site(double *site, double *S, long int spos); // 12 reals for su3
long int ipow(long int base, long int exp); // exponentiation by squaring
long int mod(long int value, long int period); // modulo function
void matmul(double *z, double *x, double *y); // product of two su(3) matrices
void cmmul(double *z, double c, double *x); // scalar times su(3) matrix
void adjoint(double *xc, double *x); // get conjugate transpose
void print_link(double *link); // prlong int gauge link, needs cosmetic attention
void trace(double *result, double *matrix); // trace of su(3) matrix
void det(double *result, double *matrix); // determinant of su(3) matrix
void transform(double *Up, double *U, double *Omega); // gauge transformation
void plaquette(double *plaq, double *U, long int spos, long int mu, long int nu); // plaquette
void all_plaquettes(double *U, double *a_plaq, double *s_plaq, double *t_plaq);
void cloverleaf(double *clover, double *U, long int spos, long int mu, long int nu); // clover
void top_charge_density(double *topo_charge_density, double *U, long int spos); // at one site
void top_charge(double *topo_charge, double *U);
void directional_sum(double *lattice, long int *directions, 
                     double *out_array, long int out_array_size);
void staples(double *staple, double *U, long int spos, long int mu, long int nu); // staple
void link_trace(double *ltr, double *U);
void gauge_action(double *action, double *U);
void rnd_config(double *U);
void rnd_gauge_tr(double *Omega);

#endif // !SU3LIB_H
