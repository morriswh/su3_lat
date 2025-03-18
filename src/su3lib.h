#ifndef SU3LIB_H
#define SU3LIB_H

// declaration of functions
void get_index(int i); // returns index from position in config array
void get_link(double *link, double *U, int pos); // 12 reals for su3
void get_site(double *site, double *S, int spos); // 12 reals for su3
int ipow(int base, int exp); // exponentiation by squaring
int mod(int value, int period); // modulo function
void matmul(double *z, double *x, double *y); // product of two su(3) matrices
void cmmul(double *z, double c, double *x); // scalar times su(3) matrix
void adjoint(double *xc, double *x); // get conjugate transpose
void print_link(double *link); // print gauge link, needs cosmetic attention
void trace(double *result, double *matrix); // trace of su(3) matrix
void det(double *result, double *matrix); // determinant of su(3) matrix
void transform(double *Up, double *U, double *Omega); // gauge transformation
void plaquette(double *plaq, double *U, int spos, int mu, int nu); // plaquette
void all_plaquettes(double *U, double *a_plaq, double *s_plaq, double *t_plaq);
void cloverleaf(double *clover, double *U, int spos, int mu, int nu); // clover
void top_charge_density(double *topo_charge_density, double *U, int spos); // at one site
void top_charge(double *topo_charge, double *U);
void directional_sum(double *lattice, int *directions, 
                     double *out_array, int out_array_size);
void staples(double *staple, double *U, int spos, int mu, int nu); // staple
void link_trace(double *ltr, double *U);
void gauge_action(double *action, double *U);
void rnd_config(double *U);
void rnd_gauge_tr(double *Omega);
int read_nersc(double *U, char *config, char *header, int *header_length);
void write_nersc(double *U, char *out_file, char *header, int header_length);
void find_value(char *value, char *label, int label_length, char *header, int header_length);
void make_header(char *header, int *header_length, double *U, char *label);
uint32_t checksum(double *U);

#endif // !SU3LIB_H
