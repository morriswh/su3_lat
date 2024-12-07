#ifndef SU3LIB_H
#define SU3LIB_H

// declaration of functions
void get_index(int i); // returns index from position in config array
void get_link(double *link, double *U, int pos); // 12 reals for su3
void get_site(double *site, double *S, int spos); // 12 reals for su3
int ipow(int base, int exp); // exponentiation by squaring
int mod(int value, int period); // modulo function
void matmul(double *z, double *x, double *y); // product of two su(3) matrices
void adjoint(double *xc, double *x); // get conjugate transpose
void print_link(double *link); // print gauge link, needs cosmetic attention
void trace(double *result, double *matrix); // trace of su(3) matrix
void transform(double *Up, double *U, double *Omega); // gauge transformation
void plaquette(double *plaq, double *U, int spos, int mu, int nu); // plaquette
void staples(double *staple, double *U, int spos, int mu, int nu); // staple
void gauge_action(double *action, double *U);
void rnd_config(double *U);
void rnd_gauge_tr(double *Omega);
void read_nersc(double *U, char *config, char *header, int *header_length);
void write_nersc(double *U, char *out_file, char *header, int header_length);

#endif // !SU3LIB_H
