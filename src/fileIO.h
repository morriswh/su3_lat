#ifndef FILEIO_H
#define FILEIO_H

int read_nersc(double *U, char *config, char *header, int *header_length);
void write_nersc(double *U, char *out_file, char *header, int header_length);
void find_nersc_value(char *value, char *label, int label_length, char *header, int header_length);
void make_nersc_header(char *header, int *header_length, double *U, char *label);

int read_scidac(double *U, char *config, char *header, int *header_length);
/*void write_scidac(double *U, char *out_file, char *header, int header_length);*/
/*void make_scidac_header(char *header, int *header_length, double *U, char *label);*/

int checksum(double *U);

#endif // !FILEIO_H
