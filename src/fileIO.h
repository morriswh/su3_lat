#ifndef FILEIO_H
#define FILEIO_H

int32_t read_nersc(double *U, char *config, char *header, int32_t *header_length);
void write_nersc(double *U, char *out_file, char *header, int32_t header_length);
void find_nersc_value(char *value, char *label, int32_t label_length, char *header, int32_t header_length);
void make_nersc_header(char *header, int32_t *header_length, double *U, char *label);

int32_t read_scidac(double *U, char *config, char *header, int32_t *header_length);
/*void write_scidac(double *U, char *out_file, char *header, int64_t header_length);*/
/*void make_scidac_header(char *header, int64_t *header_length, double *U, char *label);*/

uint32_t checksum(double *U);

#endif // !FILEIO_H
