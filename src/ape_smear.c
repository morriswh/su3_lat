#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <inttypes.h>
// #include <gmp.h>
// #include <mpfr.h>
#include <math.h>
#include <time.h>
#include "sizes.h"
#include "su3lib.h"
#include "apelib.h"

#define LENGTH(X) (sizeof X / sizeof X[0])

int main(int argc, char *argv[]) {

    // omp_set_num_threads(16);

    struct timespec start, now, finish;
    double elapsed;
    
    
    if (argc < 4) {
        printf("Input should look like: \n");
        printf("%s config_number config_path out_dir\n", argv[0]);
        exit(1);
    }

    char *config_number = argv[1];
    char *config = argv[2];
    char *out_dir = argv[3];

    int header_length = 0;
    double alpha = 0.3;
    int smear_dir[] = { 0, 1, 1, 0 };
    int smearing_steps = 10;

    int path_length;
    char *in_file, *out_file;
    int i = 0;
    int k = 0;

    // for in_file
    path_length = snprintf(NULL, 0, config);
    in_file = (char*)malloc(path_length+1);
    sprintf(in_file, config);

    // for out_file
    path_length = snprintf(NULL, 0, "%s/T64L32_E5_ape%d_a0p3_NERSC.%s", out_dir, i+1, config_number);
    out_file = (char*)malloc(path_length+1);
    sprintf(out_file, "%s/T64L32_E5_ape%d_a0p3_NERSC.%s", out_dir, i+1, config_number);

    struct stat buffer;
    // while (stat(out_file, &buffer) == 0 && i < smearing_steps) {
    for (i = 0; i < smearing_steps; i++) {

        path_length = snprintf(NULL, 0, "%s/T64L32_E5_ape%d_a0p3_NERSC.%s", out_dir, i+1, config_number);
        out_file = (char*)realloc(out_file, path_length+1);
        sprintf(out_file, "%s/T64L32_E5_ape%d_a0p3_NERSC.%s", out_dir, i+1, config_number);

        if ( stat(out_file, &buffer) == 0 ) { // if exists
            k = i + 1;
            in_file = (char*)realloc(in_file, path_length+1);
            memcpy(in_file, out_file, sizeof(char)*(path_length+1));
        }

    }

    double *U, *V;
    U = (double*)calloc(CLENGTH, sizeof(double));
    V = (double*)calloc(CLENGTH, sizeof(double));
    bzero(U, CLENGTH*sizeof(double));

    char *header;
    header = (char*)malloc(600*sizeof(char));
    bzero(header, 600*sizeof(char));

    read_nersc(U, in_file, header, &header_length);

    clock_gettime(CLOCK_MONOTONIC, &start);

    double top_charge_array[smearing_steps+1];
    double top_charge_complex[2];
    top_charge(top_charge_complex, U);
    top_charge_array[k] = top_charge_complex[0];
    clock_gettime(CLOCK_MONOTONIC, &now);
    elapsed = (now.tv_sec - start.tv_sec);
    elapsed += (now.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("[ %lf ] Top. charge: %lf\n", elapsed, top_charge_complex[0]);

    int j;
    for (j = k; j < smearing_steps; j++) {
        ape_smear(V, U, alpha, smear_dir);
        top_charge(top_charge_complex, V);
        top_charge_array[j+1] = top_charge_complex[0];
        clock_gettime(CLOCK_MONOTONIC, &now);
        elapsed = (now.tv_sec - start.tv_sec);
        elapsed += (now.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("[ %lf ] Top. charge: %lf\n", elapsed, top_charge_complex[0]);
        if ( j%50 == 49 ) {
            path_length = snprintf(NULL, 0, "%s/T64L32_E5_ape%d_a0p3_NERSC.%s", out_dir, j+1, config_number);
            out_file = (char*)realloc(out_file, path_length+1);
            make_header(header, &header_length, V);
            sprintf(out_file, "%s/T64L32_E5_ape%d_a0p3_NERSC.%s", out_dir, j+1, config_number);
            write_nersc(V, out_file, header, header_length);
            // double aplaq, splaq, tplaq;
            // all_plaquettes(V, &aplaq, &splaq, &tplaq);
            // printf( "{Plaquette} %1.15f {Spatial} %1.15f "
	           //   "{Temporal} %1.15f \n\n" , aplaq , splaq , tplaq ) ; 
        }
        memcpy(U, V, sizeof(double)*CLENGTH);
    }

    // Write out top charges

    path_length = snprintf(NULL, 0, "%s/top_charge_s%d.%s", out_dir, k+1, config_number);
    char *tc_out;
    tc_out = (char*)malloc(sizeof(char)*(path_length+1));
    sprintf(tc_out, "%s/top_charge_s%d.%s", out_dir, k, config_number);

    FILE *tc_file = fopen(tc_out, "wb");
    fwrite(top_charge_array, sizeof(double), LENGTH(top_charge_array), tc_file);
    fclose(tc_file);

    free(U);
    free(V);
    free(header);
    free(in_file);
    free(out_file);
    free(tc_out);

    
    clock_gettime(CLOCK_MONOTONIC, &finish);
    
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("Elapsed time: %lfs\n", elapsed);

    return 0;
}
