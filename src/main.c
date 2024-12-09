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

    struct timespec start, finish;
    double elapsed;
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    
    if (argc == 1) {
        printf("Input should look like: \n");
        printf("%s config_number config_path out_dir\n", argv[0]);
        exit(1);
    }

    char *config_number = argv[1];
    char *config = argv[2];
    char *out_dir = argv[3];

    int header_length = 0;
    double alpha = 0.3;
    int smearing_steps = 20;

    int path_length;
    char *in_file, *out_file;
    int i = 0;

    // for in_file
    path_length = snprintf(NULL, 0, config);
    in_file = (char*)malloc(path_length+1);
    sprintf(in_file, config);

    // for out_file
    path_length = snprintf(NULL, 0, "%s/T64L32_E5_ape%d_a0p3_NERSC.%s", out_dir, i+1, config_number);
    out_file = (char*)malloc(path_length+1);
    sprintf(out_file, "%s/T64L32_E5_ape%d_a0p3_NERSC.%s", out_dir, i+1, config_number);

    struct stat buffer;
    while (stat(out_file, &buffer) == 0 && i < smearing_steps) {
        free(in_file);
        in_file = (char*)malloc(path_length+1);
        // sprintf(in_file, "%s", out_file);
        memcpy(in_file, out_file, sizeof(char)*(path_length+1));

        i++;

        path_length = snprintf(NULL, 0, "%s/T64L32_E5_ape%d_a0p3_NERSC.%s", out_dir, i+1, config_number);
        free(out_file);
        out_file = (char*)malloc(path_length+1);
        sprintf(out_file, "%s/T64L32_E5_ape%d_a0p3_NERSC.%s", out_dir, i+1, config_number);
    }

    double *U, *V;
    U = (double*)calloc(CLENGTH, sizeof(double));
    V = (double*)calloc(CLENGTH, sizeof(double));

    char *header;
    header = (char*)malloc(600*sizeof(char));

    read_nersc(U, in_file, header, &header_length);
    int j;
    for (j = i; j < smearing_steps+1; j++) {
        ape_smear_yz(V, U, alpha);
        if ( j%4 == 3 && j > 0) {
            free(out_file);
            path_length = snprintf(NULL, 0, "%s/T64L32_E5_ape%d_a0p3_NERSC.%s", out_dir, j+1, config_number);
            out_file = (char*)malloc(path_length+1);
            sprintf(out_file, "%s/T64L32_E5_ape%d_a0p3_NERSC.%s", out_dir, j+1, config_number);
            write_nersc(V, out_file, header, header_length);
        }

        // free(in_file);
        // in_file = (char*)malloc(path_length+1);
        // memcpy(in_file, out_file, sizeof(char)*(path_length+1));

        memcpy(U, V, sizeof(double)*CLENGTH);
    }

    free(U);
    free(V);
    free(header);
    free(in_file);
    free(out_file);

    
    clock_gettime(CLOCK_MONOTONIC, &finish);
    
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("Elapsed time: %lfs\n", elapsed);

    return 0;
}
