#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <inttypes.h>
// #include <gmp.h>
// #include <mpfr.h>
#include <math.h>
#include <time.h>
#include "sizes.h"
#include "matops.h"
#include "su3lib.h"
#include "apelib.h"

#define LENGTH(X) (sizeof X / sizeof X[0])

int main(int argc, char *argv[]) {

    // omp_set_num_threads(16);

    // if (argc == 1) {
    //     printf("Input should look like: \n");
    //     printf("%s config_number config_path\n", argv[0]);
    //     exit(1);
    // }

    struct timespec start, now, finish;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &start);

    char *config_number = argv[1];
    char *config = argv[2];


    double *U;
    U = (double*)calloc(CLENGTH, sizeof(double));
    memset(U, 0, CLENGTH*sizeof(double));

    char *header;
    header = (char*)malloc(1000*sizeof(char));
    memset(header, '\0', 1000*sizeof(char));

    int header_length = 0;
    read_nersc(U, config, header, &header_length);

    make_header(header, &header_length, U);
    printf("%s", header);


    free(U);
    free(header);
    
    clock_gettime(CLOCK_MONOTONIC, &finish);
    
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("Elapsed time: %lfs\n", elapsed);

    return 0;
}
