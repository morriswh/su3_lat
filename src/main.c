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

enum dim{mu_x, mu_y, mu_z, mu_t};

int main(int argc, char *argv[]) {

    // omp_set_num_threads(16);

    struct timespec start, finish;
    double elapsed;
    
    if (argc == 1) {
        printf("Input should look like: \n");
        printf("%s config_number config_path\n", argv[0]);
        exit(1);
    }

    char *config_number = argv[1];
    char *config = argv[2];

    int header_length = 0;

    int path_length;
    char *in_file;
    int i = 0;

    // for in_file
    path_length = snprintf(NULL, 0, config);
    in_file = (char*)malloc(path_length+1);
    sprintf(in_file, config);

    double *U;
    U = (double*)calloc(CLENGTH, sizeof(double));

    char *header;
    header = (char*)malloc(600*sizeof(char));

    double tcharge[X];
    memset(tcharge, 0, sizeof(tcharge));

    read_nersc(U, in_file, header, &header_length);

    
    clock_gettime(CLOCK_MONOTONIC, &start);
    top_charge(tcharge, U);
    clock_gettime(CLOCK_MONOTONIC, &finish);
    printf("%lf\n", tcharge[0]);
    printf("%lf\n", tcharge[1]);

    free(U);
    free(header);
    free(in_file);
    
    
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("Elapsed time: %lfs\n", elapsed);

    return 0;
}
