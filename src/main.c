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

    if (argc == 1) {
        printf("Input should look like: \n");
        printf("%s config_number config_path\n", argv[0]);
        exit(1);
    }

    struct timespec start, now, finish;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &start);

    char *config_number = argv[1];
    char *config = argv[2];
    
    
    double *U;
    U = (double*)calloc(CLENGTH, sizeof(double));
    memset(U, 0, CLENGTH*sizeof(double));

    char *header;
    header = (char*)malloc(600*sizeof(char));
    memset(header, 0, 600*sizeof(char));

    int header_length = 0;
    read_nersc(U, config, header, &header_length);

    double alpha = 0.3;

    // make coefficients
    double one_minus_alpha_matrix_sc = 1 - alpha;
    double alpha_matrix_div_six_sc = alpha/6.0;

    double one_minus_alpha_matrix[X], alpha_matrix_div_six[X];

    bzero(one_minus_alpha_matrix, sizeof(one_minus_alpha_matrix));
    one_minus_alpha_matrix[0] = 1 - alpha;
    one_minus_alpha_matrix[8] = 1 - alpha;
    one_minus_alpha_matrix[16] = 1 - alpha;

    bzero(alpha_matrix_div_six, sizeof(alpha_matrix_div_six));
    alpha_matrix_div_six[0] = alpha/6.0;
    alpha_matrix_div_six[8] = alpha/6.0;
    alpha_matrix_div_six[16] = alpha/6.0;

    double mat[X];
    double sca[X], mmu[X];
    get_link(mat, U, 0);
    matmul(mmu, one_minus_alpha_matrix, mat);
    cmmul(sca, one_minus_alpha_matrix_sc, mat);
    print_link(mmu);
    print_link(sca);




    // double aplaq, splaq, tplaq;
    //
    // all_plaquettes(U, &aplaq, &splaq, &tplaq);
    //
    // printf( "{Plaquette} %1.15f {Spatial} %1.15f "
	   //   "{Temporal} %1.15f \n\n" , aplaq , splaq , tplaq ) ; 

    free(U);
    free(header);
    
    clock_gettime(CLOCK_MONOTONIC, &finish);
    
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("Elapsed time: %lfs\n", elapsed);

    return 0;
}
