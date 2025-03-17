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

void stoia(int *output, char *input);
void swap(double **left, double **right);
void compute_and_write_topological_charge_density(double *U, char *out_dir,
                                                  size_t step, char *alpha_str,
                                                  char *config_number);

int main(int argc, char *argv[]) {

    // omp_set_num_threads(16);

    
    if (argc < 9 || argc > 9) {
        printf("Input should look like: \n");
        printf("%s config_number config_path out_dir alpha smear_directions " 
               "smearing_steps keep_interval tcharge_keep_interval\n", argv[0]);
        exit(1);
    }

    char *config_number = argv[1];
    char *config = argv[2];
    char *out_dir = argv[3];
    char *alpha_str = argv[4];
    double alpha = atof(argv[4]); //0.3;
    int smear_dir[4]; // = { 0, 1, 1, 0 }; 
    stoia(smear_dir, argv[5]);
    int smearing_steps = atoi(argv[6]); //= 400;
    int keep_interval = atoi(argv[7]);
    int tc_keep_interval = atoi(argv[8]);

    // replace the '.' with 'p' in alpha for file names
    size_t a = 0;
    while (alpha_str[a] != '\0') {
        if (alpha_str[a] == '.') {
            alpha_str[a] = 'p';
            break;
        }
        a++;
    }

    int header_length = 0;

    int path_length;
    char *in_file, *out_file;
    size_t i = 0;
    size_t k = 0;

    // for in_file
    path_length = snprintf(NULL, 0, config);
    in_file = (char*)malloc(path_length+1);
    sprintf(in_file, config);

    // for out_file
    path_length = snprintf(NULL, 0, "%s/T64L32_E5_ape%d_a%s_NERSC.%s", 
                           out_dir, i+1, alpha_str, config_number);
    out_file = (char*)malloc(path_length+1);
    sprintf(out_file, "%s/T64L32_E5_ape%d_a%s_NERSC.%s", 
            out_dir, i+1, alpha_str, config_number);

    struct stat buffer;
    // while (stat(out_file, &buffer) == 0 && i < smearing_steps) {
    for (i = 0; i < smearing_steps; i++) {

        path_length = snprintf(NULL, 0, "%s/T64L32_E5_ape%d_a%s_NERSC.%s", 
                               out_dir, i+1, alpha_str, config_number);
        out_file = (char*)realloc(out_file, path_length+1);
        sprintf(out_file, "%s/T64L32_E5_ape%d_a%s_NERSC.%s", 
                out_dir, i+1, alpha_str, config_number);

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

    char read_cksm[20];
    char label[] = "CHECKSUM\0";
    find_value(read_cksm, label, 8, header, header_length);
    printf("Read checksum from header = %s\n", read_cksm);
    uint32_t read_cks = (uint32_t)strtol(read_cksm, NULL, 16);

    uint32_t cks = checksum(U);
    char cksm[20];
    snprintf(cksm, 20, "%x", cks);
    printf("Computed checksum = %s\n", cksm);

    uint32_t cks_diff = read_cks-cks;
    printf("Difference = %u\n", cks_diff);
    if (cks_diff != 0) {
        printf("checksum failed!\n");
        exit(1);
    }

    double begin = omp_get_wtime();

    double top_charge_array[smearing_steps+1];
    memset(top_charge_array, 0, sizeof(double)*(smearing_steps+1));
    double top_charge_complex[2];
    top_charge(top_charge_complex, U);
    top_charge_array[k] = top_charge_complex[0];

    double now = omp_get_wtime();
    double elapsed = now - begin;
    printf("[ %lf ] Top. charge: %.15lf\n", elapsed, top_charge_complex[0]);

    double *W1 = U;
    double *W2 = V;

    for (size_t j = k; j < smearing_steps; j++) {
        printf("Before smearing\n");
        ape_smear_fast(W2, W1, alpha, smear_dir);
        printf("Successfully smeared\n");

        if ( j%tc_keep_interval == tc_keep_interval-1 ) {
            top_charge(top_charge_complex, W2);
            top_charge_array[j+1] = top_charge_complex[0];

            now = omp_get_wtime();
            elapsed = now - begin;
            printf("[ %lf ] Top. charge: %.15lf\n", 
                   elapsed, top_charge_complex[0]);

            compute_and_write_topological_charge_density(W2, out_dir, 
                                                         j+1, alpha_str, 
                                                         config_number);
        }

        printf("Computed topological charge\n");

        if ( j%keep_interval == keep_interval-1 ) {
            path_length = snprintf(NULL, 0, "%s/T64L32_E5_ape%d_a%s_NERSC.%s", 
                                   out_dir, j+1, alpha_str, config_number);
            out_file = (char*)realloc(out_file, path_length+1);
            char label[30];
            snprintf(label, 30*sizeof(char), "T%dL%d_ape_a%s_dir%d%d%d%d", 
                     T, L, alpha_str, 
                     smear_dir[0], smear_dir[1], smear_dir[2], smear_dir[3]);
            make_header(header, &header_length, V, label);
            sprintf(out_file, "%s/T64L32_E5_ape%d_a%s_NERSC.%s", 
                    out_dir, j+1, alpha_str, config_number);
            write_nersc(V, out_file, header, header_length);
        }

        /*memcpy(U, V, sizeof(double)*CLENGTH);*/
        swap(&W1, &W2);
    }
    
    double end = omp_get_wtime();
    double time_spent = end - begin;
    printf("Elapsed time: %lfs\n", time_spent);

    // Write out top charges

    path_length = snprintf(NULL, 0,  "%s/top_charge_s%d_ape_a%s_NERSC.%s", 
                           out_dir, k, alpha_str, config_number);
    char *tc_out;
    tc_out = (char*)malloc(sizeof(char)*(path_length+1));
    sprintf(tc_out, "%s/top_charge_s%d_ape_a%s_NERSC.%s", 
            out_dir, k, alpha_str, config_number);

    FILE *tc_file = fopen(tc_out, "wb");
    fwrite(top_charge_array, sizeof(double), LENGTH(top_charge_array), tc_file);
    fclose(tc_file);

    free(U);
    free(V);
    free(header);
    free(in_file);
    free(out_file);
    free(tc_out);

    return 0;
}

void stoia(int *output, char *input) {

    int index = 0;

    int start = 0;
    int end = -1;
    do {
        end++;
        char tmps[5];
        if (input[end] == ' ' || input[end] == '\0') {
            strncpy(tmps, input + start, sizeof(char)*(end-start));
            tmps[end-start] = '\0';
            output[index] = atoi(tmps);
            index++;
            start = end+1;
        }
    } while (input[end] != '\0');

}

void swap(double **left, double **right) {

    double *swap = *left;
    *left = *right;
    *right = swap;
}

void compute_and_write_topological_charge_density(double *U, char *out_dir,
                                                  size_t step, char *alpha_str,
                                                  char *config_number) {

    // first perform the computation

    double *top_charge_density_array;
    top_charge_density_array = (double*)calloc(LVOLUME, sizeof(double));

    #pragma omp parallel for
    for (size_t spos = 0; spos < LVOLUME; spos++) {
        double tcharge_density[2] = {0};
        top_charge_density(tcharge_density, U, spos);
        memcpy(&top_charge_density_array[spos],
               &tcharge_density[0],sizeof(double));
    }

    // next prepare file name

    int path_length = snprintf(NULL, 0, 
                               "%s/top_charge_density_st%d_ape_a%s_NERSC.%s", 
                               out_dir, step, alpha_str, config_number);
    char *tcd_out;
    tcd_out = (char*)malloc(sizeof(char)*(path_length+1));
    sprintf(tcd_out, "%s/top_charge_density_st%d_ape_a%s_NERSC.%s", 
            out_dir, step, alpha_str, config_number);


    // then write to file

    FILE *tcd_file = fopen(tcd_out, "wb");
    fwrite(top_charge_density_array, sizeof(double), 
           LVOLUME, tcd_file);
    fclose(tcd_file);

    // finally free memory

    free(top_charge_density_array);
    free(tcd_out);

}

