#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <sys/resource.h>
// #include <unistd.h>
#include <inttypes.h>
// #include <gmp.h>
// #include <mpfr.h>
#include <math.h>
#include <time.h>
#include "sizes.h"
#include "su3lib.h"

typedef struct {
    uint32_t *indices; // store lattice index for unordered file write
    double **measurements; // array of measurements, probably 6 doubles
    size_t number_of_measurements;
    size_t size;
    size_t block_size;
} measurement_list;

/*typedef struct {*/
/*    int t;*/
/*    int b;*/
/*    int number_of_sources;*/
/*    double *data; // array of measurements, probably 6 doubles*/
/*} averaged_data;*/

// create a struct containing dynamic list of measurements
measurement_list *create_full_measurement_list(size_t block_size);
// delete the list
void delete_measurement_list(measurement_list *list);
void source_grid(int *pos, int source, int num_sources);
measurement_list *read_data_from_file_in_order(char **inputs);
double *average_sources(measurement_list *list);
void write_data_to_tsv_file(double *data, char **inputs);
double getPhi(int *shape);
void stoia(int *output, char *input);

int main(int argc, char *argv[])
{
    if (argc < 7 || argc > 7) {
        printf("Input should look like: \n");
        printf("%s num_sources b shape config_number"
               " binary_file data_dir\n", argv[0]);
        exit(1);
    }

    int numSources = atoi(argv[1]);
    int b = atoi(argv[2]);
    int shape[4];
    stoia(shape, argv[3]);

    char *config_number = argv[4];
    char *binary_path = argv[5];
    char *data_dir = argv[6];

    double averaged_data[2] = {0};
    measurement_list *ordered_data = read_data_from_file_in_order(argv);

    for (int s = 0; s<numSources; s++) {
        int lattice_s;
        source_grid(&lattice_s, s, numSources);
        if (ordered_data->indices[lattice_s] == -1) {
            printf("missing data!\n");
            exit(1);
        }
        for (int i = 0; i<2; i++) {
            averaged_data[i] += ordered_data->measurements[lattice_s][i];
        }
    }

    delete_measurement_list(ordered_data);

    for (int i = 0; i<2; i++) {
        averaged_data[i] /= numSources;
        /*printf("%lf ", averaged_data[i]);*/
    }
    printf("\n");

    write_data_to_tsv_file(averaged_data, argv);


    return EXIT_SUCCESS;
}

measurement_list *create_full_measurement_list(size_t block_size) {
    measurement_list *list = malloc(sizeof(measurement_list));
    if (list != NULL) {
        list->number_of_measurements = block_size;
        list->size = block_size;
        list->block_size = block_size;
        list->indices = malloc(sizeof(uint32_t)*block_size);
        list->measurements = malloc(sizeof(double*)*block_size);
        if (list->indices == NULL || list->measurements == NULL) {
            free(list);
            return NULL;
        }
    }
    return list;
} 

void delete_measurement_list(measurement_list *list) {
    for (size_t i = 0; i < list->number_of_measurements; i++) {
        if (list->indices[i] != -1) {
            free(list->measurements[i]);
        }
    }
    free(list->indices);
    free(list->measurements);
    free(list);
}

void source_grid(int *pos, int source, int num_sources) {

    const int little_L = (int)(ceil(pow((num_sources/2),0.25)));
    const int little_T = 2*little_L;
    const int fact_L = L/little_L;
    const int fact_T = T/little_T;

    int src[4];

    src[0] = fact_L*(source%little_L);
    src[1] = fact_L*((source/little_L)%little_L);
    src[2] = fact_L*(((source/little_L)/little_L)%little_L);
    src[3] = fact_T*((((source/little_L)/little_L)/little_L)%little_T);

    *pos = src[3]*LCUBE + src[2]*LSQUARE + src[1]*L + src[0];
}

measurement_list *read_data_from_file_in_order(char **inputs) {

    measurement_list *list = create_full_measurement_list(LVOLUME);
    char *data_path = inputs[5];

    memset(list->indices, -1, LVOLUME*sizeof(uint32_t));

    FILE *data_file = fopen(data_path, "rb");
    if (data_file != NULL) {
        size_t num_read;
        for (int i=0; i<LVOLUME; i++) {
            num_read = 0;
            uint32_t index;
            double measurement[2];
            num_read += fread(&index, sizeof(uint32_t), 1, data_file);
            num_read += fread(measurement, sizeof(double), 2, data_file);
            if (num_read == 3) {
                list->indices[index] = index;
                list->measurements[index] = malloc(2*sizeof(double));
                memcpy(list->measurements[index],
                       measurement,2*sizeof(double));
            } else {
                break;
            }
        }

        fclose(data_file);
    }

    return list;
}

void write_data_to_tsv_file(double *data, char **inputs) {

    int shape[4];
    stoia(shape, inputs[3]);

    double phi = getPhi(shape);

    // code to write to file
    char data_path[200];
    sprintf(data_path, "%s/s%sc%s.tsv", 
            inputs[6], inputs[1], inputs[4]);


    struct stat buffer;
    int result = stat(data_path, &buffer);

    FILE *tsv_file;
    tsv_file = fopen(data_path, "a");

    if (result != 0) {
        fprintf(tsv_file, "# b\tphi\tleg1\tleg2\tleg3\tleg4\tbutterfly_real\tbutterfly_imag\t\n");
    }

    fprintf(tsv_file, "%s\t%lf\t%d\t%d\t%d\t%d\t%lf\t%lf\t\n", 
            inputs[2], 
            phi, shape[0], shape[1], shape[2], shape[3], data[0], data[1]);

    fclose(tsv_file);
    
}

double getPhi(int *shape) {
    double dotAB = (double) (shape[0]*shape[2] + shape[1]*shape[3]);
    double magA = (double) sqrt(shape[0]*shape[0] + shape[1]*shape[1]);
    double magB = (double) sqrt(shape[2]*shape[2] + shape[3]*shape[3]);
    double ratio = dotAB/(magA*magB);
    return acos(ratio);
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












