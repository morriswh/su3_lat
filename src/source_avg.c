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
/*    long int t;*/
/*    long int b;*/
/*    long int number_of_sources;*/
/*    double *data; // array of measurements, probably 6 doubles*/
/*} averaged_data;*/

// create a struct containing dynamic list of measurements
measurement_list *create_full_measurement_list(size_t block_size);
// delete the list
void delete_measurement_list(measurement_list *list);
void source_grid(long int *pos, long int source, long int num_sources);
measurement_list *read_data_from_file_in_order(char **inputs);
double *average_sources(measurement_list *list);
void write_data_to_tsv_file(double *data, char **inputs);

long int main(long int argc, char *argv[])
{
    if (argc < 9 || argc > 9) {
        printf("Input should look like: \n");
        printf("%s n ratio num_sources b t config_number"
               " binary_file data_dir\n", argv[0]);
        exit(1);
    }

    long int n = atoi(argv[1]);
    double ratio = atof(argv[2]);
    long int numSources = atoi(argv[3]);
    long int b = atoi(argv[4]);
    long int t = atoi(argv[5]);
    char *config_number = argv[6];
    char *binary_path = argv[7];
    char *data_dir = argv[8];

    double averaged_data[6] = {0};
    measurement_list *ordered_data = read_data_from_file_in_order(argv);

    for (long int s = 0; s<numSources; s++) {
        long int lattice_s;
        source_grid(&lattice_s, s, numSources);
        if (ordered_data->indices[lattice_s] == -1) {
            printf("missing data!\n");
            exit(1);
        }
        for (long int i = 0; i<6; i++) {
            averaged_data[i] += ordered_data->measurements[lattice_s][i];
        }
    }

    delete_measurement_list(ordered_data);

    for (long int i = 0; i<6; i++) {
        averaged_data[i] /= numSources;
        printf("%lf ", averaged_data[i]);
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

void source_grid(long int *pos, long int source, long int num_sources) {

    const long int little_L = (int)(ceil(pow((num_sources/2),0.25)));
    const long int little_T = 2*little_L;
    const long int fact_L = L/little_L;
    const long int fact_T = T/little_T;

    long int src[4];

    src[0] = fact_L*(source%little_L);
    src[1] = fact_L*((source/little_L)%little_L);
    src[2] = fact_L*(((source/little_L)/little_L)%little_L);
    src[3] = fact_T*((((source/little_L)/little_L)/little_L)%little_T);

    *pos = src[3]*LCUBE + src[2]*LSQUARE + src[1]*L + src[0];
}

measurement_list *read_data_from_file_in_order(char **inputs) {

    measurement_list *list = create_full_measurement_list(LVOLUME);
    char *data_path = inputs[7];

    memset(list->indices, -1, LVOLUME*sizeof(uint32_t));

    FILE *data_file = fopen(data_path, "rb");
    if (data_file != NULL) {
        size_t num_read;
        for (long int i=0; i<LVOLUME; i++) {
            num_read = 0;
            uint32_t index;
            double measurement[6];
            num_read += fread(&index, sizeof(uint32_t), 1, data_file);
            num_read += fread(measurement, sizeof(double), 6, data_file);
            if (num_read == 7) {
                list->indices[index] = index;
                list->measurements[index] = malloc(6*sizeof(double));
                memcpy(list->measurements[index],
                       measurement,6*sizeof(double));
            } else {
                break;
            }
        }

        fclose(data_file);
    }

    return list;
}

void write_data_to_tsv_file(double *data, char **inputs) {

    // code to write to file
    char data_path[200];
    sprintf(data_path, "%s/n%ss%sc%s.tsv", 
            inputs[8], inputs[1], inputs[3], inputs[6]);


    struct stat buffer;
    long int result = stat(data_path, &buffer);

    FILE *tsv_file;
    tsv_file = fopen(data_path, "a");

    if (result != 0) {
        fprintf(tsv_file, "# r\tb\tt\tnum.real\tnum.imag\t"
                "denA.real\tdenA.imag\tdenB.real\tdenB.imag\t\n");
    }

    fprintf(tsv_file, "%s\t%s\t%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t\n", 
            inputs[2], inputs[4], inputs[5], 
            data[0], data[1], data[2], data[3], data[4], data[5]);

    fclose(tsv_file);
    
}











