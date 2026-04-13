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
    int *indices; // store lattice index for unordered file write
    double **measurements; // array of measurements, probably 2 doubles
    size_t number_of_measurements;
    size_t size;
    size_t block_size;
} measurement_list;

/*typedef struct {*/
/*    int t;*/
/*    int b;*/
/*    int number_of_sources;*/
/*    double *data; // array of measurements, probably 2 doubles*/
/*} averaged_data;*/

// create a struct containing dynamic list of measurements
measurement_list *create_full_measurement_list(size_t block_size);
// delete the list
void delete_measurement_list(measurement_list *list);
void source_grid(int *pos, int source, int num_sources);
measurement_list *read_data_from_file_in_order(char **inputs);
double *average_sources(measurement_list *list);
void write_data_to_tsv_file(double *data, char **inputs);

int main(int argc, char *argv[])
{
    if (argc < 10 || argc > 10) {
        printf("Input should look like: \n");
        printf("%s n ratio num_sources num_used_sources b t config_number"
               " binary_file data_dir\n", argv[0]);
        exit(1);
    }

    int n = atoi(argv[1]);
    double ratio = atof(argv[2]);
    int num_sources = atoi(argv[3]);
    int num_used_sources = atoi(argv[4]);
    int b = atoi(argv[5]);
    int t = atoi(argv[6]);
    char *config_number = argv[7];
    char *binary_path = argv[8];
    char *data_dir = argv[9];

    double averaged_data[2] = {0};
    measurement_list *ordered_data = read_data_from_file_in_order(argv);

    int rem = num_sources%num_used_sources;
    if (rem != 0) {
        printf("I can't handle cases where the number of used sources doesn't divide the number of sources with zero remainder for now\n");
        exit(0);
        // printf("Used sources does not divide number of sources cleanly\n");
        // printf("Adding %d to number of used sources\n");
        // num_used_sources += rem;
    }

    int stride = num_sources/num_used_sources;

    for (int s = 0; s<num_sources; s+=stride) {
        printf("%d\n", s);
        int lattice_s;
        source_grid(&lattice_s, s, num_sources);
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
        averaged_data[i] /= num_used_sources;
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
        list->indices = malloc(sizeof(int)*block_size);
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
    char *data_path = inputs[8];

    memset(list->indices, -1, LVOLUME*sizeof(int));

    FILE *data_file = fopen(data_path, "rb");
    if (data_file != NULL) {
        size_t num_read;
        for (int i=0; i<LVOLUME; i++) {
            num_read = 0;
            int index;
            double measurement[2];
            num_read += fread(&index, sizeof(int), 1, data_file);
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

    // code to write to file
    char data_path[200];
    sprintf(data_path, "%s/n%ss%sus%sc%s.tsv", 
            inputs[9], inputs[1], inputs[3], inputs[4], inputs[7]);


    struct stat buffer;
    int result = stat(data_path, &buffer);

    FILE *tsv_file;
    tsv_file = fopen(data_path, "a");

    if (result != 0) {
        fprintf(tsv_file, "# r\tb\tt\tnum.real\tnum.imag\t\n");
    }

    fprintf(tsv_file, "%s\t%s\t%s\t%lf\t%lf\t\n", 
            inputs[2], inputs[4], inputs[5], data[0], data[1]);

    fclose(tsv_file);
    
}











