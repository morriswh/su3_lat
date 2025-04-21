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
#include "fileIO.h"

#define LENGTH(X) (sizeof X / sizeof X[0])

enum dim{mu_x, mu_y, mu_z, mu_t};

// Creating the structure
typedef struct {
    uint32_t *indices; // store lattice index for unordered file write
    double **measurements; // array of measurements, probably 2 doubles
    size_t number_of_measurements;
    size_t size;
    size_t block_size;
} measurement_list;

//// declaration of functions
void segment(double *loop, double *U, int *pos, int mu, int comp);
void diagonal_segment(double *loop, double *U, int *pos, 
                      int mu_i, int mu_j, int icomp, int jcomp);
int compare(int x, int y, int xsteps, int ysteps); // for "smart" pathing
void stoia(int *output, char *input);
// create a struct containing dynamic list of measurements
measurement_list *create_measurement_list(size_t block_size);
// append a new measurement with lattice index to the list
int add_measurement_to_list(measurement_list *list, int index, 
                            double *measurement);
// delete the list
void delete_measurement_list(measurement_list *list);
void source_grid(int *pos, int source, int num_sources);
void compute_observables(double *trace_array, int spos, double *U, int b, int *shape);
double *read_configuration(char *config);
void delete_configuration(double *V);
double *random_gauge_transform(double *V);
double *copy_configuration(double *V);
int *setup_parallization(int *sizes_extent, int num_sources);
void create_info_file(char **inputs);
void write_data_to_file(measurement_list *list, char **inputs);
measurement_list *read_data_from_file(char **inputs);
uint32_t *get_existing_sources(size_t *number_of_sources, char **inputs);
int in_source_list(uint32_t source, uint32_t *source_list, size_t size);

// void build_string(char *out, char *in[], int length);

int main(int argc, char *argv[]) {

    if (argc < 7 || argc > 7) {
        printf("Input should look like: \n");
        printf("%s num_sources b shape config_number"
               " config_path data_dir\n", argv[0]);
        exit(1);
    }

    int nSources = atoi(argv[1]);
    int b = atoi(argv[2]);
    int shape[4];
    stoia(shape, argv[3]);

    char *config_number = argv[4];
    char *config = argv[5];
    char *data_dir = argv[6];

    // clock_t begin = clock();
    double begin = omp_get_wtime();

    double *V = read_configuration(config);

    create_info_file(argv); // if it doesn't already exist

    size_t existantSourceNumber; 
    uint32_t *existantSources = 
        get_existing_sources(&existantSourceNumber, argv);

    // setup parallelization
    int sizes_extent;
    int *sizes = setup_parallization(&sizes_extent, nSources);

    /*printf("%d\n", omp_get_num_threads());*/

    #pragma omp parallel for
    for (int parallel_s = 0; parallel_s < sizes_extent; parallel_s++) {

        // faster to copy the configuration for each thread, but obviously
        // uses a lot more memory
        double *U = copy_configuration(V);
        measurement_list *data_list = create_measurement_list(16);

        for (int s = sizes[parallel_s]; s < sizes[parallel_s+1]; s++) {

            // map source index to lattice index
            int lattice_s;
            source_grid(&lattice_s, s, nSources);
            // check if source has already been measured
            if (in_source_list(lattice_s, existantSources, 
                               existantSourceNumber) != 0) {
                continue;
            }

            double trace_array[2];

            // here is where the physics happens
            compute_observables(trace_array, lattice_s, U, b, shape);

            if(add_measurement_to_list(data_list, 
                                       lattice_s, trace_array) != 0) {
                printf("failed to add measurement to list!\n");
                printf("this is probably due to memory exhaustion.\n");
                exit(1);
            }

            if (data_list->number_of_measurements > (size_t)20000) {
                #pragma omp critical
                write_data_to_file(data_list, argv);
                delete_measurement_list(data_list);
                data_list = create_measurement_list(16);
            }
        }

        #pragma omp critical
        write_data_to_file(data_list, argv);
        delete_measurement_list(data_list);
        delete_configuration(U);
    }
    #pragma omp barrier

    free(sizes);
    free(existantSources);
    delete_configuration(V);

    double end = omp_get_wtime();
    double time_spent = end - begin;
    printf("Elapsed time: %lfs\n", time_spent);

    return 0;
}

void segment(double *loop, double *U, int *pos, int mu, int comp) {

    int steps = abs(comp);
    int dir = comp/abs(comp);

    int spos = *pos;

    int step = ipow(L,mu);
    int period = L*step;
    if (mu == 3) {
        period = 2*period;
    }

    double seg[X], link[X];
    memset(seg, 0, sizeof(seg));
    seg[0] = 1;  
    seg[8] = 1;  
    seg[16] = 1;

    if (dir == -1) {
        spos = spos - mod(spos, period) + mod(spos - steps*step, period);
    }

    for (int i = 0; i < steps; i++) {
        get_link(link, U, D*spos+mu);
        matmul(seg, seg, link);
        spos = spos - mod(spos, period) + mod(spos + step, period);
    }

    if (dir == -1) {
        adjoint(seg, seg);
        spos = spos - mod(spos, period) + mod(spos - steps*step, period);
    }

    // update loop and position
    matmul(loop, loop, seg);
    *pos = spos;
    
}

void diagonal_segment(double *loop, double *U, int *pos, 
                      int mu_i, int mu_j, int icomp, int jcomp) {

    // if it's just a straight line, call segment instead to avoid div0
    if (icomp == 0) {
        segment(loop, U, pos, mu_j, jcomp);
    }
    else if (jcomp == 0) {
        segment(loop, U, pos, mu_i, icomp);
    }
    else {

    int isteps = abs(icomp);
    int jsteps = abs(jcomp);
    int idir = icomp/abs(icomp);
    int jdir = jcomp/abs(jcomp);

    int spos = *pos;

    int istep = ipow(L,mu_i);
    int jstep = ipow(L,mu_j);
    int iperiod = L*istep;
    int jperiod = L*jstep;
    if (mu_i == 3) {
        iperiod = 2*iperiod;
    }
    if (mu_j == 3) {
        jperiod = 2*jperiod;
    }

    double seg[X], link[X];
    memset(seg, 0, sizeof(seg));
    seg[0] = 1;  
    seg[8] = 1;  
    seg[16] = 1;

    int i = 1;
    int j = 1;
    int mu;
    for (int k = 0; k < isteps+jsteps; k++) {
        // choose values based on result of compare function
        // want to get best approximation of diagonal line
        // need to see if there's a better
        int result = compare(i, j, isteps, jsteps);
        int mu, step, period, dir;
        if (result == 1) { // mu_j direction
            mu = mu_j;
            step = jstep;
            period = jperiod;
            dir = jdir;
            j++;
        } else { // mu_i direction
            mu = mu_i;
            step = istep;
            period = iperiod;
            dir = idir;
            i++;
        }

        if (dir == 1) {
            get_link(link, U, D*spos+mu);
            matmul(seg, seg, link);
        }
        spos = spos - mod(spos, period) + mod(spos + step*dir, period);
        if (dir == -1) {
            get_link(link, U, D*spos+mu);
            adjoint(link,link);
            matmul(seg, seg, link);
        }
    }


    // update loop and position
    matmul(loop, loop, seg);
    *pos = spos;
    }
    
}

int compare(int x, int y, int xsteps, int ysteps) {

    if (xsteps==0 || ysteps==0) {
        printf("You were about to divide by zero. Don't worry, I won't tell anyone.\n");
        exit(1);
    }

    double xratio = (double)((x-0.5)/xsteps);
    double yratio = (double)((y-0.5)/ysteps);

    return xratio > yratio;
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

measurement_list *create_measurement_list(size_t block_size) {
    measurement_list *list = malloc(sizeof(measurement_list));
    if (list != NULL) {
        list->number_of_measurements = 0;
        list->size = block_size;
        list->block_size = block_size;
        list->indices = (uint32_t*)malloc(sizeof(uint32_t*)*block_size);
        list->measurements = (double**)malloc(sizeof(double**)*block_size);
        if (list->indices == NULL || list->measurements == NULL) {
            free(list);
            return NULL;
        }
    }
    return list;
} 

int add_measurement_to_list(measurement_list *list, int index, 
                            double *measurement) {
    size_t num_measurements = list->number_of_measurements;

    if (num_measurements >= list->size) {
        size_t new_size = list->size + list->block_size;
        void *new_indices = realloc(list->indices, sizeof(uint32_t)*new_size);
        void *new_measurements = 
            realloc(list->measurements, sizeof(double*)*new_size);
        if (new_indices == NULL || new_measurements == NULL) {
            return 1;
        } else {
            list->size = new_size;
            list->indices = (uint32_t*)new_indices;
            list->measurements = (double**)new_measurements;
        }
    }

    list->indices[num_measurements] = index;
    /*list->measurements[num_measurements] = mm;*/
    list->measurements[num_measurements] = malloc(2*sizeof(double));
    memcpy(list->measurements[num_measurements],measurement,2*sizeof(double));
    ++list->number_of_measurements;

    return 0;
}

void delete_measurement_list(measurement_list *list) {
    for (size_t i = 0; i < list->number_of_measurements; i++) {
        free(list->measurements[i]);
    }
    free(list->indices);
    free(list->measurements);
    free(list);
}

void source_grid(int *pos, int source, int num_sources) {

    const int little_L = (int)(ceil(pow((num_sources/2),0.25)));
    /*const int little_L = pow((num_sources/2),0.25)+1;*/
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

void compute_observables(double *trace_array, int spos, double *U, int b, int *shape){
    // loop code goes here
        double loop[X];
        memset(loop, 0, sizeof(loop));
        loop[0] = 1;  
        loop[8] = 1;  
        loop[16] = 1;

        diagonal_segment(loop, U, &spos, 1, 2, shape[0], shape[1]);
        segment(loop, U, &spos, 0, b);
        diagonal_segment(loop, U, &spos, 1, 2, -shape[0], -shape[1]);
        diagonal_segment(loop, U, &spos, 1, 2, shape[2], shape[3]);
        segment(loop, U, &spos, 0, -b);
        diagonal_segment(loop, U, &spos, 1, 2, -shape[2], -shape[3]);

        trace(&trace_array[0], loop);
}

double *read_configuration(char *config) {

    double *V;
    V = (double*)calloc(CLENGTH, sizeof(double));

    char *header;
    header = (char*)malloc(600*sizeof(char));
    bzero(header, 600*sizeof(char));
    int header_length = 0;

    if (read_nersc(V, config, header, &header_length)!=0) {
        printf("failed to read configuration!\n");
        exit(1);
    }

    free(header);

    return V;

}

void delete_configuration(double *V) {
    free(V);
}

double *random_gauge_transform(double *V) {
    // for gauge invariance check
    
    double *Om;
    Om = (double*)calloc(CLENGTH, sizeof(double));
    rnd_gauge_tr(Om);
    transform(V, V, Om);
    return Om;
}

double *copy_configuration(double *V) {
    double *U;
    U = (double*)calloc(CLENGTH, sizeof(double));
    memcpy(U, V, sizeof(double)*CLENGTH);
    return U;
}

int *setup_parallization(int *sizes_extent, int num_sources) {
    // Set up partition of sources for parallelization
    int num_threads = omp_get_max_threads();
    int *sizes;
    sizes = (int*)calloc(num_threads+1, sizeof(int));
    sizes[0] = 0; // This is sacrosanct
    if (num_sources < num_threads || num_threads == 1) {
        *sizes_extent = 1;
        sizes[1] = num_sources;
    } else {
        *sizes_extent = num_threads;

        int high = num_sources/(*sizes_extent) + 1;
        int low = num_sources/(*sizes_extent);
        int num_high = num_sources%*sizes_extent;
        for (int i = 0; i < num_high; i++) {
            sizes[i+1] = sizes[i] + high;
        }
        for (int i = num_high; i < *sizes_extent; i++) {
            sizes[i+1] = sizes[i] + low;
        }
    }
    return sizes;
}

void create_info_file(char **inputs) {

    // make an info file
    char info_path[200];
    sprintf(info_path, "%s/info_src%s", inputs[6], inputs[1]);

    struct stat buffer;
    if (stat(info_path, &buffer) != 0) {
        FILE *info_file;
        info_file = fopen(info_path, "w");

        /*fprintf(info_file, header);*/
        fprintf(info_file, "Layout: (uint32_t)(2*double)\n");
        fprintf(info_file, "b=%s\n",inputs[2]);
        fprintf(info_file, "t=%s\n",inputs[3]);

        fclose(info_file);
    }
}

void write_data_to_file(measurement_list *list, char **inputs) {

    int shape[4];
    stoia(shape, inputs[3]);

    // code to write to file
    char data_path[200];
    sprintf(data_path, "%s/b%ssh%d_%d_%d_%dc%s.bin", inputs[6], inputs[2], 
            shape[0], shape[1], shape[2], shape[3], 
            inputs[4]);

    FILE *data_file;
    data_file = fopen(data_path, "ab");
    for (int i = 0; i < list->number_of_measurements; i++) {
        fwrite(&list->indices[i], sizeof(uint32_t), 1, data_file);
        fwrite(list->measurements[i], sizeof(double), 2, data_file);
    }
    fclose(data_file);

}

measurement_list *read_data_from_file(char **inputs) {

    int shape[4];
    stoia(shape, inputs[3]);

    measurement_list *list = create_measurement_list(16);
    // code to read from file
    char data_path[200];
    sprintf(data_path, "%s/b%ssh%d_%d_%d_%dc%s.bin", inputs[6], inputs[2], 
            shape[0], shape[1], shape[2], shape[3], 
            inputs[4]);

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
                add_measurement_to_list(list, index, measurement);
            } else {
                break;
            }
        }

        fclose(data_file);
    }
    return list;
}

uint32_t *get_existing_sources(size_t *number_of_sources, char **inputs) {

    measurement_list *list = read_data_from_file(inputs);
    uint32_t *existant_sources = 
        malloc(sizeof(*existant_sources) * list->number_of_measurements);
    /*if (existant_sources == NULL) {*/
    /*    return 0;*/
    /*}*/

    // don't copy pointer, copy actual memory
    // we want a new location
    for (int i = 0; i < list->number_of_measurements; i++) {
        existant_sources[i] = list->indices[i];
    }
    *number_of_sources = list->number_of_measurements;
    delete_measurement_list(list);

    return existant_sources;
}

int in_source_list(uint32_t source, uint32_t *source_list, size_t size) {

    for (size_t i = 0; i < size; i++) {
        if (source == source_list[i]) {
            return 1;
        }
    }

    return 0;
}


