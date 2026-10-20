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
#include <limits.h>
// #include <gmp.h>
// #include <mpfr.h>
#include <math.h>
#include <time.h>
#include "sizes.h"
#include "su3lib.h"
#include "fileIO.h"

#define LENGTH(X) (sizeof X / sizeof X[0])
#define MAX_SOURCES LVOLUME

enum dim{mu_x, mu_y, mu_z, mu_t};

// Creating the structure
typedef struct {
    uint32_t *indices; // store lattice index for unordered file write
    double *measurements; // array of measurements, flat
    size_t number_of_measurements;
    size_t size;
} measurement_list;

// Struct for input parameters
typedef struct {
    char config[200];
    char config_number[10];
    char data_dir[200];
    char averaged_data_dir[200];
    int32_t num_sources;
    int32_t shape[4];
    int32_t r_direction[4];
    char (*r_list)[16];
    int32_t r_size;
    char (*t_list)[16];
    int32_t t_size;
    char (*b_list)[16];
    int32_t b_size;
    int32_t n; // iterations of (1-H0/2n)
    double w; // wilson parameter

    double *U;
} input_parameters;

// Struct for loop directions
typedef struct {
    const int32_t mu_e;
    const int32_t mu_r;
    const int32_t mu_b;
    int32_t width;
    int64_t bstep;
    int64_t bperiod;
    int64_t estep;
    int64_t eperiod;
} loop_parameters;

//// declaration of functions
// contains all operations including Wilson/mass term, also this has an overall opposite sign from the paper
void hamiltonian(double coef, double r, double w, double *result, double *layer, double *Ulayer, int32_t width, int32_t lpos);
// compute LMN type propagator [add reference]
void tpropLMN(double *prop, double *U, int32_t n, int64_t spos, int32_t mu, 
		 int32_t width, int32_t mu_e, int32_t esteps, int32_t edir, double ratio, double w);
// compute transverse links
void trLinks(double *prop, double *U, int64_t spos, int32_t mu, int32_t width, int32_t bmu, int32_t bsteps, int32_t bdir);

void source_grid(int32_t *pos, int32_t *zpos, int32_t mu_r, int32_t source, int32_t num_sources);
void setup_loop_parameters(loop_parameters *loop);
void compute_observables(double *trace_array, double *prop_num, int32_t spos, const int32_t zspos, 
        const input_parameters *input, const double ratio, const int32_t t, const int32_t b, loop_parameters *loop );
double *get_unit_configuration();
double *read_configuration(char *config);
void delete_configuration(double *V);
double *random_gauge_transform(double *V);
double *copy_configuration(double *V);
void create_info_file(char **inputs); // broken for now 
void create_data_path(char *data_path, input_parameters *input, char *r_str, char *t_str, char *b_str);
void get_existing_sources(size_t *number_of_sources, uint32_t *existant_sources, char *data_path);
int32_t in_source_list(uint32_t source, uint32_t *source_list, size_t size);
void do_run_for_config(input_parameters *input);
////functions for handling inputs
input_parameters *read_input_file_to_struct(char *file_name);
void delete_input_struct(input_parameters *input);
void get_contents(char *out, char* in);
void get_string_array(char (**a)[16], int32_t *size, char *s);
int32_t ssl2sa_alloc(char *s, char **a);
int32_t ssl2ia_fixed(char *s, int32_t *a, size_t size);
void print_input_struct_info(input_parameters *input);

//// functions for handling and saving results of calculation
// TODO: create a clear_measurement_list() function to avoid constant allocator calls
measurement_list *create_measurement_list(size_t block_size);
int64_t add_measurement_to_list(measurement_list *list, int64_t index, 
                            double *measurement);
void delete_measurement_list(measurement_list *list);
void write_data_to_file(measurement_list *list, char *data_path);
measurement_list *read_data_from_file(char *data_path);

//// functions for source averaging
void do_source_average(input_parameters *input); // FIX
measurement_list *read_data_from_file_in_order(char *data_path);
measurement_list *create_full_measurement_list(size_t block_size);
void write_data_to_tsv_file(double *data, char *r_str, char *b_str, char *t_str, input_parameters *input);

int main(int32_t argc, char *argv[]) {
    if (argc != 2) {
        printf("There should only be one argument here. That is the input file.\n");
        exit(1);
    }

    double begin = omp_get_wtime();

    input_parameters *input = read_input_file_to_struct(argv[1]);
    print_input_struct_info(input);

    do_run_for_config(input);
    do_source_average(input);

    // uncomment for gauge invariance check
    printf("Redoing run after random gauge transformation.\n");
    random_gauge_transform(input->U);

    do_run_for_config(input);
    do_source_average(input);

    delete_configuration(input->U);

    // double middle = omp_get_wtime();
    // double time_spent = middle - begin;
    // printf("Measurement elapsed time: %lfs\n", time_spent);

    /* do the source average here */
    // do_source_average(input);

    delete_input_struct(input);

    double end = omp_get_wtime();
    // time_spent = end - middle;
    // printf("Source average elapsed time : %lfs\n", time_spent);
    double time_spent = end - begin;
    printf("Total elapsed time: %lfs\n", time_spent);

    return 0;
}

void do_run_for_config(input_parameters *input) {
    // allocate memory for this here so that so we're not doing it every iteration
    uint32_t *existant_sources;
    existant_sources = (uint32_t*)malloc(MAX_SOURCES * sizeof(uint32_t));

    // setup directions for evolution, transverse, etc.
    loop_parameters loop = {
        .mu_e = mu_z, 
        .mu_r = mu_t, 
        .mu_b = mu_x
    };
    setup_loop_parameters(&loop);

    // there is nothing wrong with a triple for loop when the situation demands it
    // don't be ashamed
    for (int32_t i = 0; i < input->r_size; i++) { // loop through ratio values
    for (int32_t j = 0; j < input->b_size; j++) { // loop through b_perp values
    for (int32_t k = 0; k < input->t_size; k++) { // loop through time values
        char *ratioStr = input->r_list[i];
        char *bStr = input->b_list[j];
        char *tStr = input->t_list[k];
        
        double ratio = atof(input->r_list[i]);
        int32_t b = atoll(input->b_list[j]);
        int32_t t = atoll(input->t_list[k]);

        char data_path[250];
        create_data_path(data_path, input, ratioStr, tStr, bStr); 

        size_t existantSourceNumber; 
        get_existing_sources(&existantSourceNumber, existant_sources, data_path);

        #pragma omp parallel
        {
            // separate list for each thread, everything dumps to the same file
            measurement_list *data_list = 
                create_measurement_list(input->num_sources/omp_get_max_threads() + 1);

            // propagator, this is just an Lx1 array of 3x3 matrices, we don't need anything bigger for this
            double prop_num[loop.width*X];

            #pragma omp for schedule(dynamic, 4)
            for (int32_t s = 0; s < input->num_sources; s++) {
                // map source index to lattice index
                int32_t lattice_s, lattice_zs;
                source_grid(&lattice_s, &lattice_zs, loop.mu_r, s, input->num_sources);

                // check if source has already been measured
                if (in_source_list(lattice_s, existant_sources, 
                                   existantSourceNumber) != 0) {
                    continue;
                }

                double trace_array[NUM_DATAS];

                // here is where the physics happens
                compute_observables(trace_array, prop_num, lattice_s, lattice_zs, input, ratio, t, b, &loop);

                // print result of compute_observables
                // printf("t=%d, b=%d, thread=%d, src=%d, %lf, %lf\n", t, b, parallel_s, s, trace_array[0], trace_array[1]);

                if(add_measurement_to_list(data_list, 
                                           lattice_s, trace_array) != 0) {
                    printf("failed to add measurement to list!\n");
                    printf("this is probably due to memory exhaustion.\n");
                    exit(1);
                }

                if (data_list->number_of_measurements > (size_t)20000) {
                    #pragma omp critical
                    write_data_to_file(data_list, data_path);
                    delete_measurement_list(data_list);
                    data_list = create_measurement_list(16);
                }
            }
            #pragma omp critical

            write_data_to_file(data_list, data_path);
            delete_measurement_list(data_list);
        }
    } // triple
    } // for
    } // loop

    free(existant_sources);
    // free(sizes);
}

void hamiltonian(double coef, double r, double w, double *result, double *layer, double *Ulayer, int32_t width, int32_t lpos) {
    // this is not really a hamiltonian after what I've done to it

    lpos = mod(lpos,width); //mod it just in case
    int32_t blpos = mod(lpos-1, width); //bwd site pos
    int32_t flpos = mod(lpos+1, width); //fwd site pos

    double fsite[X], bsite[X];

    matmul_na_nc(fsite, &layer[X*flpos], &Ulayer[X*lpos]);
    matmul_nc(bsite, &layer[X*blpos], &Ulayer[X*blpos]);
    double *site = &layer[X*lpos]; // this can simply be a pointer since it doesn't multiply links

    // I have decided to just put both derivatives in the same function. 
    for (size_t i = 0; i < (size_t)X/2; i++) {
        size_t j = 2*i;
        //                     fwd-bwd derivative                              laplacian
        result[j]   = coef*( r*0.5*(fsite[j+1] - bsite[j+1]) - fabs(r)*w*0.5*(fsite[j]   + bsite[j]   - 2*site[j]  ));
        result[j+1] = coef*(-r*0.5*(fsite[j]   - bsite[j]  ) - fabs(r)*w*0.5*(fsite[j+1] + bsite[j+1] - 2*site[j+1]));
        // the index swapping in the fwd-bwd derivative accounts for the multiplication by imaginary number i
        // the overall minus sign is correct, laplacian is semi-negative definite
    }
}

void tpropLMN(double *prop, double *U, int32_t n, int64_t spos, int32_t mu, 
        int32_t width, int32_t mu_e, int32_t esteps, int32_t edir, double ratio, 
        double w) {

    // this function is based on a result in arXiv:hep-lat/9205007

    int64_t estep = ipow(L,mu_e); // since T is the last dimension, it won't effect the step size for any direction
    int64_t eperiod = (mu_e == mu_t ? T : L)*estep;

    int64_t step = ipow(L,mu); // mu=0,1,2,3 // this one stays as L, it is the step size
    int64_t period = width*step;

    int64_t rpos = mod(spos/step, width);
    // pull the WL direction back to the edge so
    // we don't have to worry about boundaries
    spos = spos - mod(spos, period) + mod(spos - step*rpos, period);

    double helper[width*X];
    double Ulayer[width*X];

    // fill Ulayer for first use
    int32_t l;
    for (l = 0; l < width; l++) {
        memcpy(&Ulayer[X*l], &U[X*(D*(spos+step*l)+mu)], sizeof(double)*X);
    }

    double coef = (double)1.0/(2*n);
    double coef2 = 0.5*coef;

    // also initiliaze arrays for iterative updates
    double update[X];

    // iterators
    int32_t i, j, k, tn;
    for (i = 0; i != esteps; i++) {

        // (1-H0/2n)^n |_tprev
        tn = n;
        while (tn--) {
            for (j = 0; j < width; j++) { // scan across (x,y,z,t) direction
                hamiltonian(coef, ratio, w, update, prop, Ulayer, width, j); // fwd-bwd derivative
                for (k = 0; k < X; k++) {
                    helper[X*j+k] = update[k];
                }
            }

            // apply the hamiltonian again to the result of the previous application to get H^2
            for (j = 0; j < width; j++) { // scan across (x,y,z,t) direction
                hamiltonian(coef2, ratio, w, update, helper, Ulayer, width, j); // fwd-bwd derivative
                for (k = 0; k < X; k++) {
                    prop[X*j+k] += update[k] - helper[X*j+k];
                }
            }
        }

        // multiply time links if edir = +1
        // U^+(tprev)
        if (edir == 1) {
            for(j = 0; j < width; j++) {
                int64_t pos = X*(D*(spos+step*j)+mu_e);
                matmul_nc(&helper[X*j],&prop[X*j],&U[pos]);
            }
        }

        // prepare spos for next evol slice
    	  spos = spos - mod(spos, eperiod) + mod(spos+estep*edir, eperiod);

        // multiply time links if edir = -1
        // U^+(tprev)
        if (edir == -1) {
            for(j = 0; j < width; j++) {
                int64_t pos = X*(D*(spos+step*j)+mu_e);
                matmul_na_nc(&helper[X*j],&prop[X*j],&U[pos]);
            }
        }

        // update Ulayer for next loop
        for (l = 0; l < width; l++) {
            memcpy(&Ulayer[X*l], &U[X*(D*(spos+step*l)+mu)], sizeof(double)*X);
        }

        // (1-H0/2n)^n |_t
        tn = n;
        while (tn--) {
            // TODO: this section and the similar one before could be abstracted, and probably made slightly faster
            for (j = 0; j < width; j++) { // scan across x direction
                hamiltonian(coef, ratio, w, update, helper, Ulayer, width, j); // Hamiltonian
                for (k = 0; k < X; k++) {
                    prop[X*j+k] = update[k];
                }
            }

            // apply the hamiltonian again to the result of the previous application to get H^2
            for (j = 0; j < width; j++) { // scan across (x,y,z,t) direction
                hamiltonian(coef2, ratio, w, update, prop, Ulayer, width, j); // fwd-bwd derivative
                for (k = 0; k < X; k++) {
                    helper[X*j+k] += update[k] - prop[X*j+k];
                }
            }
        }
        memcpy(prop, helper, sizeof(double)*width*X);
    }
}

void trLinks(double *prop, double *U, int64_t spos, int32_t mu, int32_t width, 
        int32_t bmu, int32_t bsteps, int32_t bdir) {

    int64_t step = ipow(L,mu);
    int64_t period = width*step;
    int64_t rpos = mod(spos/step, width);

    int64_t bstep = ipow(L,bmu);
    int64_t bperiod = L*bstep; // this will not work if we want the time direction
                          // to be our bdir but why would we want that?

    // pull the WL direction back to the edge so
    // we don't have to worry about boundaries
    spos = spos - mod(spos, period) + mod(spos - step*rpos, period);

    double layer[width*X]; // chocolate cake

    // iterators
    int32_t i, j;
    for (i = 0; i != bsteps; i++) { 

        memcpy(layer, prop, sizeof(double)*width*X);

        if (bdir == 1) {
            for(j = 0; j < width; j++) {
                int64_t pos = X*(D*(spos+step*j)+bmu);
                matmul(&layer[X*j],&layer[X*j],&U[pos]);
            }
        }

        // prepare spos for next step
        spos = spos - mod(spos, bperiod) + mod(spos + bstep*bdir, bperiod);

        if (bdir == -1) {
            for(j = 0; j < width; j++) {
                int64_t pos = X*(D*(spos+step*j)+bmu);
                matmul_na(&layer[X*j],&layer[X*j],&U[pos]);
            }
        }

        memcpy(prop, layer, sizeof(double)*width*X);
    }
}

measurement_list *create_measurement_list(size_t size) {
    measurement_list *list = malloc(sizeof(measurement_list));
    if (list != NULL) {
        list->number_of_measurements = 0;
        list->size = size;
        list->indices = malloc(size*sizeof(uint32_t));
        list->measurements = malloc(NUM_DATAS*size*sizeof(double));
        if (list->indices == NULL || list->measurements == NULL) {
            // QSTN: do we need to free list->indices and list->measurements if they are already NULL?
            free(list->indices);
            free(list->measurements);
            free(list);
            return NULL;
        }
    }
    return list;
} 

int64_t add_measurement_to_list(measurement_list *list, int64_t index, 
                            double *measurement) {
    size_t num_measurements = list->number_of_measurements;
    if (num_measurements >= list->size) {
        size_t new_size = 2*list->size;
        void *new_indices = realloc(list->indices, new_size*sizeof(uint32_t));
        if (new_indices == NULL) return 1;
        void *new_measurements = 
            realloc(list->measurements, NUM_DATAS*new_size*sizeof(double));
        if (new_measurements == NULL) return 1;
        list->size = new_size;
        list->indices = (uint32_t*)new_indices;
        list->measurements = (double*)new_measurements;
    }
    list->indices[num_measurements] = index;
    memcpy(&list->measurements[NUM_DATAS*num_measurements], measurement, NUM_DATAS*sizeof(double));
    ++list->number_of_measurements;
    return 0;
}

void delete_measurement_list(measurement_list *list) {
    free(list->indices);
    free(list->measurements);
    free(list);
}

void source_grid(int32_t *pos, int32_t *zpos, int32_t mu_r, int32_t source, int32_t num_sources) {
    // TODO: rework this
    int32_t src[] = {0, 0, 0, 0};

    if (num_sources != 1) {
        const int32_t little_L = (int32_t)(ceil(pow((num_sources/2),0.25)));
        /*const int32_t little_L = pow((num_sources/2),0.25)+1;*/
        const int32_t little_T = 2*little_L;
        const int32_t fact_L = L/little_L;
        const int32_t fact_T = T/little_T;

        src[0] = fact_L*(source%little_L);
        src[1] = fact_L*((source/little_L)%little_L);
        src[2] = fact_L*(((source/little_L)/little_L)%little_L);
        src[3] = fact_T*((((source/little_L)/little_L)/little_L)%little_T);
    }

    if (zpos != NULL) *zpos = src[mu_r]; // for propagator which is one dimensional array in our case
    *pos = src[3]*LCUBE + src[2]*LSQUARE + src[1]*L + src[0];
}

void setup_loop_parameters(loop_parameters *loop) {
    loop->width = (loop->mu_r == mu_t ? T : L);
    loop->bstep = ipow(L,loop->mu_b);
    loop->bperiod = L*loop->bstep; // period for b direction, not generalized to loop->mu_b=mu_t
    loop->estep = ipow(L,loop->mu_e);
    loop->eperiod = (loop->mu_e == mu_t ? T : L)*loop->estep;
}

void compute_observables(double *trace_array, double *prop_num, int32_t spos, const int32_t zspos, 
        const input_parameters *input, const double ratio, const int32_t t, const int32_t b, loop_parameters *loop ) {

    double *U = input->U;

    // int64_t srcp = X*spos;
    int32_t zsrcp = X*zspos;
    double source[] = {1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0};

    int64_t bstep = loop->bstep;
    int64_t bperiod = loop->bperiod; // period for b direction, not generalized to loop->mu_b=mu_t
    int64_t estep = loop->estep;
    int64_t eperiod = loop->eperiod;

    int32_t bdir; // edir is direction of evolution, usually t or z direction

    memset(prop_num, 0, sizeof(double)*loop->width*X);
    memcpy(&prop_num[zsrcp],source,sizeof(double)*X);

    tpropLMN(prop_num, U, input->n, spos, loop->mu_r, loop->width, loop->mu_e, 
            t, input->shape[0], ratio*input->r_direction[0], input->w);
    spos = spos - mod(spos, eperiod) + mod(spos+estep*t*input->shape[0], eperiod);

    trLinks(prop_num, U, spos, loop->mu_r, loop->width, loop->mu_b, b, bdir=1);
    spos = spos - mod(spos, bperiod) + mod(spos+bstep*b*bdir, bperiod);

    tpropLMN(prop_num, U, input->n, spos, loop->mu_r, loop->width, loop->mu_e, 
            t, input->shape[1], ratio*input->r_direction[1], input->w);
    spos = spos - mod(spos, eperiod) + mod(spos+estep*t*input->shape[1], eperiod);

    tpropLMN(prop_num, U, input->n, spos, loop->mu_r, loop->width, loop->mu_e, 
            t, input->shape[2], ratio*input->r_direction[2], input->w);
    spos = spos - mod(spos, eperiod) + mod(spos+estep*t*input->shape[2], eperiod);

    trLinks(prop_num, U, spos, loop->mu_r, loop->width, loop->mu_b, b, bdir=-1);
    spos = spos - mod(spos, bperiod) + mod(spos+bstep*b*bdir, bperiod);

    tpropLMN(prop_num, U, input->n, spos, loop->mu_r, loop->width, loop->mu_e, 
            t, input->shape[3], ratio*input->r_direction[3], input->w);

    trace(&trace_array[0], &prop_num[zsrcp]);

    // // print out the sum trace along the z direction
    // double tra[2];
    // double tsum_real = 0.0;
    // double tsum_imag = 0.0;
    // for (int32_t i = 0; i < width; ++i) {
    //     trace(&tra[0], &prop_num[X*i]);
    //     tsum_real += tra[0]; 
    //     tsum_imag += tra[1]; 
    // }
    // printf("sum trace: %lf, %lf\n", tsum_real, tsum_imag);
}

double *get_unit_configuration() {
    double *V;
    V = (double*)malloc(CLENGTH*sizeof(double));

    unit_config(V);

    return V;
}

double *read_configuration(char *config) {
    double *V;
    // CHECK
    V = (double*)malloc(CLENGTH*sizeof(double));

    char header[1000];
    bzero(header, 1000*sizeof(char));
    int32_t header_length = 0;

    if (read_nersc(V, config, header, &header_length)!=0) {
        printf("failed to read configuration!\n");
        exit(1);
    }

    return V;
}

void delete_configuration(double *V) {
    free(V);
}

double *random_gauge_transform(double *V) {
    // for gauge invariance check
    
    double *Om;
    Om = (double*)malloc(OLENGTH*sizeof(double));
    rnd_gauge_tr(Om);
    transform(V, V, Om);
    return Om;
}

double *copy_configuration(double *V) {
    double *U;
    U = (double*)malloc(CLENGTH*sizeof(double));
    memcpy(U, V, sizeof(double)*CLENGTH);
    return U;
}

void create_info_file(char **inputs) {
    // TODO: change how we do this. It's basically useless right now
    // make an info file
    char info_path[200];
    sprintf(info_path, "%s/info_src%s", inputs[8], inputs[3]);

    struct stat buffer;
    if (stat(info_path, &buffer) != 0) {
        FILE *info_file;
        info_file = fopen(info_path, "a");

        /*fprintf(info_file, header);*/
        fprintf(info_file, "Layout: (uint32_t)(NUM_DATAS*double)\n");
        fprintf(info_file, "cols:\tnum.real\tnum.imag\n");
        fprintf(info_file, "n=%s\n",inputs[1]);
        fprintf(info_file, "r=%s\n",inputs[2]);
        fprintf(info_file, "b=%s\n",inputs[4]);
        fprintf(info_file, "t=%s\n",inputs[5]);
        fprintf(info_file, "c=%s\n",inputs[6]);
        fprintf(info_file, "\n");
        /*fprintf(info_file, "index = b + ncols * t\n");*/

        fclose(info_file);
    }
}

void create_data_path(char *data_path, input_parameters *input, char *r_str, char *t_str, char *b_str) {
    char data_sub_dir[215];
    sprintf(data_sub_dir, "%s/r%sb%st%s", input->data_dir, 
            r_str, b_str, t_str);
    mkdir(data_sub_dir, 0755);

    sprintf(data_path, "%s/c%s.bin", data_sub_dir, input->config_number);
}

void write_data_to_file(measurement_list *list, char *data_path) {
    // code to write to file

    FILE *data_file;
    data_file = fopen(data_path, "ab");
    for (size_t i = 0; i < list->number_of_measurements; i++) {
        fwrite(&list->indices[i], sizeof(uint32_t), 1, data_file);
        fwrite(&list->measurements[NUM_DATAS*i], sizeof(double), NUM_DATAS, data_file);
    }
    fclose(data_file);
}

measurement_list *read_data_from_file(char *data_path) {
    measurement_list *list = create_measurement_list(16);
    printf("Looking for file: %s\n", data_path);
    FILE *data_file = fopen(data_path, "rb");
    if (data_file != NULL) {
        size_t num_read;
        for (int64_t i=0; i<LVOLUME; i++) {
            num_read = 0;
            uint32_t index;
            double measurement[NUM_DATAS];
            num_read += fread(&index, sizeof(uint32_t), 1, data_file);
            num_read += fread(measurement, sizeof(double), NUM_DATAS, data_file);
            if (num_read == NUM_DATAS+1) {
                add_measurement_to_list(list, index, measurement);
            } else {
                break;
            }
        }

        fclose(data_file);
    }
    else {
        printf("Data file not found. A new data file will be created.\n");
    }
    return list;
}

void get_existing_sources(size_t *number_of_sources, uint32_t *existant_sources, char *data_path) {
    measurement_list *list = read_data_from_file(data_path);

    // don't copy pointer, copy actual memory
    // we want a new location
    for (size_t i = 0; i < list->number_of_measurements; i++) {
        existant_sources[i] = list->indices[i];
    }
    *number_of_sources = list->number_of_measurements;
    delete_measurement_list(list);
}

int32_t in_source_list(uint32_t source, uint32_t *source_list, size_t size) {
    for (size_t i = 0; i < size; i++) {
        if (source == source_list[i]) {
            return 1;
        }
    }

    return 0;
}

//// functions for handling program input
input_parameters *read_input_file_to_struct(char *file_name) {
    // TODO: update this so it better handles bad input files and extra spaces
    const size_t NUM_PARAMS = 12;

    char *PARAM_LABELS[] = {
        "CONFIG_FILE",
        "CONFIG_NUM",
        "DATA_DIR",
        "AVERAGED_DATA_DIR",
        "NUM_SRC",
        "SHAPE",
        "R_SHAPE",
        "T_LIST",
        "B_LIST",
        "R_LIST",
        "N_PARAM",
        "W_PARAM"
    };
    size_t PARAM_SIZES[] = {11, 10, 8, 17, 7, 5, 7, 6, 6, 6, 7, 7 };

    char PARAM_STRINGS[NUM_PARAMS][200];

    input_parameters *input;
    input = malloc(sizeof(input_parameters));

    char line[200];

    FILE *infile = fopen(file_name, "r");
    if (infile == NULL) {
        printf("Cannot find input file. Exiting.\n");
        exit(1);
    }

    while( fgets(line, 200, infile) != NULL) {
        int res;
        for (size_t i = 0; i < NUM_PARAMS; ++i) {
            res = memcmp(line, PARAM_LABELS[i], PARAM_SIZES[i]);
            if (res == 0) {
                char *str = &line[PARAM_SIZES[i]+3];
                get_contents(PARAM_STRINGS[i], str);
                continue;
            }
        }
    }

    size_t index = 0;

    memcpy(input->config, PARAM_STRINGS[index], sizeof(input->config));
    memcpy(input->config_number, PARAM_STRINGS[++index], sizeof(input->config_number));
    memcpy(input->data_dir, PARAM_STRINGS[++index], sizeof(input->data_dir));
    memcpy(input->averaged_data_dir, PARAM_STRINGS[++index], sizeof(input->averaged_data_dir));

    input->num_sources = atoi(PARAM_STRINGS[++index]);

    ssl2ia_fixed(PARAM_STRINGS[++index], input->shape, 4);
    ssl2ia_fixed(PARAM_STRINGS[++index], input->r_direction, 4);

    // input->t_size = ssl2sa_alloc(PARAM_STRINGS[++index], NULL);
    // input->t_list = malloc(input->t_size*sizeof(char*));
    // ssl2sa_alloc(PARAM_STRINGS[index], input->t_list);
    get_string_array(&input->t_list, &input->t_size, PARAM_STRINGS[++index]);

    // input->b_size = ssl2sa_alloc(PARAM_STRINGS[++index], NULL);
    // input->b_list = malloc(input->b_size*sizeof(char*));
    // ssl2sa_alloc(PARAM_STRINGS[index], input->b_list);
    get_string_array(&input->b_list, &input->b_size, PARAM_STRINGS[++index]);

    // input->r_size = ssl2sa_alloc(PARAM_STRINGS[++index], NULL);
    // input->r_list = malloc(input->r_size*sizeof(char*));
    // ssl2sa_alloc(PARAM_STRINGS[index], input->r_list);
    get_string_array(&input->r_list, &input->r_size, PARAM_STRINGS[++index]);

    input->n = atoi(PARAM_STRINGS[++index]);
    input->w = atof(PARAM_STRINGS[++index]);

    fclose(infile);

    // TODO: maybe have this separate
    input->U = read_configuration(input->config);

    mkdir(input->data_dir, 0755);
    mkdir(input->averaged_data_dir, 0755);

    return input;
}

void delete_input_struct(input_parameters *input) {
    // for (int32_t i = 0; i < input->r_size; i++) {
    //     free(input->r_list[i]);
    // }
    // for (int32_t i = 0; i < input->b_size; i++) {
    //     free(input->b_list[i]);
    // }
    // for (int32_t i = 0; i < input->t_size; i++) {
    //     free(input->t_list[i]);
    // }
    free(input->r_list);
    free(input->b_list);
    free(input->t_list);
    free(input);
}

void get_contents(char *out, char* in) {
    size_t i = 0;
    size_t j = 0;
    while (in[i] != '\n') {
        if (in[i] != '"') {
            out[j] = in[i];
            j++;
        }
        i++;
    }
    out[j] = '\0';
}

void get_string_array(char (**a)[16], int32_t *size, char *s) {
    int32_t len = 8;
    *a = (char(*)[16])malloc(len * sizeof **a);
    *size = 0;

    int32_t i = 0;
    int32_t ri = 0;
    while (1) {
        if ( s[i] == ' ' || s[i] == '\0') {
            int32_t entry_len = i-ri;
            if (entry_len >= 16) {
                printf("One or more inputs to T_LIST, B_LIST or R_LIST larger than 15 digits. This is pointless.\n");
                exit(1);
            }
            strncpy((*a)[*size], s+ri, entry_len*sizeof(char));
            (*a)[*size][i-ri] = '\0';
            ri = i + 1;
            *size += 1;
            if (*size == len) {
                void *new_array = (char(*)[16])realloc(*a, 2*len * sizeof **a);
                if (new_array == NULL) {
                    printf("realloc returned NULL array in get_string_array()\n");
                    exit(1);
                }
                *a = new_array;
                len *= 2;
            }
        }
        if (s[i] == '\0') break;
        i++;
    }
}

int32_t ssl2sa_alloc(char *s, char **a) {
    int32_t c = 0;
    if (a == NULL) {
        int32_t i = 0;
        while (1) {
            if ( s[i] == ' ' || s[i] == '\0') {
                c++;
            }
            if (s[i] == '\0') break;
            i++;
        }
    }
    else {
        int32_t i = 0;
        int32_t ri = 0;
        while (1) {
            if ( s[i] == ' ' || s[i] == '\0') {
                *(a+c) = malloc((i-ri+1)*sizeof(char));
                strncpy(*(a+c), s+ri, i-ri);
                ri = i + 1;
                c++;
            }
            if (s[i] == '\0') break;
            i++;
        }
    }
    return c;
}

int32_t ssl2ia_fixed(char *s, int32_t *a, size_t size) {
    size_t c = 0;
    size_t i = 0;
    size_t ri = 0;
    while (1) {
        if ( s[i] == ' ' || s[i] == '\0') {
            a[c] = atoi(&s[ri]); 
            ri = i+1;
            c++;
        }
        if (s[i] == '\0' || c == size) {
            break;
        }
        i++;
    }
    return c;
}

void print_input_struct_info(input_parameters *input) {
    printf("Doing run with following inputs:\n");
    printf("Configuration file: %s\n", input->config);
    printf("Configuration number: %s\n", input->config_number);
    printf("Data write directory: %s\n", input->data_dir);
    printf("Source averaged data write directory: %s\n", input->averaged_data_dir);
    printf("Number of sources: %d\n", input->num_sources);
    printf("Auxiliary propagator t directions: %d, %d, %d, %d\n", 
            input->shape[0], input->shape[1], input->shape[2], input->shape[3]);
    printf("Auxiliary propagator r directions: %d, %d, %d, %d\n", 
            input->r_direction[0], input->r_direction[1], input->r_direction[2], input->r_direction[3]);
    printf("Ratio inputs: ");
    for (int32_t i = 0; i < input->r_size; i++) {
        printf("%s", input->r_list[i]);
        if (i < input->r_size-1) {
            printf(", ");
        }
        else {
            printf("\n");
        }
    }
    printf("Time inputs: ");
    for (int32_t i = 0; i < input->t_size; i++) {
        printf("%s", input->t_list[i]);
        if (i < input->t_size-1) {
            printf(", ");
        }
        else {
            printf("\n");
        }
    }
    printf("bT inputs: ");
    for (int32_t i = 0; i < input->b_size; i++) {
        printf("%s", input->b_list[i]);
        if (i < input->b_size-1) {
            printf(", ");
        }
        else {
            printf("\n");
        }
    }
    printf("Hamiltonian iterations, n: %d\n", input->n);
    printf("Wilson parameter, w: %lf\n", input->w);
}

//// functions for source averaging
void do_source_average(input_parameters *input) {
    for (int32_t i = 0; i < input->r_size; i++) { // loop through ratio values
    for (int32_t j = 0; j < input->b_size; j++) { // loop through b_perp values
    for (int32_t k = 0; k < input->t_size; k++) { // loop through time values
                                                  //
        char *r_str = input->r_list[i];
        char *b_str = input->b_list[j];
        char *t_str = input->t_list[k];

        char data_path[250];
        create_data_path(data_path, input, r_str, t_str, b_str); 

        double averaged_data[NUM_DATAS];
        memset(averaged_data, 0, NUM_DATAS*sizeof(double));

        measurement_list *ordered_data = read_data_from_file_in_order(data_path);

        for (int64_t s = 0; s < input->num_sources; ++s) {
            // TODO: lattice_zs is not needed here, maybe rework this
            int32_t lattice_s;
            source_grid(&lattice_s, NULL, 0, s, input->num_sources);
            if (ordered_data->indices[lattice_s] == UINT32_MAX) {
                printf("missing data!\n");
                exit(1);
            }
            for (size_t i = 0; i < NUM_DATAS; i++) {
                averaged_data[i] += ordered_data->measurements[NUM_DATAS*lattice_s+i];
            }
        }

        delete_measurement_list(ordered_data);

        for (size_t i = 0; i<NUM_DATAS; i++) {
            averaged_data[i] /= (double)input->num_sources;
            // printf("%lf ", averaged_data[i]);
        }
        // printf("\n");

        write_data_to_tsv_file(averaged_data, r_str, b_str, t_str, input);
    } // triple
    } // for
    } // loop
}

measurement_list *create_full_measurement_list(size_t size) {
    measurement_list *list = malloc(sizeof(measurement_list));
    if (list != NULL) {
        list->number_of_measurements = size;
        list->size = size;
        list->indices = malloc(size*sizeof(uint32_t));
        list->measurements = malloc(NUM_DATAS*size*sizeof(double));
        if (list->indices == NULL || list->measurements == NULL) {
            free(list->indices);
            free(list->measurements);
            free(list);
            return NULL;
        }
    }
    return list;
} 

measurement_list *read_data_from_file_in_order(char *data_path) {
    measurement_list *list = create_full_measurement_list(LVOLUME);
    memset(list->indices, UINT32_MAX, LVOLUME*sizeof(uint32_t));
    FILE *data_file = fopen(data_path, "rb");
    if (data_file != NULL) {
        size_t num_read;
        for (int64_t i=0; i<LVOLUME; i++) {
            num_read = 0;
            uint32_t index;
            double measurement[NUM_DATAS];
            num_read += fread(&index, sizeof(uint32_t), 1, data_file);
            num_read += fread(measurement, sizeof(double), NUM_DATAS, data_file);
            if (num_read == NUM_DATAS+1) {
                list->indices[index] = index;
                memcpy(&list->measurements[NUM_DATAS*index],
                       measurement,2*sizeof(double));
            } else {
                break;
            }
        }
        fclose(data_file);
    }
    return list;
}

void write_data_to_tsv_file(double *data, char *r_str, char *b_str, char *t_str, input_parameters *input) {
    // code to write to file
    char averaged_data_path[215];
    sprintf(averaged_data_path, "%s/c%s.tsv", 
            input->averaged_data_dir, input->config_number);

    struct stat buffer;
    int64_t result = stat(averaged_data_path, &buffer);

    FILE *tsv_file;
    tsv_file = fopen(averaged_data_path, "a");

    if (result != 0) {
        fprintf(tsv_file, "# r\tb\tt\treal\timag\t\n");
    }

    fprintf(tsv_file, "%s\t%s\t%s\t%lf\t%lf\t\n", 
            r_str, b_str, t_str, data[0], data[1]);

    fclose(tsv_file);
}

