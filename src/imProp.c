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
    double **measurements; // array of measurements, probably 6 doubles
    size_t number_of_measurements;
    size_t size;
    size_t block_size;
} measurement_list;

//// declaration of functions
// compute forward-backward derivative
void fbdelta(double *result, double *layer, double *Ulayer, int lpos);
// compute LMN type propagator [add reference]
void tpropLMNnew(double *prop, double *U, int n, int spos, int mu, int tsteps, 
                 int tdir, double ratio);
// compute transverse links
void trLinksnew(double *prop, double *U, int spos, int mu, int bmu, 
                int bsteps, int bdir);
// create a struct containing dynamic list of measurements
measurement_list *create_measurement_list(size_t block_size);
// append a new measurement with lattice index to the list
int add_measurement_to_list(measurement_list *list, int index, 
                            double *measurement);
// delete the list
void delete_measurement_list(measurement_list *list);
void source_grid(int *pos, int source, int num_sources);
void compute_observables(double *trace_array, int spos, double *U, int n, 
                         double ratio, int t, int b,
                         double *prop_num, double *prop_den);
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
int ssl2sa(char *s, char **a); // space-separated list to string array

// void build_string(char *out, char *in[], int length);

int main(int argc, char *argv[]) {

    if (argc < 9 || argc > 9) {
        printf("Input should look like: \n");
        printf("%s <n> <ratio list> num_sources <b list> <t list> config_number"
               " config_path data_dir\n", argv[0]);
        exit(1);
    }

    // clock_t begin = clock();
    double begin = omp_get_wtime();

    int n = atoi(argv[1]);

    int rSize = ssl2sa(argv[2], NULL);
    char **ratioList = malloc(rSize*sizeof(char*));
    ssl2sa(argv[2], ratioList);

    int nSources = atoi(argv[3]);

    int bSize = ssl2sa(argv[4], NULL);
    char **bList = malloc(bSize*sizeof(char*));
    ssl2sa(argv[4], bList);

    int tSize = ssl2sa(argv[5], NULL);
    char **tList = malloc(tSize*sizeof(char*));
    ssl2sa(argv[5], tList);

    char *config_number = argv[6];
    char *config = argv[7];
    char *data_dir = argv[8];

    double *V = read_configuration(config);


    // change triple for-loop to single loop later
    for (int i = 0; i < rSize; i++) {
    for (int j = 0; j < bSize; j++) {
    for (int k = 0; k < tSize; k++) {

        char *ratioStr = ratioList[i];
        char *bStr = bList[j];
        char *tStr = tList[k];
        char *inargs[] = { argv[0], argv[1], ratioStr, argv[3], bStr, tStr, 
                            argv[6], argv[7], argv[8] };
        
        double ratio = atof(ratioStr);
        int b = atoll(bStr);
        int t = atoll(tStr);
        printf("Calculating for %lf, %d, %d\n", ratio, b, t);

        printf("here1\n");
        create_info_file(inargs); // if it doesn't already exist
        printf("here2\n");

        size_t existantSourceNumber; 
        uint32_t *existantSources = 
            get_existing_sources(&existantSourceNumber, inargs);

        // setup parallelization // move this outside loops, global variable
        int sizes_extent;
        int *sizes = setup_parallization(&sizes_extent, nSources);

        #pragma omp parallel for
        for (int parallel_s = 0; parallel_s < sizes_extent; parallel_s++) {

            // can't remember why i was making copies, but it doesn't seem to help
            // and uses a lot more memory
            // double *U = copy_configuration(V);
            measurement_list *data_list = 
                create_measurement_list(nSources/sizes_extent+1);

            // propagators
            double *prop_num = (double*)calloc(OLENGTH, sizeof(double));
            double *prop_den = (double*)calloc(OLENGTH, sizeof(double));

            for (int s = sizes[parallel_s]; s < sizes[parallel_s+1]; s++) {

                // map source index to lattice index
                int lattice_s;
                source_grid(&lattice_s, s, nSources);
                // check if source has already been measured
                if (in_source_list(lattice_s, existantSources, 
                                   existantSourceNumber) != 0) {
                    continue;
                }

                double trace_array[NUM_DATAS];

                // here is where the physics happens
                compute_observables(trace_array, lattice_s, V, n, ratio, t, b, 
                                    prop_num, prop_den);

                if(add_measurement_to_list(data_list, 
                                           lattice_s, trace_array) != 0) {
                    printf("failed to add measurement to list!\n");
                    printf("this is probably due to memory exhaustion.\n");
                    exit(1);
                }

                if (data_list->number_of_measurements > (size_t)20000) {
                    #pragma omp critical
                    write_data_to_file(data_list, inargs);
                    delete_measurement_list(data_list);
                    data_list = create_measurement_list(16);
                }
            }

            free(prop_num);
            free(prop_den);

            #pragma omp critical
            write_data_to_file(data_list, inargs);
            delete_measurement_list(data_list);
            // delete_configuration(U);
        }
        #pragma omp barrier

        free(sizes);
        free(existantSources);

    } // triple
    } // for
    } // loop

    delete_configuration(V);

    // Clean this up
    for (int i = 0; i < rSize; i++) {
        free(ratioList[i]);
    }
    for (int j = 0; j < bSize; j++) {
        free(bList[j]);
    }
    for (int k = 0; k < tSize; k++) {
        free(tList[k]);
    }
    free(ratioList);
    free(bList);
    free(tList);

    double end = omp_get_wtime();
    double time_spent = end - begin;
    printf("Elapsed time: %lfs\n", time_spent);

    return 0;
}

void fbdelta(double *result, double *layer, double *Ulayer, int lpos) {

    lpos = mod(lpos,L); //mod it just in case

    int blpos = mod(lpos-1, L); //bwd site pos
    int flpos = mod(lpos+1, L); //fwd site pos

    double fsite[X], bsite[X];
    double flink[X], blink[X];

    memcpy(fsite,&layer[X*flpos],sizeof(fsite));
    memcpy(bsite,&layer[X*blpos],sizeof(bsite));
    memcpy(flink,&Ulayer[X*lpos],sizeof(flink));
    adjoint(flink, flink);
    memcpy(blink,&Ulayer[X*blpos],sizeof(blink));

    matmul(fsite, fsite, flink);
    matmul(bsite, bsite, blink);

    for (int i = 0; i < X; i++) {
        result[i] = fsite[i] - bsite[i];
    }

    // coolest possible way to multiply by a scalar
    double coef = 0.5;
    double coefm[X];
    memset(coefm, 0, sizeof(coefm));
    coefm[0] = coef;
    coefm[8] = coef;
    coefm[16] = coef;

    matmul(result, coefm, result);
}

void tpropLMNnew(double *prop, double *U, int n, int spos, int mu, int tsteps, 
                 int tdir, double ratio) {

    // this function is based on a result in arXiv:hep-lat/9205007
    //
    // need to add a description of variables and constants

    int step = ipow(L,mu); // mu=0,1,2,3
    int period = L*step;
    int rpos = mod(spos/step, L);
    // pull the WL direction back to the edge so
    // we don't have to worry about boundaries
    spos = spos - mod(spos, period) + mod(spos - step*rpos, period);

    // allocate memory for layer and helper array
    // double *layer, *helper, *Ulayer;
    // layer = (double*)cal  loc(L*X, sizeof(double));
    // helper = (double*)cal  loc(L*X, sizeof(double));
    // Ulayer = (double*)cal  loc(L*X, sizeof(double));

    double layer[L*X];
    double helper[L*X];
    double Ulayer[L*X];

    // fill Ulayer for first use
    int l;
    for (l = 0; l < L; l++) {
        memcpy(&Ulayer[X*l], &U[X*(D*(spos+step*l)+mu)], sizeof(double)*X);
    }


    // build coefficient matrix, factor of 'i' accounted for 
    // since odd indexes represent the imaginary part
    double coef = (double)tdir * ratio/(2*n);
    double coefm[X];
    memset(coefm, 0, sizeof(coefm));
    coefm[1] = coef;  // 1, 9, 17 correspond to the imaginary
    coefm[9] = coef;  // part on the diagonal
    coefm[17] = coef;

    // also initiliaze arrays for iterative updates
    double fbder[X], update[X];
    double Ut[X];

    // iterators
    int i, j, k, tn;
    for (i = 0; i != tsteps; i++) {

        for (l = 0; l < L; l++) {
            // step variable determines which space direction we are using
            // for our Wilson line direction
            memcpy(&layer[X*l], &prop[X*(spos+step*l)], sizeof(double)*X);
        }
        memcpy(helper, layer, sizeof(double)*L*X);

        // (1-H0/2)^n |_tprev
        tn = n;
        while (tn--) {
            for (j = 0; j < L; j++) { // scan across (x,y,z) direction
                fbdelta(fbder, layer, Ulayer, j); // fwd-bwd derivative
                matmul(update,coefm,fbder);
                for (k = 0; k < X; k++) {
                    helper[X*j+k] += update[k];
                }
            }
            memcpy(layer, helper, sizeof(double)*L*X);
        }

        // multiply time links if tdir = +1
        // U^+(tprev)
        if (tdir == 1) {
            for(j = 0; j < L; j++) {
                // X*D*spos+X*3+X*D*j
                get_link(Ut,U,D*(spos+step*j)+mu_t);
                matmul(&layer[X*j],&layer[X*j],Ut);
            }
        }

        // prepare spos for next time slice
        spos = mod(spos + tdir*LCUBE, LVOLUME);

        // multiply time links if tdir = -1
        // U^+(tprev)
        if (tdir == -1) {
            for(j = 0; j < L; j++) {
                get_link(Ut,U,D*(spos+step*j)+mu_t);
                adjoint(Ut,Ut);
                matmul(&layer[X*j],&layer[X*j],Ut);
            }
        }

        // update helper array for next loop
        memcpy(helper, layer, sizeof(double)*L*X);

        // update Ulayer for next loop
        for (l = 0; l < L; l++) {
            memcpy(&Ulayer[X*l], &U[X*(D*(spos+step*l)+mu)], sizeof(double)*X);
        }

        // (1-H0/2)^n |_t
        tn = n;
        while (tn--) {
            for (j = 0; j < L; j++) { // scan across x direction
                fbdelta(fbder, layer, Ulayer, j); // fwd-bwd derivative
                matmul(update,coefm,fbder);
                for (k = 0; k < X; k++) {
                    helper[X*j+k] += update[k];
                }
            }
            memcpy(layer, helper, sizeof(double)*L*X);
        }

        // update propagator
        for (l = 0; l < L; l++) {
            // step variable determines which space direction we are using
            // for our Wilson line direction
            memcpy(&prop[X*(spos+step*l)], &layer[X*l], sizeof(double)*X);
        }
    }

    // free(layer);
    // free(helper);
    // free(Ulayer);
}


void trLinksnew(double *prop, double *U, int spos, int mu, int bmu, 
                int bsteps, int bdir) {

    int step = ipow(L,mu);
    int period = L*step;
    int rpos = mod(spos/step, L);
    int bstep = ipow(L,bmu);
    int bperiod = L*bstep; // this will not work if we want the time direction
                          // to be our bdir but why would we want that?

    // pull the WL direction back to the edge so
    // we don't have to worry about boundaries
    spos = spos - mod(spos, period) + mod(spos - step*rpos, period);

    // double *layer; // chocolate cake
    // layer = (double*)calloc(L*X, sizeof(double));

    double layer[L*X]; // chocolate cake

    double Ub[X];

    // iterators
    int i, j, l;
    for (i = 0; i != bsteps; i++) { 

        for (l = 0; l < L; l++) {
            // step variable determines which space direction we are using
            // for our Wilson line direction
            memcpy(&layer[X*l], &prop[X*(spos+step*l)], sizeof(double)*X);
        }

        if (bdir == 1) {
            for(j = 0; j < L; j++) {
                get_link(Ub,U,D*(spos+step*j)+bmu);
                matmul(&layer[X*j],&layer[X*j],Ub);
            }
        }

        // prepare spos for next step
        spos = spos - mod(spos, bperiod) + mod(spos + bstep*bdir, bperiod);

        if (bdir == -1) {
            for(j = 0; j < L; j++) {
                get_link(Ub,U,D*(spos+step*j)+bmu);
                adjoint(Ub,Ub);
                matmul(&layer[X*j],&layer[X*j],Ub);
            }
        }

        // update propagator
        for (l = 0; l < L; l++) {
            // step variable determines which space direction we are using
            // for our Wilson line direction
            memcpy(&prop[X*(spos+step*l)], &layer[X*l], sizeof(double)*X);
        }
    }

    // free(layer);
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
    /*double mm[6];*/
    /*memcpy(mm, measurement, 6*sizeof(double));*/

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
    list->measurements[num_measurements] = malloc(6*sizeof(double));
    memcpy(list->measurements[num_measurements],measurement,6*sizeof(double));
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

void compute_observables(double *trace_array, int spos, double *U, int n, 
                         double ratio, int t, int b, 
                         double *prop_num, double *prop_den){

    int srcp = X*spos;
    double source[] = {1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0};
    // choose directions for loop
    int mu_r = mu_z;
    int mu_b = mu_x;
    if (mu_b == mu_t) {
        printf("Invalid direction. Exiting.\n");
        exit(1);
    }
    int step = ipow(L,mu_b);
    int period = L*step; //only works if transverse direction is not time
    int tdir, bdir;

    memset(prop_num, 0, sizeof(double)*OLENGTH);
    memcpy(&prop_num[srcp],source,sizeof(double)*X);

    // Numerator
    tpropLMNnew(prop_num, U, n, spos, mu_r, t, tdir=1, ratio);
    spos = mod(spos + LCUBE*t*tdir, LVOLUME);
    trLinksnew(prop_num, U, spos, mu_r, mu_b, b, bdir=1);
    spos = spos - mod(spos, period) + mod(spos+step*b*bdir, period);
    tpropLMNnew(prop_num, U, n, spos, mu_r, t, tdir=-1, ratio);
    // mem copy to RH denominator, 
    // because first half of loop is the same as numerator 
    memcpy(prop_den, prop_num, sizeof(double)*OLENGTH);
    // Finish numerator
    spos = mod(spos + LCUBE*t*tdir, LVOLUME);
    tpropLMNnew(prop_num, U, n, spos, mu_r, t, tdir=1, -ratio);
    spos = mod(spos + LCUBE*t*tdir, LVOLUME);
    trLinksnew(prop_num, U, spos, mu_r, mu_b, b, bdir=-1);
    spos = spos - mod(spos, period) + mod(spos+step*b*bdir, period);
    tpropLMNnew(prop_num, U, n, spos, mu_r, t, tdir=-1, -ratio);

    // trace(trnum, &prop_num[srcp]);
    trace(&trace_array[0], &prop_num[srcp]);
    // printf("Trace: (%.12e) + i(%.12e)\n", trnum[0], trnum[1]);

    // RH denominator

    // set correct bperp position
    bdir = 1;
    spos = spos - mod(spos, period) + mod(spos+step*b*bdir, period);

    spos = mod(spos + LCUBE*t*tdir, LVOLUME);
    tpropLMNnew(prop_den, U, n, spos, mu_r, t, tdir=-1, ratio);
    spos = mod(spos + LCUBE*t*tdir, LVOLUME);
    trLinksnew(prop_den, U, spos, mu_r, mu_b, b, bdir=-1);
    spos = spos - mod(spos, period) + mod(spos+step*b*bdir, period);
    tpropLMNnew(prop_den, U, n, spos, mu_r, t, tdir=1, ratio);

    // trace(trdenA, &prop_den[srcp]);
    trace(&trace_array[2], &prop_den[srcp]);
    // printf("Trace: (%.12e) + i(%.12e)\n", trdenA[0], trdenA[1]);

    //LH denominator
    memset(prop_den, 0, sizeof(double)*OLENGTH);
    memcpy(&prop_den[srcp],source,sizeof(double)*X);

    spos = mod(spos + LCUBE*t*tdir, LVOLUME);

    tpropLMNnew(prop_den, U, n, spos, mu_r, t, tdir=-1, -ratio);
    spos = mod(spos + LCUBE*t*tdir, LVOLUME);
    trLinksnew(prop_den, U, spos, mu_r, mu_b, b, bdir=1);
    spos = spos - mod(spos, period) + mod(spos+step*b*bdir, period);
    tpropLMNnew(prop_den, U, n, spos, mu_r, t, tdir=1, -ratio);
    spos = mod(spos + LCUBE*t*tdir, LVOLUME);
    tpropLMNnew(prop_den, U, n, spos, mu_r, t, tdir=1, -ratio);
    spos = mod(spos + LCUBE*t*tdir, LVOLUME);
    trLinksnew(prop_den, U, spos, mu_r, mu_b, b, bdir=-1);
    spos = spos - mod(spos, period) + mod(spos+step*b*bdir, period);
    tpropLMNnew(prop_den, U, n, spos, mu_r, t, tdir=-1, -ratio);

    trace(&trace_array[4], &prop_den[srcp]);
}

double *read_configuration(char *config) {

    double *V;
    V = (double*)calloc(CLENGTH, sizeof(double));

    char header[1000];
    bzero(header, 600*sizeof(char));
    int header_length = 0;

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
    sprintf(info_path, "%s/info_src%s", inputs[8], inputs[3]);

    struct stat buffer;
    if (stat(info_path, &buffer) != 0) {
        FILE *info_file;
        info_file = fopen(info_path, "w");

        /*fprintf(info_file, header);*/
        fprintf(info_file, "Layout: (uint32_t)(6*double)\n");
        fprintf(info_file, "cols:\tnum.real\tnum.imag\tdenA.real\t"
                           "denA.imag\tdenB.real\tdenB.imag\n");
        fprintf(info_file, "n=%s\n",inputs[1]);
        fprintf(info_file, "r=%s\n",inputs[2]);
        fprintf(info_file, "b=%s\n",inputs[4]);
        fprintf(info_file, "t=%s\n",inputs[5]);
        fprintf(info_file, "c=%s\n",inputs[6]);
        /*fprintf(info_file, "index = b + ncols * t\n");*/

        fclose(info_file);
    }
}

void write_data_to_file(measurement_list *list, char **inputs) {
    // code to write to file
    char data_path[200];
    sprintf(data_path, "%s/n%sr%sb%st%sc%s.bin", inputs[8], inputs[1], 
            inputs[2], inputs[4], inputs[5], inputs[6]);

    FILE *data_file;
    data_file = fopen(data_path, "ab");
    for (int i = 0; i < list->number_of_measurements; i++) {
        fwrite(&list->indices[i], sizeof(uint32_t), 1, data_file);
        fwrite(list->measurements[i], sizeof(double), 6, data_file);
    }
    fclose(data_file);

}

measurement_list *read_data_from_file(char **inputs) {
    measurement_list *list = create_measurement_list(16);
    // code to write to file
    char data_path[200];
    sprintf(data_path, "%s/n%sr%sb%st%sc%s.bin", inputs[8], inputs[1], 
            inputs[2], inputs[4], inputs[5], inputs[6]);

    FILE *data_file = fopen(data_path, "rb");
    if (data_file != NULL) {
        size_t num_read;
        for (int i=0; i<LVOLUME; i++) {
            num_read = 0;
            uint32_t index;
            double measurement[6];
            /*double mm[6];*/
            num_read += fread(&index, sizeof(uint32_t), 1, data_file);
            num_read += fread(measurement, sizeof(double), 6, data_file);
            /*memcpy(mm, measurement, 6*sizeof(double));*/
            if (num_read == 7) {
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

int ssl2sa(char *s, char **a) {
    int c = 0;
    if (a == NULL) {
        int i = 0;
        while (s[i] != '\0') {
            if (s[i] == ' ') {
                c++;
            }
            i++;
        }
    }
    else {
        int i = 0;
        int ri = 0;
        while (s[i] != '\0') {
            if (s[i] == ' ') {
                *(a+c) = calloc(i-ri+1, sizeof(char));
                strncpy(*(a+c), s+ri, i-ri);
                ri = i + 1;
                c++;
            }
            i++;
        }
    }
    return c;
}
