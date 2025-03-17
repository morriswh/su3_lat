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

#define LENGTH(X) (sizeof X / sizeof X[0])

enum dim{mu_x, mu_y, mu_z, mu_t};

// declarations
void segment(double *loop, double *U, int *pos, int mu, int comp);
void diagonal_segment(double *loop, double *U, int *pos, 
                      int mu_i, int mu_j, int icomp, int jcomp,
                      int path_index);
int compare(int x, int y, int xsteps, int ysteps); // for "smart" pathing
void source_grid(int *pos, int source, int num_sources); // commenting this doesn't give an implicit declaration error. why?
void stoia(int *output, char *input);

const int paths[][7] = {
    {0, 0, 0, 1, 1, 1, 1},
    {1, 1, 1, 1, 0, 0, 0},
    {1, 0, 1, 0, 1, 0, 1},
    {1, 1, 0, 0, 1, 1, 0},
    {0, 1, 1, 0, 0, 1, 1}
};

int main(int argc, char *argv[]) {

    if (argc == 1) {
        printf("Input should look like: \n");
        printf("%s config_number config_path data_dir num_sources bmin_bmax shape onrails\n", argv[0]);
        exit(1);
    }

    char *config_number = argv[1];
    char *config = argv[2];
    char *data_dir = argv[3];
    int nSources = atoi(argv[4]); // max sources: 2097152
    int bmm[2];
    stoia(bmm, argv[5]);
    int shape[4];
    stoia(shape, argv[6]);
    int onrails = atoi(argv[7]);


    char *header;
    header = (char*)malloc(600*sizeof(char));
    bzero(header, 600*sizeof(char));
    int header_length = 0;

    double *U;
    U = (double*)calloc(CLENGTH, sizeof(double));
    read_nersc(U, config, header, &header_length);

    // // for gauge invariance check
    //
    // double *Om;
    // Om = (double*)calloc(CLENGTH, sizeof(double));
    // rnd_gauge_tr(Om);
    // transform(U, U, Om);

    char data_path[200];
    sprintf(data_path, "%s/c%s.bin", data_dir, config_number);

    // struct stat buffer;
    // if (stat(data_path, &buffer) == 0) {
    //     printf("%s already exists. Skipping.\n", data_path);
    //     exit(1);
    // }

    clock_t begin = clock();

    int num_paths = 5;
    if (onrails == 1) {
        num_paths = 1;
    }

    // #pragma omp parallel for
    for( int bi = bmm[0]; bi < bmm[1]; bi++) {
        for (int pi = 0; pi < num_paths; pi++) {
            for (int s = 0; s < nSources; s++) {
                printf("What\n");

                // butterfly loop with one wing along z and one wing along y direction
                // declare and then initialize
                double loop[X];
                memset(loop, 0, sizeof(loop));
                loop[0] = 1;  
                loop[8] = 1;  
                loop[16] = 1;

                int pos;
                source_grid(&pos, s, nSources);

                // argument 'pi' is meangingless if we only move along grid directions
                diagonal_segment(loop, U, &pos, 1, 2, shape[0], shape[1], pi);
                segment(loop, U, &pos, 0, bi);
                diagonal_segment(loop, U, &pos, 1, 2, -shape[0], -shape[1], pi);
                diagonal_segment(loop, U, &pos, 1, 2, shape[2], shape[3], pi);
                segment(loop, U, &pos, 0, -bi);
                diagonal_segment(loop, U, &pos, 1, 2, -shape[2], -shape[3], pi);

                double loop_trace;
                trace(&loop_trace, loop);

                FILE *data_file;
                data_file = fopen(data_path, "ab+");
                fwrite(&loop_trace, sizeof(double), 1, data_file);
                fclose(data_file);

            }
        }
    }

    free(U);
    free(header);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
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
                      int mu_i, int mu_j, int icomp, int jcomp,
                      int path_index) {

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
        // int result = compare(i, j, isteps, jsteps);
        int result = paths[path_index][k];
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

void source_grid(int *pos, int source, int num_sources) {

    const int little_L = (int)(ceil(pow((nSources/2),0.25)));
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


