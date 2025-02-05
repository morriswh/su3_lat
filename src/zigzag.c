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
void segment(double *loop, double *U, int *pos, int mu, int steps, int dir);
void diagonal_segment(double *loop, double *U, int *pos, 
                      int mu_i, int mu_j, int isteps, int jsteps,
                      int idir, int jdir);
int compare(int x, int y, int xsteps, int ysteps);

int main(int argc, char *argv[]) {

    int nSources = 2097152; // nSources < V

    if (argc == 1) {
        printf("Input should look like: \n");
        printf("%s config_number config_path data_dir\n", argv[0]);
        exit(1);
    }

    char *config_number = argv[1];
    char *config = argv[2];
    char *data_dir = argv[3];

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

    // #pragma omp parallel for
    for (int s = 0; s < nSources; s++) {

        // butterfly loop with one wing along z and one wing along y direction
        // declare and then initialize
        double loop[X];
        memset(loop, 0, sizeof(loop));
        loop[0] = 1;  
        loop[8] = 1;  
        loop[16] = 1;

        int pos = s;

        diagonal_segment(loop, U, &pos, 1, 2, 5, 0, 1, 1);
        segment(loop, U, &pos, 0, 3, 1);
        diagonal_segment(loop, U, &pos, 1, 2, 5, 0, -1, -1);
        diagonal_segment(loop, U, &pos, 1, 2, 0, 5, 1, -1);
        segment(loop, U, &pos, 0, 3, -1);
        diagonal_segment(loop, U, &pos, 1, 2, 0, 5, -1, 1);

        double loop_trace;
        trace(&loop_trace, loop);

        FILE *data_file;
        data_file = fopen(data_path, "ab+");
        fwrite(&loop_trace, sizeof(double), 1, data_file);
        fclose(data_file);

    }

    free(U);
    free(header);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Elapsed time: %lfs\n", time_spent);

    return 0;
}

void segment(double *loop, double *U, int *pos, int mu, int steps, int dir) {

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
                      int mu_i, int mu_j, int isteps, int jsteps,
                      int idir, int jdir) {

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


int compare(int x, int y, int xsteps, int ysteps) {

    double xratio = (double)((x-0.5)/xsteps);
    double yratio = (double)((y-0.5)/ysteps);

    return xratio > yratio;
}
