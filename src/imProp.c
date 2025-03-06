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

#define LENGTH(X) (sizeof X / sizeof X[0])

enum dim{mu_x, mu_y, mu_z, mu_t};

// declaration of functions
void fbdelta(double *result, double *layer, double *Ulayer, int lpos);
void tpropLMNnew(double *prop, double *U, int n, int spos, int mu, int tsteps, int tdir, double ratio);
void trLinksnew(double *prop, double *U, int spos, int mu, int bmu, int bsteps, int bdir);
// void build_string(char *out, char *in[], int length);

int main(int argc, char *argv[]) {

    if (argc < 7 || argc > 7) {
        printf("Input should look like: \n");
        printf("%s n ratio num_sources config_number config_path data_dir\n", argv[0]);
        exit(1);
    }

    int n = atoi(argv[1]);
    double ratio = atof(argv[2]);
    int nSources = atoi(argv[3]);
    char *config_number = argv[4];
    char *config = argv[5];
    char *data_dir = argv[6];

    char *header;
    header = (char*)malloc(600*sizeof(char));
    bzero(header, 600*sizeof(char));
    int header_length = 0;

    double *V;
    V = (double*)calloc(CLENGTH, sizeof(double));
    read_nersc(V, config, header, &header_length);

    char read_cksm[20];
    char label[] = "CHECKSUM\0";
    find_value(read_cksm, label, 8, header, header_length);
    printf("Read checksum from header = %s\n", read_cksm);
    uint32_t read_cks = (uint32_t)strtol(read_cksm, NULL, 16);

    uint32_t cks = checksum(V);
    char cksm[20];
    snprintf(cksm, 20, "%x", cks);
    printf("Computed checksum = %s\n", cksm);

    uint32_t cks_diff = read_cks-cks;
    printf("Difference = %u\n", cks_diff);
    if (cks_diff != 0) {
        printf("checksum failed!\n");
        exit(1);
    }

    // clock_t begin = clock();
    double begin = omp_get_wtime();

    // // for gauge invariance check
    //
    // double *Om;
    // Om = (double*)calloc(CLENGTH, sizeof(double));
    // rnd_gauge_tr(Om);
    // transform(U, U, Om);
    
    // srand(4743);

    // make a grid of sources, stuck with this way for now
    // int nSources = 2097152; // nSources < V
    /// Assuming lattice has form T = 2*L.
    int little_L = (int)(ceil(pow((nSources/2),0.25)));
    int little_T = 2*little_L;
    int fact_L = L/little_L;
    int fact_T = T/little_T;

    int tmin = 0, tmax = 11;
    int bmin = 3, bmax = 3;
    // double trnum[2], trdenA[2], trdenB[2];

    // make an info file
    char info_path[200];
    sprintf(info_path, "%s/info", data_dir);

    struct stat buffer;
    if (stat(info_path, &buffer) != 0) {
        FILE *info_file;
        info_file = fopen(info_path, "w");

        fprintf(info_file, header);
        fprintf(info_file, "cols:\tnum.real\tnum.imag\tdenA.real\tdenA.imag\tdenB.real\tdenB.imag\n");
        fprintf(info_file, "bmin=%d\n",bmin);
        fprintf(info_file, "bmax=%d\n",bmax);
        fprintf(info_file, "tmin=%d\n",tmin);
        fprintf(info_file, "tmax=%d\n",tmax);
        fprintf(info_file, "n=%d\n",n);
        fprintf(info_file, "r=%lf\n",ratio);
        fprintf(info_file, "index = b + ncols * t\n");

        fclose(info_file);
    }

    // Create directory structure first, to avoid race condition 
    // in directory creation with parallel code
    for (int i = 0; i < T; i++) {
        char tsubdir[200];
        sprintf(tsubdir, "%s/%d", data_dir, i);
        struct stat st = {0};
        if (stat(tsubdir, &st) == -1) {
            if (mkdir(tsubdir, 0755) == -1) {
                printf("Error: %s\n", strerror(errno));
            }
        }
    }

    // Set up partition of sources for parallelization
    int num_threads = omp_get_max_threads();
    int *sizes;
    sizes = (int*)calloc(num_threads+1, sizeof(int));
    int sizes_extent;
    sizes[0] = 0; // This is sacrosanct
    if (nSources < num_threads || num_threads == 1) {
        sizes_extent = 1;
        sizes[1] = nSources;
    } else {
        sizes_extent = num_threads;

        int high = nSources/sizes_extent + 1;
        int low = nSources/sizes_extent;
        int num_high = nSources%sizes_extent;
        for (int i = 0; i < num_high; i++) {
            sizes[i+1] = sizes[i] + high;
        }
        for (int i = num_high; i < sizes_extent; i++) {
            sizes[i+1] = sizes[i] + low;
        }
    }

    double trace_array_average[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };


    #pragma omp parallel for
    for (int parallel_s = 0; parallel_s < sizes_extent; parallel_s++) {

        double *U;
        U = (double*)calloc(CLENGTH, sizeof(double));
        memcpy(U, V, sizeof(double)*CLENGTH);

        // propagators
        double *prop_num;
        prop_num = (double*)calloc(OLENGTH, sizeof(double));

        double *prop_den;
        prop_den = (double*)calloc(OLENGTH, sizeof(double));

        // stuff for writing data to file
        int file_index;
        double *trace_array;
        int data_length = 6*(bmax-bmin+1)*(tmax-tmin+1); // 3 complex observables -> 6 doubles
        trace_array = (double*)calloc(data_length, sizeof(double));

        for (int s = sizes[parallel_s]; s < sizes[parallel_s+1]; s++) {
            file_index = 0;

            // src[3] = fact_T*(s%little_T);
            // src[2] = fact_L*((s/little_T)%little_L);
            // src[1] = fact_L*(((s/little_T)/little_L)%little_L);
            // src[0] = fact_L*((((s/little_T)/little_L)/little_L)%little_L);

            int src[4];
            src[0] = fact_L*(s%little_L);
            src[1] = fact_L*((s/little_L)%little_L);
            src[2] = fact_L*(((s/little_L)/little_L)%little_L);
            src[3] = fact_T*((((s/little_L)/little_L)/little_L)%little_T);

            int lattice_s = src[0] + L*src[1] + LSQUARE*src[2] + LCUBE*src[3];

            char tsubdir[200];
            sprintf(tsubdir, "%s/%d", data_dir, src[3]);

            char data_path[200];
            // sprintf(data_path, "%s/run_r%s_%s_src%i-%i-%i-%i.bin", tsubdir, argv[2], config_number, src[0], src[1], src[2], src[3]);
            sprintf(data_path, "%s/r%s_c%s_s%i.bin", tsubdir, argv[2], config_number, lattice_s);

            // continue;

            // struct stat buffer;
            if (stat(data_path, &buffer) == 0) {
                printf("%s already exists. Skipping.\n", data_path);
                continue;
            }

            // printf("On thread: %d\n", parallel_s);
            // printf("On source: (%d, %d, %d, %d)\n", src[0], src[1], src[2], src[3]);

            int spos, srcp;
            double source[] = {1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0};
            // choose directions for loop
            int mu_r = mu_z;
            int mu_b = mu_x;
            int step = ipow(L,mu_b);
            int period = L*step; //only works if transverse direction is not time
            int tdir, bdir;

            for (int t = tmin; t < tmax+1; t++) {
                for (int b = bmin; b < bmax+1; b++) {

                    spos = src[0]+L*(src[1]+L*(src[2]+L*src[3]));
                    srcp = X*spos;

                    memset(prop_num, 0, sizeof(double)*OLENGTH);
                    memcpy(&prop_num[srcp],source,sizeof(double)*X);

                    // Numerator
                    tpropLMNnew(prop_num, U, n, spos, mu_r, t, tdir=1, ratio);
                    spos = mod(spos + LCUBE*t*tdir, LVOLUME);
                    trLinksnew(prop_num, U, spos, mu_r, mu_b, b, bdir=1);
                    spos = spos - mod(spos, period) + mod(spos+step*b*bdir, period);
                    tpropLMNnew(prop_num, U, n, spos, mu_r, t, tdir=-1, ratio);
                    // mem copy to RH denominator, because first half of loop is the same as numerator 
                    memcpy(prop_den, prop_num, sizeof(double)*OLENGTH);
                    // Finish numerator
                    spos = mod(spos + LCUBE*t*tdir, LVOLUME);
                    tpropLMNnew(prop_num, U, n, spos, mu_r, t, tdir=1, -ratio);
                    spos = mod(spos + LCUBE*t*tdir, LVOLUME);
                    trLinksnew(prop_num, U, spos, mu_r, mu_b, b, bdir=-1);
                    spos = spos - mod(spos, period) + mod(spos+step*b*bdir, period);
                    tpropLMNnew(prop_num, U, n, spos, mu_r, t, tdir=-1, -ratio);

                    // trace(trnum, &prop_num[srcp]);
                    trace(&trace_array[file_index], &prop_num[srcp]);
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
                    trace(&trace_array[file_index+2], &prop_den[srcp]);
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

                    // trace(trdenB, &prop_den[srcp]);
                    trace(&trace_array[file_index+4], &prop_den[srcp]);
                    // printf("Trace: (%.12e) + i(%.12e)\n", trdenB[0], trdenB[1]);

                    // fprintf(data_file, "%d\t%d\t%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n", 
                    //         index++, t, b, trnum[0], trnum[1], trdenA[0], trdenA[1],
                    //         trdenB[0], trdenB[1]);

                    // increase file_index by the number of columns
                    file_index += 6;

                }
            }

            // write the file all at once so it doesn't just sit there open
            FILE *data_file;
            data_file = fopen(data_path, "wb");
            fwrite(trace_array, sizeof(double), data_length, data_file);
            fclose(data_file);

        }
        free(U);
        free(prop_num);
        free(prop_den);
        free(trace_array);
    }

    free(V);
    free(sizes);
    // free(Om);
    free(header);

    // clock_t end = clock();
    // double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
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

void tpropLMNnew(double *prop, double *U, int n, int spos, int mu, int tsteps, int tdir, double ratio) {

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
    double *layer, *helper, *Ulayer;
    layer = (double*)calloc(L*X, sizeof(double));
    helper = (double*)calloc(L*X, sizeof(double));
    Ulayer = (double*)calloc(L*X, sizeof(double));

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
            // printf("tn=%d\n", tn);
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

    free(layer);
    free(helper);
    free(Ulayer);
}


void trLinksnew(double *prop, double *U, int spos, int mu, int bmu, int bsteps, int bdir) {

    int step = ipow(L,mu);
    int period = L*step;
    int rpos = mod(spos/step, L);
    int bstep = ipow(L,bmu);
    int bperiod = L*bstep; // this will not work if we want the time direction to
                          // be our bdir but why would we want that?

    // pull the WL direction back to the edge so
    // we don't have to worry about boundaries
    spos = spos - mod(spos, period) + mod(spos - step*rpos, period);

    double *layer; // chocolate cake
    layer = (double*)calloc(L*X, sizeof(double));

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

    free(layer);
}


