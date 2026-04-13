#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <inttypes.h>
#include <time.h>
#include "sizes.h"
#include "fileIO.h"
#include "su3lib.h"

int read_nersc(double *U, char *config, char *header, int *header_length) {
    // need to figure out the issue with fread and doubles for the NERSC config
    // for now the union method works I guess
    union {
        long int retriever;
        double element;
    } package; 

    FILE *config_file = fopen(config, "rb");
    int safety = 2000;
    char check_end[] = "END_HEADER\n";
    *header_length = 0;
    for (int i = 0; i < safety; i++) {
        header[i] = fgetc(config_file);
        *header_length += 1;
        if ( i > 10) {
            if (strncmp(check_end, &header[i-10], sizeof(check_end)) == 0) {
                break;
            }
        }
    }

    // read in configuration data
    for (size_t i = 0; i != CLENGTH; i++) {
        fread(&package.retriever, sizeof(package.retriever), 1, config_file);
        package.retriever = be64toh(package.retriever);
        U[i] = package.element;
    }
    fclose(config_file);

    char read_cksm[20];
    char label[] = "CHECKSUM\0";
    find_nersc_value(read_cksm, label, 8, header, *header_length);
    printf("Read checksum from header = %s\n", read_cksm);
    int read_cks = (int)strtol(read_cksm, NULL, 16);

    int cks = checksum(U);
    char cksm[20];
    snprintf(cksm, 20, "%x", cks);
    printf("Computed checksum = %s\n", cksm);

    int cks_diff = read_cks-cks;
    printf("Difference = %u\n", cks_diff);
    if (cks_diff != 0) {
        printf("checksum failed!\n");
        return 1;
    }

    return 0;
}

void write_nersc(double *U, char *out_file, char *header, int header_length) {
    union {
        long int sender;
        double element;
    } package; 

    FILE *config_file = fopen(out_file, "wb");
    fwrite(header, sizeof(char), header_length, config_file);
    // fwrite(U, sizeof(double), CLENGTH, config_file);

    for (size_t i = 0; i != CLENGTH; i++) {
        package.element = U[i];
        package.sender = htobe64(package.sender);
        fwrite(&package.sender, sizeof(package.sender), 1, config_file);
    }
    fclose(config_file);
}

void make_nersc_header(char *header, int *header_length, double *U, char* label) {
    // First gather all the information
    // Ultimately this should take a structure as an argument with
    // all the relavant information

    int cks = checksum(U);
    char cksm[20];
    snprintf(cksm, 20, "%x", cks);

    double tr;
    link_trace(&tr, U);
    char ltr[30];
    snprintf(ltr, 30*sizeof(char), "%.14lf", tr);

    double tplaq, splaq, aplaq;
    all_plaquettes(U, &aplaq, &splaq, &tplaq);
    char plaq[30];
    snprintf(plaq, 30*sizeof(char), "%.14lf", aplaq);

    /*char label[30];*/
    /*snprintf(label, 30*sizeof(char), "T%dL%d_ape", T, L);*/

    // char *user = getlogin();
    char user[10];
    snprintf(user, 10*sizeof(char), "wayne");

    time_t t = time(NULL);
    struct tm tms = *localtime(&t);
    char today[80];
    snprintf(today, 80*sizeof(char), 
             "%d-%02d-%02d %02d:%02d:%02d", 
             tms.tm_year + 1900, tms.tm_mon + 1, 
             tms.tm_mday, tms.tm_hour, tms.tm_min, tms.tm_sec);

    char *inputs[] = {
        cksm,
        ltr,
        plaq,
        label,
        user,
        today
    };

    // Now make the header
    char new_header[1000];
    // char line[100];
    int old_index = 0, new_index = 0;
    // char line_end[] = "END_HEADER\n";
    char *lines[] = {
        "CHECKSUM\0",
        "LINK_TRACE\0",
        "PLAQUETTE\0",
        "ENSEMBLE_LABEL\0",
        "CREATOR_MACHINE\0",
        "CREATION_DATE\0"
    };
    while( old_index < *header_length ){
        new_header[new_index] = header[old_index];
        for (int j = 0; j < 6; j++ ) {
            int prop_size = 0;
            while (lines[j][prop_size] != '\0') {
                prop_size++;
            };
            if (strncmp(&header[old_index], lines[j], prop_size*sizeof(char)) == 0) {
                memcpy(&new_header[new_index], lines[j], prop_size*sizeof(char));
                new_index += prop_size;
                memcpy(&new_header[new_index], " = ", 3*sizeof(char));
                new_index += 3;
                int i = 0;
                while (inputs[j][i] != '\0') {
                    new_header[new_index] = inputs[j][i];
                    i++;
                    new_index++;
                };
                new_header[new_index] = '\n';
                while (header[old_index] != '\n') {
                    old_index++;
                };
            };
        };
        new_index++;
        old_index++;
    };
    new_header[new_index] = '\0';
    // // bzero(header, sizeof(char)*(*header_length));
    *header_length = new_index;
    memcpy(header, new_header, sizeof(char)*(*header_length+1));
}

void find_nersc_value(char *value, char *label, int label_length, char *header, int header_length) {
    for (int i = 0; i < header_length; i++) {
        if (strncmp(&header[i], label, label_length*sizeof(char)) == 0) {
            i = i + label_length + 3; // the 3 is for ' = '
            int j = 0;
            while (header[i] != '\n') {
                value[j] = header[i];
                i++;
                j++;
            }
            value[j] = '\0';
            break;
        }
    }
}

int read_scidac(double *U, char *config, char *header, int *header_length) {
    // this is nothing, what is this, fix this
    union {
        long int retriever;
        double element;
    } package; 

    FILE *config_file = fopen(config, "rb");
    int safety = 2000;
    char check_end[] = "END_HEADER\n";
    *header_length = 0;
    for (int i = 0; i < safety; i++) {
        header[i] = fgetc(config_file);
        *header_length += 1;
        if ( i > 10) {
            if (strncmp(check_end, &header[i-10], sizeof(check_end)) == 0) {
                break;
            }
        }
    }

    // read in configuration data
    for (size_t i = 0; i != CLENGTH; i++) {
        fread(&package.retriever, sizeof(package.retriever), 1, config_file);
        package.retriever = be64toh(package.retriever);
        U[i] = package.element;
    }
    fclose(config_file);

    char read_cksm[20];
    char label[] = "CHECKSUM\0";
    find_nersc_value(read_cksm, label, 8, header, *header_length);
    printf("Read checksum from header = %s\n", read_cksm);
    int read_cks = (int)strtol(read_cksm, NULL, 16);

    int cks = checksum(U);
    char cksm[20];
    snprintf(cksm, 20, "%x", cks);
    printf("Computed checksum = %s\n", cksm);

    int cks_diff = read_cks-cks;
    printf("Difference = %u\n", cks_diff);
    if (cks_diff != 0) {
        printf("checksum failed!\n");
        return 1;
    }

    return 0;

}

int checksum(double *U) {
    // crc checksum adapted from Hudspith's GLU code

    int res = 0;
    for (size_t i = 0; i < CLENGTH; i++) {
        int *buf = (int*)(U + i);
        res += *buf + *(buf+1);
    }

    return res;
}

