#ifndef _SABRE_H
#define _SABRE_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <zlib.h>
#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include <libgen.h>
#include <time.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef PROGRAM_NAME
#define PROGRAM_NAME "sabre"
#endif

#ifndef AUTHORS
#define AUTHORS "\n\
	         \n      Nikhil Joshi, UC Davis Bioinformatics Core\
                 \n      Kirill Tsyganov, Monash Bioinformatics Platform\
                 \n"
#endif

//https://semver.org/
#define VERSION_MAJOR 0
#define VERSION_MINOR 3
#define VERSION_PATCH 1

#define MAX_BARCODE_LENGTH 100
#define MAX_FILENAME_LENGTH 200

//NOTE more info on struc and typedef
//https://stackoverflow.com/questions/1675351/typedef-struct-vs-struct-definitions
typedef struct actl_bc_cnt {
    char *bc;
    int cnt;
    struct actl_bc_cnt *next;
} actl_bc_cnt;

typedef struct listel_p {
	//TODO this is future implementation
	//char **bc;
	//typedef struct barcode_t *b
	// this is a pointer to an array of pointers
	char **bc;
	char *bc_grp;
	int num_records;
	FILE* bcfile1;
	FILE* bcfile2;
        struct actl_bc_cnt *actl_bc_cnt;
	struct listel_p *next;
} barcode_data_t;

//TODO this is for future implementation
//typdef strcut {
//  char *bc;
//  barcode_data_t *bacrcode_data;
//} barcode_t

typedef struct {
    gzFile fq1_fd;
    gzFile fq2_fd;
    FILE* unassigned1_fd;
    FILE* unassigned2_fd;
    FILE* umis_2_short_fd;
    int mismatch;
    int combine;
    int umi;
    int paired;
    int min_umi_len;
    int max_5prime_crop;
    int no_comment;
} param_t;

typedef struct {
    int num_unknown;
    int total;
} metrics_t;

typedef struct {
    int id;
    const param_t* params;
    barcode_data_t* curr;
    metrics_t* metrics;
    pthread_mutex_t *in_lock, *out_lock;
    pthread_cond_t *cv;
    volatile int *line_num;  // Pointer to line number.  Access to this protected by in_lock
    volatile int *out_line_num;  // Pointer to line number output.  Protected by out_lock
} thread_data_t;


#endif /*_SABRE_H*/
