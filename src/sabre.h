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
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

#ifndef PROGRAM_NAME
#define PROGRAM_NAME "sabre"
#endif

#ifndef AUTHORS
#define AUTHORS "\n\
	         \n      Nikhil Joshi, UC Davis Bioinformatics Core\
                 \n      Kirill Tsyganov, Monash Bioinformatics Platform\
                 \n"
#endif

#ifndef VERSION
//https://semver.org/
#define VERSION_MAJOR 0
#define VERSION_MINOR 3
#define VERSION_PATCH 1
#endif

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
	char* bc;
	int num_records;
	gzFile bcfile1;
	gzFile bcfile2;
        struct actl_bc_cnt *actl_bc_cnt;
	struct listel_p *next;
} barcode_data_t;

typedef struct {
    gzFile fq1_fd; //gzFile is file descriptor - fd //https://en.wikipedia.org/wiki/File_descriptor
    gzFile fq2_fd;
    gzFile unassigned1_fd; 
    gzFile unassigned2_fd; 
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
