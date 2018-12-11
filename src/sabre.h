#ifndef SABRE_H
#define SABRE_H

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

/* Options drawn from GNU's coreutils/src/system.h */
/* These options are defined so as to avoid conflicting with option
values used by commands */
enum {
  GETOPT_HELP_CHAR = (CHAR_MIN - 2),
  GETOPT_VERSION_CHAR = (CHAR_MIN - 3)
};
#define GETOPT_HELP_OPTION_DECL \
"help", no_argument, NULL, GETOPT_HELP_CHAR
#define GETOPT_VERSION_OPTION_DECL \
"version", no_argument, NULL, GETOPT_VERSION_CHAR
#define case_GETOPT_HELP_CHAR(Usage_call) \
case GETOPT_HELP_CHAR: \
Usage_call(EXIT_SUCCESS); \
break;
#define case_GETOPT_VERSION_CHAR(Program_name, Version, Authors) \
case GETOPT_VERSION_CHAR: \
fprintf(stdout, "%s version %0.3f\nCopyright (c) 2011 The Regents " \
"of University of California, Davis Campus.\n" \
"%s is free software and comes with ABSOLUTELY NO WARRANTY.\n"\
"Distributed under the MIT License.\n\nWritten by %s\n", \
Program_name, Version, Program_name, Authors); \
exit(EXIT_SUCCESS); \
break;
// end code drawn from system.h

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

#define barcode_destroy(barcode_data_t *bc_data)\
            gzclose(bc_data->bcfile1)\
            gzclose(bc_data->bcfile2)\
            free(bc_data->bc)

/* Function Prototypes */
int paired_main (int argc, char *argv[]);

typedef struct barcodes_t {
    char *bc;
    int cnts;
} barcodes_t;

typedef struct umis_t {
    char *bc;
    int len;
    int cnts;
} umis_t;

KSEQ_INIT(gzFile, gzread)

typedef struct {
    gzFile fq1_fd=NULL; //gzFile is file descriptor - fd //https://en.wikipedia.org/wiki/File_descriptor
    gzFile fq2_fd=NULL;
    FILE* umis_2_short_fd=NULL;
    int mismatch=0;
    int threads=4;
    int combine = -1;
    int umi = -1;
    int paired = -1;
    int min_umi_len=0;
    int max_5prime_crop=0;
    int no_comment = -1;
} param_t;

typedef struct {
    const *param_t params;	
    const *barcode_data_t curr;
    int id;
    pthread_mutex_t *in_lock, *out_lock;
    pthread_cond_t *cv;
    volatile int *line_num;  // Pointer to line number.  Access to this protected by in_lock
    kseq_t fqrec1;
    kseq_t fqrec2;
} thread_data_t;

#define param_destroy(param_t *param)\
            gzclose(param->fq1_fd)\
            gzclose(param->fq2_fd)\
            fclose(param->umis_2_short_fd)

typedef struct {
    int num_unknown=0;
    int total=0;
} metrics_t;
//This is needed if compilling with -std=c99, read below for more
//https://stackoverflow.com/questions/26284110/strdup-confused-about-warnings-implicit-declaration-makes-pointer-with
char *strdup(const char*);
char *strndup(const char *s, size_t n);

const char * _mkdir (const char *dir);
int chk_bc_mtch(const char *orig_bc, const char *orig_read, size_t mismatch, int max_5prime_crop);
void get_fqread(char **fqread, kseq_t *fqrec, char *barcode, char *umi_idx, int no_comment, int n_crop);
void get_merged_fqread(char **fqread, kseq_t *fqrec1, kseq_t *fqrec2, char *barcode, char *umi_idx, int no_comment, int n_crop);
void get_bc_fn(char **bcout_fn, char *s_name, char *barcode, int read_type);

typedef struct {
    int id;
    pthread_mutex_t *in_lock, *out_lock;
    pthread_cond_t *cv;
    volatile int *line_num;  // Pointer to line number.  Access to this protected by in_lock
    volatile int *out_line_num;  // Pointer to line number output.  Protected by out_lock
    kseq_t *fqrec1;
    kseq_t *fqrec2;
} thread_data;

void* demult_runner(void *arg);

// usage.c
void usage(int status);
void little_story(int status);
void version(int status);

#endif /*SABRE_H*/
