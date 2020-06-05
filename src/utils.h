#ifndef UTILS_H
#define UTILS_H

#include "sabre.h"
#include "fastq.h"

typedef struct barcodes_t {
    char *bc;
    int cnts;
} barcodes_t;

typedef struct umis_t {
    char *bc;
    int len;
    int cnts;
} umis_t;

typedef struct {
    int mismatches;
    int cropped;
} match_ret_t;

//This is needed if compilling with -std=c99, read below for more
//https://stackoverflow.com/questions/26284110/strdup-confused-about-warnings-implicit-declaration-makes-pointer-with
// char *strdup(const char*);
// char *strndup(const char *s, size_t n);

const char * _mkdir (const char *dir);
match_ret_t chk_bc_mtch(const char *orig_bc, const char *orig_read, size_t max_mismatch, int max_5prime_crop);
void get_fqread(char *fqread, fq_rec_t *fq_rec, const char *barcode, char *umi_idx, int no_comment, int n_crop);
void get_merged_fqread(char *fqread, fq_rec_t *fq_rec1, fq_rec_t *fq_rec2, const char *barcode, char *umi_idx, int no_comment, int n_crop);
void get_bc_fn(char **bcout_fn, char *s_name, char *barcode, int read_type, int gz);

void set_default_params(param_t *params);
void params_destroy(param_t *params);

#endif /*UTILS_H*/
