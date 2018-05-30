#ifndef UTILS_H
#define UTILS_H

#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

//This is needed if compilling with -std=c99, read below for more
//https://stackoverflow.com/questions/26284110/strdup-confused-about-warnings-implicit-declaration-makes-pointer-with
char *strdup(const char*);
char *strndup(const char *s, size_t n);

const char * _mkdir (const char *dir);
int chk_bc_mtch(const char *orig_bc, const char *orig_read, size_t mismatch, int max_5prime_crop);
void get_fqread(char **fqread, kseq_t *fqrec, char *barcode, char *umi_idx, int no_comment, int n_crop);
void get_merged_fqread(char **fqread, kseq_t *fqrec1, kseq_t *fqrec2, char *barcode, char *umi_idx, int no_comment, int n_crop);
void get_bc_fn(char **bcout_fn, char *s_name, char *barcode, int read_type);

#endif /*UTILS_H*/
