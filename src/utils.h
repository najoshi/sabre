#ifndef UTILS_H
#define UTILS_H

#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

//This is needed if compilling with -std=c99, read below for more
//https://stackoverflow.com/questions/26284110/strdup-confused-about-warnings-implicit-declaration-makes-pointer-with
char *strdup(const char*);

int strncmp_with_mismatch (const char *orig_bc, const char *orig_read, size_t mismatch, int max_5prime_crop);
const char * _mkdir (const char *dir);
char * get_fqread(kseq_t *fqrec, char *barcode, int no_comment, int remove_seq);

#endif /*UTILS_H*/
