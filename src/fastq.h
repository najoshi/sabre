#ifndef _FASTQ_H
#define _FASTQ_H

#include "sabre.h"

#define LINE_SIZE 600

typedef struct {
    char name[LINE_SIZE];
    char *comment;
    char seq[LINE_SIZE];
    char other[LINE_SIZE];
    char qual[LINE_SIZE];
    char *r2;
} fq_rec_t;

int get_line(gzFile fq_fd, char *line, int buff);
int get_fq_read(fq_rec_t *fq_rec, gzFile fq_fd);

#endif /*_FASTQ_H*/
