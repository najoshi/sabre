#ifndef FASTQ_H
#define FASTQ_H

#include "sabre.h"

#define LINE_SIZE 600

typedef struct {
    char name[LINE_SIZE];
    char *comment;
    char seq[LINE_SIZE];
    char other[LINE_SIZE];
    char qual[LINE_SIZE];
} fq_rec_t;

void init_fq_rec(fq_rec_t *fq_rec);
int get_line(gzFile fq_fd, char *line, int buff);
int get_fq_rec(fq_rec_t *fq_rec, gzFile fq_fd);

#endif /*FASTQ_H*/
