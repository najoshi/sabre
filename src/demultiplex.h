#ifndef DEMULTIPLEX_H
#define DEMULTIPLEX_H

#include "sabre.h"
#include "utils.h"
#include "fastq.h"

#define MAX_READ_SIZE 2048
#define MAX_READ_NUMBER 1000
#define MAX_READ_BUFFER MAX_READ_SIZE*MAX_READ_NUMBER

void* demult_runner(void *arg);

#endif /*DEMULTIPLEX_H*/
