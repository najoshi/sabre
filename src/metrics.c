#include "utils.h"

/*
   gcc -Wall -O2 -std=c99 -o metrics metrics.c -lz
 */

#define BARCODE_ARRAY 1000000

//int chk_bc_arr(const char *arr, char *bc);

/*
   WISDOM char **arr and char *arr[] mean the same things
   but later one more informative
 */

int chk_bc_arr(barcodes_t *arr, char *bc) {

    int i = 0;
    while(arr[i].bc != 0) {

        if(strcmp(bc, arr[i].bc) == 0) {
            return i;
        }

        i += 1;
    }

    arr[i].bc = strdup(bc);

    // error handling
    if(i == BARCODE_ARRAY-1) {
        fprintf(stderr, "ERROR: gone too far in the array\n");
        exit(1);
    }

    // could also do initialisation of the next element here
    // arr[i+1]=0;
    return i;
}

int chk_umi_arr(umis_t *arr, char *bc) {

    int i = 0;
    while(arr[i].bc != 0) {

        if(strcmp(bc, arr[i].bc) == 0) {
            return i;
        }

        i += 1;
    }

    arr[i].bc = strdup(bc);
    arr[i].len = strlen(bc);

    if(i == BARCODE_ARRAY-1) {
        fprintf(stderr, "ERROR: gone too far in the array\n");
        exit(1);
    }
    return i;
}

int umi_n_cmp(const void *a1, const void *b1) {

    const umis_t* a=a1;
    const umis_t* b=b1;

    return a->cnts - b->cnts;
}

int bc_n_cmp(const void *a1, const void *b1) {

    const barcodes_t* a=a1;
    const barcodes_t* b=b1;

    return a->cnts - b->cnts;
}

int main (int argc, char *argv[]) {

    if (argc <= 2) {
        fprintf(stderr, "\n\
                         \n  Usage:   metrics <fastq_file> <mode>\
                         \n\
                         \n  Options:\
                         \n\
                         \n           mode INT     [0|1]; 0 = sample barcode, 1 = umis barcodes\
                         \n\
                         \n");
        exit (EXIT_SUCCESS);
    }

    gzFile pe1=NULL;
    pe1 = gzopen(argv[1], "r");
    kseq_t *fqrec1;
    int l1;
    int mode = atoi(argv[2]);

    fqrec1 = kseq_init(pe1);

    /*
       WISDOM these tow are actually different things
       char barcodes[BARCODE_ARRAY][BARCODE];
       char *barcodes[];
     */

    //char **barcodes = calloc(BARCODE_ARRAY, sizeof(char*));
    //int *bc_cnts = calloc(BARCODE_ARRAY, sizeof(int));

    barcodes_t *barcodes;
    barcodes = calloc(BARCODE_ARRAY, sizeof(barcodes_t));

    umis_t *umis;
    umis = calloc(BARCODE_ARRAY, sizeof(umis_t));

    /* Get reads, one at a time */
    while ((l1 = kseq_read(fqrec1)) >= 0) {

        char *last;
        char *last2;

        char *p = strtok(fqrec1->name.s, ":");

        while (p != NULL) {
            last2 = last;
            last = p;
            p = strtok(NULL, ":");
        }

        if(mode == 0) {
            int bc_idx = chk_bc_arr(barcodes, last2);
            barcodes[bc_idx].cnts += 1;
        }
        else if(mode == 1) {

            int umi_idx = chk_umi_arr(umis, last);
            umis[umi_idx].cnts += 1;
        }
    }

    if(mode == 0) {

        int n = 0;
        while(barcodes[n].bc != 0) {
            n++;
        }
        qsort(barcodes, n, sizeof(barcodes_t), bc_n_cmp);

        /*
           WISDOM this is to limit the right scope for i
           */

        int i = 0;
        while(barcodes[i].bc != 0) {
            fprintf(stdout, "%s\t%d\n", barcodes[i].bc, barcodes[i].cnts);
            i++;
        }
    }
    else if(mode == 1) {

        int n = 0;
        while(umis[n].bc != 0) {
            n++;
        }
        qsort(umis, n, sizeof(umis_t), umi_n_cmp);

        int i = 0;
        while(umis[i].bc != 0) {
            fprintf(stdout, "%s\t%d\t%d\n", umis[i].bc, umis[i].len, umis[i].cnts);
            i++;
        }
    }

    gzclose(pe1);
    kseq_destroy(fqrec1);

    return EXIT_SUCCESS;
}
//TODO use struct instead
//used qsort
//write a function to sort  e.g bubbleSort
//actually sorting function should be simpler then bubbleSort
//it should return -1, 0 or 1. OR -1 or 1.
//
//barcodes_t barcodes[size];
//
//barcodes[i].cnts +=1;
