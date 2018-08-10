#include <stdio.h>
#include <string.h>
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

int chk_bc_arr(char *arr[], char *bc) {

    int i = 0;
    while(arr[i] != 0) {
        if(strcmp(bc, arr[i]) == 0) {
            return i;
    }

        i += 1;
    }

    arr[i] = strdup(bc);

    // error handling
    if(i == BARCODE_ARRAY-1) {
        fprintf(stderr, "ERROR: gone too far in the array\n");
        exit(1);
    }

    // could also do initialisation of the next element here
    // arr[i+1]=0;
    return i;
}

int main (int argc, char *argv[]) {

    gzFile pe1=NULL;
    pe1 = gzopen(argv[1], "r");
    kseq_t *fqrec1;
    int l1;

    fqrec1 = kseq_init(pe1);

    /*
      WISDOM these tow are actually different things
      char barcodes[BARCODE_ARRAY][BARCODE];
      char *barcodes[];
     */

    char **barcodes = calloc(BARCODE_ARRAY, sizeof(char*));
    int *bc_cnts = calloc(BARCODE_ARRAY, sizeof(int));

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

        int bc_idx = chk_bc_arr(barcodes, last2);
        bc_cnts[bc_idx] += 1;
    }

    {
        /*
	  WISDOM this is to limit the right scope for i
	 */

        int i = 0;
        while(barcodes[i] != 0) {
            fprintf(stdout, "%s %d\n", barcodes[i], bc_cnts[i]);
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
//typedef struct barcodes_t {
//    char *barcodes;
//    int cnts;
//} barcodes_t;
//
//
//barcodes_t barcodes[size];
//
//barcodes[i].cnts +=1;
