#include <stdio.h>
#include <string.h>
#include "utils.h"

/*
gcc -Wall -O2 -std=c99 -o metrics metrics.c -lz
*/

#define BARCODE_ARRAY 1000000

//int chk_bc_arr(char *arr, char *bc);
//
//int chk_bc_arr(char *arr, char *bc) {
//
//    int i = 0;
//    while(arr[i] != '\0') {
//
//        if(strcmp(bc, arr[i]) == 0) {
//            return i; 
//	}
//
//        i += 1;
//    }
//
//    //for(int i; i < strlen(arr); i++) {
//    //    if(strcmp(bc, arr[i]) == 0) {
//    //        return i; 
//    //    }
//    //}
//
//    return -1;
//}

int main (int argc, char *argv[]) {

    gzFile pe1=NULL;
    pe1 = gzopen(argv[1], "r");
    kseq_t *fqrec1;
    int l1;

    fqrec1 = kseq_init(pe1);

    char barcodes[BARCODE_ARRAY];
    //barcodes[0] = '\0';
    int bc_cnts[BARCODE_ARRAY];
    //bc_cnts[0] = '\0';
    int loc = 0; 
    /* Get reads, one at a time */
    while ((l1 = kseq_read(fqrec1)) >= 0) {
        //fprintf(stdout, "%s\n", fqrec1->name.s);

	char *last;
	char *last2;

	char *p = strtok(fqrec1->name.s, ":");
        //fprintf(stdout, "check %s\n", p);
        while (p != NULL) {
            //fprintf(stdout, "%s\n", p);
	    last2 = last;
	    last = p;
            p = strtok(NULL, ":");
	}

	//fprintf(stdout, "%s %zu\n", last2, strlen(last2));
        //int bc_idx = chk_bc_arr(barcodes, last2);
        int idx = 0;
        while(barcodes[idx] != '\0') {
	    fprintf(stdout, "%s %zu\n", barcodes[idx], strlen(barcodes[idx]));
        
        }

        //if(bc_idx >= 0) {
        //    bc_cnts[bc_idx] += 1;
        //}
        //else {
        //    barcodes[loc] = last2;
        //    loc += 1;
        //}

        //fprintf(stdout, "%s\n", last2);
    }

    //for(int i; i < strlen(barcodes); i++) {
    //    fprintf(stdout, "%s %d\n", barcodes[i], bc_cnts[i]); 
    //}

    gzclose(pe1);
    kseq_destroy(fqrec1);

    return EXIT_SUCCESS;
}
