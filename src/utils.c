#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <libgen.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "sabre.h"
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include "utils.h"

// https://stackoverflow.com/questions/2336242/recursive-mkdir-system-call-on-unix/11425692
// https://stackoverflow.com/questions/7430248/creating-a-new-directory-in-c
const char * _mkdir(const char *file_path) {
    // return straigth away if a file_path is not nested file path
    if(strstr(file_path, "/") == NULL) {
        return file_path;
    }
    //TODO check if directory already exists or not
    //struct stat st = {0};
    //
    //if (stat(dir, &st) == -1) {
    //    mkdir(tmp, S_IRWXU);
    //}

    //char tmp[PATH_MAX]; // can't get this to work..
    char tmp[256];
    char *p = NULL;
    char *dirc = strdup(file_path);
    size_t len;

    const char *dir = dirname(dirc);
    snprintf(tmp, sizeof(tmp),"%s", dir);
    len = strlen(tmp);
    
    if(tmp[len - 1] == '/') {
        tmp[len - 1] = 0;
    }

    for(p = tmp + 1; *p; p++) {
        if(*p == '/') {
           *p = 0;
           mkdir(tmp, S_IRWXU);
           *p = '/';
        }
    }

    mkdir(tmp, S_IRWXU);
    free(dirc);

    return file_path;
}

//NOTE retuns zero on success
//strcmp can be used for sorting, returns pos, zero, neg
//BUT this new implementation can't be used as such just FYI 
int strncmp_with_mismatch (const char *orig_bc, const char *orig_read, size_t mismatch, int max_5prime_crop) {

    int orig_read_len = strlen(orig_read);
    int orig_bc_len = strlen(orig_bc);
    int n_crop = 0;

    if(orig_bc_len > orig_read_len) {
        fprintf (stderr, "Length of the barcode %d is greater than length of the reads %d.", orig_bc_len, orig_read_len);
        return 1;
    }

    while(n_crop <= max_5prime_crop) {

        if(n_crop > orig_read_len) {
            return 1;
        }

        int cnt = 0;
        char u1, u2;
        const char *bc = orig_bc;
        const char *read = orig_read+n_crop;
        int bc_len = orig_bc_len;
    
        while (bc_len-- > 0) {
            u1 = *bc++;
            u2 = *read++;

            if (u1 != u2) {
                cnt++;
                if (cnt > mismatch) {
                    break;
                }
            }

            if (u1 == '\0' || u2 == '\0') {
                return 0;
            }
        }

        if(cnt <= mismatch) {
            return 0;
        }

        n_crop++;
    }
    //this is in the case of error
    return 1;
}

// https://stackoverflow.com/questions/21880730/c-what-is-the-best-and-fastest-way-to-concatenate-strings
//TODO this is a fastq mystrcat function, that returns a pointer to the end of the string
char * get_fqread(kseq_t *fqrec, char *barcode, int no_comment, int remove_seq) {

    size_t fqread_size = 0;

    fqread_size += strlen(fqrec->seq.s);
    fqread_size += (strlen(fqrec->name.s)*2);
    fqread_size += strlen(fqrec->qual.s);
    fqread_size += (strlen(fqrec->comment.s)*2);
    fqread_size += 2;// header signs @ and +
    fqread_size += 2;//two colons (:)
    fqread_size += 4;//cariage returns
    fqread_size += 2;//two spaces

    char *umi = NULL;

    if(barcode[0] != '\0') {
        umi = (char*) malloc( strlen(fqrec->seq.s)-strlen(barcode) + 1 );
        strcpy(umi, (fqrec->seq.s)+strlen(barcode));
        fqread_size += strlen(umi);
    }
    
    char *fqread = (char*) malloc(fqread_size + 1);
    //makes it a zero length string
    fqread[0] = '\0';

    //@READNAME:BACRCODE:UMI
    //1st line
    strcat(fqread, "@");
    strcat(fqread, fqrec->name.s);
    //TODO later can have conditional here depending on the the structure and/or BARCODE/UMI
    if(barcode[0] != '\0') {
        strcat(fqread, ":");
        strcat(fqread, barcode);

        if(umi[0] == '\0') {
            fprintf(stderr, "Error: This shouldn't happened.\n");
            exit(EXIT_FAILURE);
        }

        strcat(fqread, ":");
        strcat(fqread, umi);
        free(umi);
    }

    if(fqrec->comment.l && no_comment == -1) {
        strcat(fqread, " ");
        strcat(fqread, fqrec->comment.s);
    }
    strcat(fqread, "\n");

    //2nd line
    if(remove_seq == 1) {
        strcat(fqread, "N");
    }
    else {
        strcat(fqread, (fqrec->seq.s)+strlen(barcode));
    }
    strcat(fqread, "\n");

    //3rd line
    strcat(fqread, "+");
    strcat(fqread, fqrec->name.s);
    if(fqrec->comment.l && no_comment == -1) {
        strcat(fqread, " ");
        strcat(fqread, fqrec->comment.s);
    }
    strcat(fqread, "\n");

    //4th line
    strcat(fqread, fqrec->qual.s);
    strcat(fqread, "\n");

    return fqread;
}
