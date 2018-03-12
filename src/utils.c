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
    char *dirc;
    dirc = strdup(file_path);
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
int strncmp_with_mismatch (const char *bc, const char *read, register size_t bc_len, register size_t mismatch, int max_5prime_crop) {

    register char u1, u2;
    int cnt=0;
    int n_crop=0;

    char *orig_bc = strdup(bc);
    char *orig_read = strdup(read);

    while(n_crop <= max_5prime_crop) {
        bc = orig_bc;
        read = orig_read+n_crop;
    
        while (bc_len-- > 0) {
            u1 = *bc++;
            u2 = *read++;

            if (u1 != u2) {
                cnt++;
                if (cnt > mismatch) {
                    //return u1 - u2;
                    break;
                }
            }
            if (u1 == '\0' || u2 == '\0') {
                return 0;
            }
        }
        n_crop++;
    }
    return 0;
}
