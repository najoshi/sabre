#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <libgen.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "sabre.h"

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
