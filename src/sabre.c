#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <limits.h>
#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include "sabre.h"

void main_usage (int status) {

    fprintf (stdout,
             "\n\
              \n  Usage: %s <command> [options]\
              \n\
              \n  Command:\
              \n\
              \n    se\tsingle-end barcode de-multiplexing\
              \n    pe\tpaired-end barcode de-multiplexing\
              \n\
              \n    --help\tto get more help\
              \n    --version\tprint current version to stdout\
              \n\
              \n  Info: Sabre is a heavy cavalry sword with a curved blade and a single cutting edge\
              \n        Not sure though if the meaning was intended by original author...\
              \n\
              \n",
              PROGRAM_NAME);

    exit (status);
}

int main (int argc, char *argv[]) {
    int retval=0;

    if (argc < 2 || (strcmp (argv[1],"pe") != 0 && strcmp (argv[1],"se") != 0 && strcmp (argv[1],"--version") != 0 && strcmp (argv[1],"--help") != 0)) {
        main_usage (EXIT_FAILURE);
    }

    if (strcmp (argv[1],"--version") == 0) {
        fprintf(stdout,
                "\n\
                 \n   %s version %0.3f\
                 \n\
                 \n   Copyright (c) 2011 The Regents of University of California, Davis Campus.\
                 \n   %s is free software and comes with ABSOLUTELY NO WARRANTY.\
                 \n   Distributed under the MIT License.\
                 \n\
                 \n   Written by %s\
                 \n\
                 \n",
                PROGRAM_NAME, VERSION, PROGRAM_NAME, AUTHORS);
            exit (EXIT_SUCCESS);
    }

    else if (strcmp (argv[1],"--help") == 0) {
        main_usage (EXIT_SUCCESS);
    }

    else if (strcmp (argv[1],"pe") == 0) {
        retval = paired_main (argc, argv);
        return (retval);
    }

    return 0;
}
