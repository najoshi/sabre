
#include "usage.h"

void usage(int status) {

    fprintf(stderr, "\n  Usage: %s [OPTIONS] -f <fastq_R1> -r <fastq_R2> -b <barcode_file>\
                     \n\
                     \n\
                     \n  Options:\
                     \n\
                     \n    Required:\
                     \n\
                     \n        -f, --fq1 FILE                Input FASTQ R1 read\
                     \n        -r, --fq2 FILE                Input FASTQ R2 reads\
                     \n        -b, --barcodes FILE           Barcodes files, one barcode per line, e.g BC\\tPREFIX\
                     \n        -w, --unassigned CHAR         Unassigned prefix\
                     \n\
                     \n    Other:\
                     \n\
                     \n        -c, --combine                 Combine R1 and R2 [NULL]\
                     \n        -u, --umi                     Indicates that umi present in the R1 read [NULL]\
                     \n        -m, --max-mismatch INT        Maximum number of mismatches allowed in a barcode [0]\
                     \n        -l, --min-umi-len INT         Minimum UMI length to keep [0]\
                     \n        -a, --max-5prime-crop INT     Maximum number of possible bases cropped from 5prime [0]\
                     \n        -n, --no-comment              Drop extra comments from FASTQ header [NULL]\
                     \n        -s, --stats FILE              Write stats to file instead of STDOUT [STDOUT]\
                     \n\
                     \n    Extras:\
                     \n\
                     \n        -t, --threads INT             specify number of threads to use [4]\
                     \n        -v, --version                 get current version\
                     \n        -h, --hel                     get help menu, exit status is zero\
                     \n        -o, --story                   little story about sabre tool\
                     \n\
                     \n",
                     PROGRAM_NAME);

    exit(status);
}

void little_story(int status) {

    fprintf(stdout, "\n\
                     \n  Little story:\
                     \n\
                     \n        Sabre is a heavy cavalry sword with a curved blade and a single cutting edge\
                     \n        Not sure though if the meaning was intended by original author...\
                     \n\
                     \n        Later on I was pointed out to me that yes of course it was intended\
                     \n        since we are cutting off adaptors..\
                     \n\
                     \n        to be continued...\
                     \n\
                     \n");

    exit(status);
}

void version(int status) {

    fprintf(stdout, "\n\
                     \n   %s\
                     \n\
		     \n   version: %d.%d.%d\
                     \n\
                     \n   Copyright (c) 2011 The Regents of University of California, Davis Campus.\
                     \n   %s is free software and comes with ABSOLUTELY NO WARRANTY.\
                     \n   Distributed under the MIT License.\
                     \n\
                     \n   Written by: %s\
                     \n\
                     \n",
                    PROGRAM_NAME,
		    VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH,
		    PROGRAM_NAME, AUTHORS);

    exit(status);
}
