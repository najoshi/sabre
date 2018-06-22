#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include "sabre.h"
#include "utils.h"

    //more about getopts http://www.informit.com/articles/article.aspx?p=175771&seqNum=3
    static struct option paired_long_options[] = {
        {"fq1", required_argument, NULL, 'f'},
        {"fq2", required_argument, NULL, 'r'},
        {"barcodes", required_argument, NULL, 'b'},
        {"unassinged1", required_argument, NULL, 'z'},
        {"unassinged2", required_argument, NULL, 'w'},
        {"combine", optional_argument, NULL, 'c'},
        {"umi", optional_argument, NULL, 'u'},
        {"max-mismatch", required_argument, 0, 'm'},
        {"min-umi-len", required_argument, 0, 'l'},
        {"max-5prime-crop", required_argument, 0, 'a'},
        {"stats", required_argument, NULL, 's'},
        {"no-comment", no_argument, 0, 'n'},
        //{"quiet", no_argument, 0, 'z'},
        {GETOPT_HELP_OPTION_DECL},
        {GETOPT_VERSION_OPTION_DECL},
        {NULL, 0, NULL, 0}
    };

void paired_usage (int status) {

    fprintf (stderr, "\n  Usage: %s pe [OPTIONS] -f <fastq_R1> -r <fastq_R2> -b <barcode_file> -u <unassigned_R1> -w <unassigned_R2>\
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
            \n",
            PROGRAM_NAME);

    exit (status);
}

int paired_main (int argc, char *argv[]) {

    //clock_t begin = clock();
    time_t start, end;
    start = time(NULL);

    gzFile pe1=NULL;
    gzFile pe2=NULL;
    kseq_t *fqrec1;
    kseq_t *fqrec2;
    int l1,l2;
    FILE* barfile = NULL;

    gzFile unknownfile1=NULL;
    gzFile unknownfile2=NULL;
    char *unknownfn1=strdup("unassigned_R1.fastq.gz");
    char *unknownfn2=strdup("unassigned_R2.fastq.gz");

    FILE* log_file=NULL;
    int optc;
    extern char *optarg;
    char *fq1=NULL;
    char *fq2=NULL;
    char *barfn=NULL;
    char s_name[MAX_FILENAME_LENGTH];
    barcode_data_paired *curr, *head, *temp;
    char barcode [MAX_BARCODE_LENGTH];
    int num_unknown=0;
    int total=0;
    int mismatch=0;

    int combine = -1;
    int umi = -1;
    int paired = -1;
    int min_umi_len=0;
    int max_5prime_crop=0;
    char *log_fn=NULL;
    int no_comment=-1;

    while (1) {
        int option_index = 0;
        //colon after a flag means should have arguments and no colon means just a flag i.e bool, no args after it
        optc = getopt_long (argc, argv, "dnucf:r:b:z:w:m:s:l:z:a:", paired_long_options, &option_index);

        if (optc == -1) break;

        switch (optc) {
            if (paired_long_options[option_index].flag != 0) break;

            case 'f':
            fq1 = (char*) malloc (strlen (optarg) + 1);
            strcpy (fq1, optarg);
            break;

            case 'r':
            fq2 = (char*) malloc (strlen (optarg) + 1);
            strcpy (fq2, optarg);
            break;

            case 'b':
            barfn = (char*) malloc (strlen (optarg) + 1);
            strcpy (barfn, optarg);
            break;

            case 'z':
            if(unknownfn1) {
                free(unknownfn1);
            }
            unknownfn1 = (char*) malloc (strlen (optarg) + 1);
            strcpy (unknownfn1, optarg);
            break;

            case 'w':
            if(unknownfn2) {
                free(unknownfn2);
            }
            unknownfn2 = (char*) malloc (strlen (optarg) + 1);
            strcpy (unknownfn2, optarg);
            break;

            case 'c':
            combine=1;
            break;

            case 'u':
            umi=1;
            break;

            case 'm':
            mismatch = atoi (optarg);
            break;

            case 's':
            log_fn = (char*) malloc (strlen (optarg) + 1);
            strcpy (log_fn, optarg);
            break;

            case 'l':
            min_umi_len = atoi (optarg);
            break;

            case 'a':
            max_5prime_crop = atoi (optarg);
            break;

            case 'n':
            no_comment = 1;
            break;

            case_GETOPT_HELP_CHAR(paired_usage);
            case_GETOPT_VERSION_CHAR(PROGRAM_NAME, VERSION, AUTHORS);

            case '?':
            paired_usage (EXIT_FAILURE);
            break;

            default:
            paired_usage (EXIT_FAILURE);
            break;
        }
    }

    if (!fq1 || !fq2 || !unknownfn1 || !unknownfn2 || !barfn) {
        paired_usage (EXIT_FAILURE);
    }

    if (!strcmp (fq1, fq2) || !strcmp (fq1, unknownfn1) || !strcmp (fq1, unknownfn2) ||
            !strcmp (fq1, barfn) || !strcmp (fq2, unknownfn1) || !strcmp (fq2, unknownfn2) ||
            !strcmp (fq2, barfn) || !strcmp (unknownfn1, unknownfn2) || !strcmp (unknownfn1, barfn) ||
            !strcmp (unknownfn2, barfn)) {

        fprintf (stderr, "Error: Duplicate input and/or output file names.\n");
        return EXIT_FAILURE;
    }

    pe1 = gzopen (fq1, "r");
    if (!pe1) {
        fprintf (stderr, "Could not open input file 1 '%s'.\n", fq1);
        return EXIT_FAILURE;
    }

    pe2 = gzopen (fq2, "r");
    if (!pe2) {
        fprintf (stderr, "Could not open input file 2 '%s'.\n", fq2);
        return EXIT_FAILURE;
    }

    unknownfile1 = gzopen(unknownfn1, "wb");
    if (!unknownfile1) {
        fprintf (stderr, "Could not open unknown output file 1 '%s'.\n", unknownfn1);
        return EXIT_FAILURE;
    }

    unknownfile2 = gzopen(unknownfn2, "wb");
    if (!unknownfile2) {
        fprintf (stderr, "Could not open unknown output file 2 '%s'.\n", unknownfn2);
        return EXIT_FAILURE;
    }

    barfile = fopen (barfn, "r");
    if (!barfile) {
        fprintf (stderr, "Could not open barcode file '%s'.\n", barfn);
        return EXIT_FAILURE;
    }

    if(fq2) {
        paired = 1;
    }

    fprintf(stderr, "\n\
            \n  Running: %s\
            \n  Command line args:\
            \n      --pe-file1 %s\
            \n      --pe-file2 %s\
            \n      --barcode-file %s\
            \n      --unknown-output1 %s\
            \n      --unknown-output2 %s\
            \n      --combine %d\
            \n      --max-mismatch %d\
            \n      --min-umi-len %d\
            \n      --max-5prime-crop %d\
            \n      --stats %s\
            \n      --no-comment %d\
            \n\
            \n  In Progess...\
            \n", PROGRAM_NAME,\
            fq1, fq2,\
            barfn,\
            unknownfn1, unknownfn2,\
            combine,\
            mismatch, min_umi_len, max_5prime_crop, log_fn, no_comment);

    char *bcout_fn1 = NULL;
    char *bcout_fn2 = NULL;
    /* Creating linked list of barcode data */
    // https://www.hackerearth.com/practice/data-structures/linked-list/singly-linked-list/tutorial/
    // where each node is represents one barcode from the barcode file
    head = NULL;
    while (fscanf (barfile, "%s%s", barcode, s_name) != EOF) {
        curr = (barcode_data_paired*) malloc(sizeof(barcode_data_paired));
        curr->bc = (char*) malloc(strlen(barcode) + 1);
        strcpy(curr->bc, barcode);

        bcout_fn1 = (char *) malloc(MAX_FILENAME_LENGTH*2);
        bcout_fn1[0] = '\0';
        get_bc_fn(&bcout_fn1, s_name, curr->bc, 1);
        //curr->bcfile1 = fopen (_mkdir(bcout_fn1), "w");
        curr->bcfile1 = gzopen(_mkdir(bcout_fn1), "wb");

    if(paired > 0 && combine < 0) {
        bcout_fn2 = (char *) malloc(MAX_FILENAME_LENGTH*2);
        bcout_fn2[0] = '\0';
        get_bc_fn(&bcout_fn2, s_name, curr->bc, 2);
        //curr->bcfile2 = fopen (_mkdir(bcout_fn2), "w");
        curr->bcfile2 = gzopen(_mkdir(bcout_fn2), "wb");
    }

        curr->num_records = 0;
        curr->next = head;
        head = curr;

    }

    free(bcout_fn1);
    free(bcout_fn2);

    fqrec1 = kseq_init (pe1);

    if(paired > 0) {
        fqrec2 = kseq_init (pe2);
    }

    /* Get reads, one at a time */
    while ((l1 = kseq_read (fqrec1)) >= 0) {

        int n_crop = 0;

        char *actl_bc = NULL;
        char *umi_idx = NULL;

        char *fqread1 = NULL;
        char *fqread2 = NULL;

        size_t fq_size = 0;

        fq_size += strlen(fqrec1->seq.s);
        fq_size += (strlen(fqrec1->name.s)*2);
        fq_size += strlen(fqrec1->qual.s);
        fq_size += (strlen(fqrec1->comment.s)*2);
        fq_size += 2;// header signs @ and +
        fq_size += 2;//two colons (:)
        fq_size += 4;//cariage returns
        fq_size += 2;//two spaces
        fq_size += 1000;//test

        if(paired > 0 || combine > 0) {
            l2 = kseq_read (fqrec2);
            if (l2 < 0) {
                fprintf (stderr, "ERROR: R2 file is shorter than R1 file. Disregarding rest of R1 file \n");
                break;
            }
            fq_size += strlen(fqrec2->seq.s);
        }

        /* Find matching barcode */
        curr = head;
        while (curr) {
            n_crop = chk_bc_mtch(curr->bc, fqrec1->seq.s, mismatch, max_5prime_crop);
            if(n_crop >= 0) {
                //found matching barcode
                actl_bc = strndup( (fqrec1->seq.s)+n_crop, strlen(curr->bc) );
                break;
            }
            curr = curr->next;
        }

        /* Write read out into barcode specific file */
        if (curr != NULL) {
            //for now assume barcode and umi are in R1 raed
            if(umi > 0) {
                const char *actl_umi_idx = (fqrec1->seq.s)+strlen(curr->bc)+n_crop;
                umi_idx = strdup(actl_umi_idx);
                umi_idx[strlen(umi_idx)-(max_5prime_crop-n_crop)] = '\0';

                fq_size += strlen(umi_idx);

                if(strlen(umi_idx) < min_umi_len) {
                    break;
                }
            }

            if(combine > 0) {
                fqread1 = (char*) malloc(fq_size);
                fqread1[0] = '\0';

                get_merged_fqread(&fqread1, fqrec1, fqrec2, actl_bc, umi_idx, no_comment, n_crop);
                gzwrite(curr->bcfile1, fqread1, strlen(fqread1));
            }
            else {
                fqread1 = (char*) malloc(fq_size + 1);
                fqread2 = (char*) malloc(fq_size + 1);

                fqread1[0] = '\0';
                fqread2[0] = '\0';

                get_fqread(&fqread1, fqrec1, actl_bc, umi_idx, no_comment, n_crop);
                gzwrite(curr->bcfile1, fqread1, strlen(fqread1));

                if(paired > 0) {
                    get_fqread(&fqread2, fqrec1, actl_bc, umi_idx, no_comment, n_crop);
                    //fprintf(curr->bcfile2, "%s", fqread2);
                    gzwrite(curr->bcfile2, fqread2, strlen(fqread2));
                    curr->num_records += 1;
                }
            }
            curr->num_records += 1;
        }
        else {
            fqread1 = (char*) malloc(fq_size + 1);
            fqread2 = (char*) malloc(fq_size + 1);

            fqread1[0] = '\0';
            fqread2[0] = '\0';

            get_fqread(&fqread1, fqrec1, NULL, NULL, no_comment, 0);
            gzwrite(unknownfile1, fqread1, strlen(fqread1));
            num_unknown += 1;

            if(paired > 0) {
                get_fqread(&fqread2, fqrec2, NULL, NULL, no_comment, 0);
                gzwrite(unknownfile2, fqread2, strlen(fqread2));
                num_unknown += 1;
            }
        }

        total += 2;

        free(fqread1);
        free(fqread2);
        free(actl_bc);
        free(umi_idx);
    }

    if (l1 < 0) {
        l2 = kseq_read (fqrec2);
        if (l2 >= 0) {
            fprintf (stderr, "Error: PE file 1 is shorter than PE file 2. Disregarding rest of PE file 2.\n");
        }
    }

    //if (!quiet) {
    //if (!log_fn) { is this better?
    if (log_fn == NULL) {
        log_file = stdout;
    }
    else {
        log_file = fopen(log_fn, "w");
    }

    fprintf (log_file, "Barcode\tN_records\tN_pairs\tP_pairs\n");
    curr = head;
    int total_pairs = total/2;

    while (curr) {

        int n_pairs = curr->num_records/2;
        float percent_pairs = (float) n_pairs/total_pairs;

        fprintf (log_file,"%s\t%d\t%d\t%.2f\n", curr->bc, curr->num_records, n_pairs, percent_pairs);

        curr = curr->next;
    }

    int unknown_pairs = num_unknown/2;
    float percent_unknown = (float) unknown_pairs/total_pairs;
    float tot_chk = (float) total_pairs/total_pairs;

    fprintf (log_file, "unassigned\t%d\t%d\t%.2f\n", num_unknown, unknown_pairs, percent_unknown);
    fprintf (log_file, "total\t%d\t%d\t%.2f\n", total, total_pairs, tot_chk);

    end = time(NULL);
    fprintf(stderr, "\n All done :) \
                     \n It took %.2f minutes\n",
                     difftime(end, start)/60);

    kseq_destroy (fqrec1);
    kseq_destroy (fqrec2);
    gzclose (pe1);
    gzclose (pe2);
    gzclose (unknownfile1);
    gzclose (unknownfile2);
    fclose (barfile);
    fclose (log_file);

    free (fq1);
    free (fq2);
    free (barfn);
    free (unknownfn1);
    free (unknownfn2);
    free (log_fn);

    curr = head;
    while (curr) {
        gzclose(curr->bcfile1);
        gzclose(curr->bcfile2);

        free (curr->bc);
        temp = curr;
        curr = curr->next;
        free (temp);
    }

    return EXIT_SUCCESS;
}
