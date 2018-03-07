#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <limits.h>
#include <zlib.h>
#include "sabre.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

//more about getopts http://www.informit.com/articles/article.aspx?p=175771&seqNum=3
static struct option paired_long_options[] = {
	{"pe-file1", required_argument, NULL, 'f'},
	{"pe-file2", required_argument, NULL, 'r'},
	{"barcode-file", required_argument, NULL, 'b'},
	{"unknown-output1", required_argument, NULL, 'u'},
	{"unknown-output2", required_argument, NULL, 'w'},
	{"both-barcodes", optional_argument, NULL, 'c'},
	{"max-mismatch", required_argument, 0, 'm'},
	{"min-umi-len", required_argument, 0, 'l'},
	{"stats", required_argument, NULL, 's'},
	{"no-comment", no_argument, 0, 'n'},
	//{"quiet", no_argument, 0, 'z'},
	{GETOPT_HELP_OPTION_DECL},
	{GETOPT_VERSION_OPTION_DECL},
	{NULL, 0, NULL, 0}
};

void paired_usage (int status) {

	fprintf (stderr, "\nUsage: %s pe -f <paired-end fastq file 1> -r <paired-end fastq file 2> -b <barcode file> -u <unknown barcode output file 1> -w <unknown barcode output file 2>\n\n\
Options:\n\
-f, --pe-file1, Input paired-end fastq file 1 (required, must have same number of records as pe2)\n\
-r, --pe-file2, Input paired-end fastq file 2 (required, must have same number of records as pe1)\n\
-b, --barcode-file, File with barcode and two output file names per line (required)\n", PROGRAM_NAME);

fprintf (stderr, "-u, --unknown-output1, Output paired-end file 1 that contains records with no barcodes found. (required)\n\
-w, --unknown-output2, Output paired-end file 2 that contains records with no barcodes found. (required)\n\
-c, --both-barcodes, Optional flag that indicates that both fastq files have barcodes.\n\
-m <n>, --max-mismatch <n>, Optional argument that is the maximum number of mismatches allowed in a barcode. Default 0.\n\
-l <n>, --min-umi-len <n>, Optional argument that is the minimum UMI length to keep. Default [0].\n\
-n, --no-comment, Optional argument to drop extra comments from FASTQ header. Default [NULL].\n\
-s <FILENAME>, --stats <FILENAME>, Optional argument to write logs into a file instead of STDOUT. Default [STDOUT].\n");

fprintf (stderr, "--quiet, don't print barcode matching info\n\
--help, display this help and exit\n\
--version, output version information and exit\n\n");

	exit (status);
}

int paired_main (int argc, char *argv[]) {

	gzFile pe1=NULL;
	gzFile pe2=NULL;
	kseq_t *fqrec1;
	kseq_t *fqrec2;
	int l1,l2;
	FILE* barfile = NULL;
	FILE* unknownfile1=NULL;
	FILE* unknownfile2=NULL;
	FILE* log_file=NULL;
	int debug=0;
	int optc;
	extern char *optarg;
	char *infn1=NULL;
	char *infn2=NULL;
	char *barfn=NULL;
	char *unknownfn1=NULL;
	char *unknownfn2=NULL;
	int both_have_barcodes=0;
	barcode_data_paired *curr, *head, *temp;
	char barcode [MAX_BARCODE_LENGTH];
	char baroutfn1 [MAX_FILENAME_LENGTH];
	char baroutfn2 [MAX_FILENAME_LENGTH];
	int num_unknown=0;
	int total=0;
	int mismatch=0;
	//int quiet=0;

	int min_umi_len=0;
	char *log_fn=NULL;
	int no_comment=-1;


	while (1) {
		int option_index = 0;
		optc = getopt_long (argc, argv, "dcf:r:b:u:w:m:s:l:n:z", paired_long_options, &option_index);

		if (optc == -1) break;

		switch (optc) {
			if (paired_long_options[option_index].flag != 0) break;

			case 'f':
				infn1 = (char*) malloc (strlen (optarg) + 1);
				strcpy (infn1, optarg);
				break;

			case 'r':
				infn2 = (char*) malloc (strlen (optarg) + 1);
				strcpy (infn2, optarg);
				break;

			case 'b':
				barfn = (char*) malloc (strlen (optarg) + 1);
				strcpy (barfn, optarg);
				break;

			case 'u':
				unknownfn1 = (char*) malloc (strlen (optarg) + 1);
				strcpy (unknownfn1, optarg);
				break;

			case 'w':
				unknownfn2 = (char*) malloc (strlen (optarg) + 1);
				strcpy (unknownfn2, optarg);
				break;

			case 'c':
				both_have_barcodes=1;
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

			case 'n':
				no_comment = 1;
				break;

			case 'z':
				//quiet=1;
				break;

			case 'd':
				debug = 1;
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

	if (!infn1 || !infn2 || !unknownfn1 || !unknownfn2 || !barfn) {
		paired_usage (EXIT_FAILURE);
	}

	if (!strcmp (infn1, infn2) || !strcmp (infn1, unknownfn1) || !strcmp (infn1, unknownfn2) ||
		!strcmp (infn1, barfn) || !strcmp (infn2, unknownfn1) || !strcmp (infn2, unknownfn2) ||
		!strcmp (infn2, barfn) || !strcmp (unknownfn1, unknownfn2) || !strcmp (unknownfn1, barfn) ||
		!strcmp (unknownfn2, barfn)) {

		fprintf (stderr, "Error: Duplicate input and/or output file names.\n");
		return EXIT_FAILURE;
	}

	pe1 = gzopen (infn1, "r");
	if (!pe1) {
		fprintf (stderr, "Could not open input file 1 '%s'.\n", infn1);
		return EXIT_FAILURE;
	}

	pe2 = gzopen (infn2, "r");
	if (!pe2) {
		fprintf (stderr, "Could not open input file 2 '%s'.\n", infn2);
		return EXIT_FAILURE;
	}

	unknownfile1 = fopen (unknownfn1, "w");
	if (!unknownfile1) {
		fprintf (stderr, "Could not open unknown output file 1 '%s'.\n", unknownfn1);
		return EXIT_FAILURE;
	}

	unknownfile2 = fopen (unknownfn2, "w");
	if (!unknownfile2) {
		fprintf (stderr, "Could not open unknown output file 2 '%s'.\n", unknownfn2);
		return EXIT_FAILURE;
	}

	barfile = fopen (barfn, "r");
	if (!barfile) {
		fprintf (stderr, "Could not open barcode file '%s'.\n", barfn);
		return EXIT_FAILURE;
	}

        fprintf(stderr, "\n\
                       \n  Running: %s\
                       \n  Command line args:\
                       \n      --pe-file1 %s\
                       \n      --pe-file2 %s\
                       \n      --barcode-file %s\
                       \n      --unknown-output1 %s\
                       \n      --unknown-output2 %s\
                       \n      --both-barcodes %d\
                       \n      --max-mismatch %d\
                       \n      --min-umi-len %d\
                       \n      --stats %s\
                       \n      --no-comment %d\
                       \n\
                       \n  In Progess...\
                       \n", PROGRAM_NAME,\
                            infn1, infn2,\
                            barfn,\
                            unknownfn1, unknownfn2,\
                            both_have_barcodes,\
                            mismatch, min_umi_len, log_fn, no_comment);


	/* Creating linked list of barcode data */
	head = NULL;
	while (fscanf (barfile, "%s%s%s", barcode, baroutfn1, baroutfn2) != EOF) {
		curr = (barcode_data_paired*) malloc (sizeof (barcode_data_paired));
		curr->bc = (char*) malloc (strlen(barcode) + 1);
		strcpy (curr->bc, barcode);

		curr->bcfile1 = fopen (_mkdir(baroutfn1), "w");
		curr->bcfile2 = fopen (_mkdir(baroutfn2), "w");
		curr->num_records = 0;

		curr->next = head;
		head = curr;
	}


	fqrec1 = kseq_init (pe1);
	fqrec2 = kseq_init (pe2);

	while ((l1 = kseq_read (fqrec1)) >= 0) {

		l2 = kseq_read (fqrec2);
		if (l2 < 0) {
			fprintf (stderr, "Error: PE file 2 is shorter than PE file 1. Disregarding rest of PE file 1.\n");
			break;
		}


		/* Go through all barcode data and check if any match to beginning of read */
		/* If it does then put read in that barcode's file, otherwise put in unknown file */
		curr = head;
		while (curr) {
			if (strncmp_with_mismatch (curr->bc, fqrec1->seq.s, strlen (curr->bc), mismatch) == 0) {
				break;
			}

			curr = curr->next;
		}


		if (curr != NULL) {
			// if UMI is shorter then 10, discard the reads
			if(strlen((fqrec1->seq.s)+strlen(curr->bc)) >= min_umi_len) {
			    //@READNAME:BACRCODE:UMI
			    fprintf (curr->bcfile1, "@%s:%s:%s", fqrec1->name.s, curr->bc, (fqrec1->seq.s)+strlen(curr->bc));
			    if (fqrec1->comment.l && no_comment == -1) fprintf (curr->bcfile1, " %s\n", fqrec1->comment.s);
			    else fprintf (curr->bcfile1, "\n");

			    //fprintf (curr->bcfile1, "%s\n", (fqrec1->seq.s)+strlen(curr->bc));
			    //This tmp hack knowning that data is single end, and R2 is simply a string of BARCODE+UMI
			    fprintf (curr->bcfile1, "N\n");

			    fprintf (curr->bcfile1, "+%s", fqrec1->name.s);
			    if (fqrec1->comment.l) fprintf (curr->bcfile1, " %s\n", fqrec1->comment.s);
			    else fprintf (curr->bcfile1, "\n");

			    fprintf (curr->bcfile1, "%s\n", (fqrec1->qual.s)+strlen(curr->bc));


			    //fprintf (curr->bcfile2, "@%s:%s", fqrec2->name.s, curr->bc);
			    //@READNAME:BACRCODE:UMI
			    fprintf (curr->bcfile2, "@%s:%s:%s", fqrec2->name.s, curr->bc, (fqrec1->seq.s)+strlen(curr->bc));
			    if (fqrec2->comment.l && no_comment == -1) fprintf (curr->bcfile2, " %s\n", fqrec2->comment.s);
			    else fprintf (curr->bcfile2, "\n");

			    if (!both_have_barcodes) fprintf (curr->bcfile2, "%s\n", fqrec2->seq.s);
			    else fprintf (curr->bcfile2, "%s\n", (fqrec2->seq.s)+strlen(curr->bc));

			    fprintf (curr->bcfile2, "+%s", fqrec2->name.s);
			    if (fqrec2->comment.l) fprintf (curr->bcfile2, " %s\n", fqrec2->comment.s);
			    else fprintf (curr->bcfile2, "\n");

			    if (!both_have_barcodes) fprintf (curr->bcfile2, "%s\n", fqrec2->qual.s);
			    else fprintf (curr->bcfile2, "%s\n", (fqrec2->qual.s)+strlen(curr->bc));

			    curr->num_records += 2;
			}
		}

		else {
			fprintf (unknownfile1, "@%s", fqrec1->name.s);
			if (fqrec1->comment.l) fprintf (unknownfile1, " %s\n", fqrec1->comment.s);
			else fprintf (unknownfile1, "\n");

			fprintf (unknownfile1, "%s\n", fqrec1->seq.s);

			fprintf (unknownfile1, "+%s", fqrec1->name.s);
			if (fqrec1->comment.l) fprintf (unknownfile1, " %s\n", fqrec1->comment.s);
			else fprintf (unknownfile1, "\n");

			fprintf (unknownfile1, "%s\n", fqrec1->qual.s);


			fprintf (unknownfile2, "@%s", fqrec2->name.s);
			if (fqrec2->comment.l) fprintf (unknownfile2, " %s\n", fqrec2->comment.s);
			else fprintf (unknownfile2, "\n");

			fprintf (unknownfile2, "%s\n", fqrec2->seq.s);

			fprintf (unknownfile2, "+%s", fqrec2->name.s);
			if (fqrec2->comment.l) fprintf (unknownfile2, " %s\n", fqrec2->comment.s);
			else fprintf (unknownfile2, "\n");

			fprintf (unknownfile2, "%s\n", fqrec2->qual.s);

			num_unknown += 2;
		}

		total += 2;
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
        
        fprintf (log_file, "\nTotal FastQ records: %d (%d pairs)\n\n", total, total/2);
        curr = head;
        while (curr) {
        	fprintf (log_file, "FastQ records for barcode %s: %d (%d pairs)\n", curr->bc, curr->num_records, curr->num_records/2);
        	curr = curr->next;
        }
        fprintf (log_file, "\nFastQ records with no barcode match: %d (%d pairs)\n", num_unknown, num_unknown/2);
        fprintf (log_file, "\nNumber of mismatches allowed: %d\n\n", mismatch);

        fprintf (stderr, "\n  All done :)!");

	kseq_destroy (fqrec1);
	kseq_destroy (fqrec2);
	gzclose (pe1);
	gzclose (pe2);
	fclose (unknownfile1);
	fclose (unknownfile2);
	fclose (barfile);
	fclose (log_file);

	free (infn1);
	free (infn2);
	free (barfn);
	free (unknownfn1);
	free (unknownfn2);
	free (log_fn);

	curr = head;
	while (curr) {
		fclose (curr->bcfile1);
		fclose (curr->bcfile2);
		free (curr->bc);
		temp = curr;
		curr = curr->next;
		free (temp);
	}

	return EXIT_SUCCESS;
}
