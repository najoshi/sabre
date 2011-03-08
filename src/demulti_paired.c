#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include <limits.h>
#include <zlib.h>
#include "sabre.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)


static struct option paired_long_options[] = {
	{"pe-file1", required_argument, 0, 'f'},
	{"pe-file2", required_argument, 0, 'r'},
	{"barcode-file", required_argument, 0, 'b'},
	{"unknown-output1", required_argument, 0, 'u'},
	{"unknown-output2", required_argument, 0, 'w'},
	{"both-barcodes", optional_argument, 0, 'h'},
	{GETOPT_HELP_OPTION_DECL},
	{GETOPT_VERSION_OPTION_DECL},
	{NULL, 0, NULL, 0}
};

void paired_usage (int status) {

	fprintf (stderr, "\nUsage: %s pe -f <paired-end fastq file 1> -r <paired-end fastq file 2> -b <barcode file> -u <unknown barcode output file 1> -w <unknown barcode output file 2>\n\n\
Options:\n\
-f, --pe-file1, Input paired-end fastq file 1 (required, must have same number of records as pe2)\n\
-r, --pe-file2, Input paired-end fastq file 2 (required, must have same number of records as pe1)\n", PROGRAM_NAME);

	fprintf (stderr, "-b, --barcode-file, File with barcodes and output file names, one per line (required)\n\
-u, --unknown-output1, Output paired-end file 1 that contains records with no barcodes found. (required)\n\
-w, --unknown-output2, Output paired-end file 2 that contains records with no barcodes found. (required)\n\
-h, --both-barcodes, Optional flag that indicates that both fastq files have barcodes.\n\
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
	int debug=0;
	int optc;
	extern char *optarg;
	char *infn1=NULL;
	char *infn2=NULL;
	char *barfn=NULL;
	char *unknownfn1=NULL;
	char *unknownfn2=NULL;
	int both_have_barcodes=0;
	barcode_data_paired *curr, *head;
	char barcode [MAX_BARCODE_LENGTH];
	char baroutfn1 [MAX_FILENAME_LENGTH];
	char baroutfn2 [MAX_FILENAME_LENGTH];


	while (1) {
		int option_index = 0;
		optc = getopt_long (argc, argv, "dhf:r:b:u:w:", paired_long_options, &option_index);

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

			case 'h':
				both_have_barcodes=1;
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


	head = NULL;
	while (fscanf (barfile, "%s%s%s", barcode, baroutfn1, baroutfn2) != EOF) {
		curr = (barcode_data_paired*) malloc (sizeof (barcode_data_paired));
		curr->bc = (char*) malloc (strlen(barcode) + 1);
		strcpy (curr->bc, barcode);

		curr->bcfile1 = fopen (baroutfn1, "w");
		curr->bcfile2 = fopen (baroutfn2, "w");

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


		curr = head;
		while (curr) {
			if (strncmp (curr->bc, fqrec1->seq.s, strlen (curr->bc)) == 0) {
				break;
			}

			curr = curr->next;
		}


		if (curr != NULL) {
			fprintf (curr->bcfile1, "@%s", fqrec1->name.s);
			if (fqrec1->comment.l) fprintf (curr->bcfile1, " %s\n", fqrec1->comment.s);
			else fprintf (curr->bcfile1, "\n");

			fprintf (curr->bcfile1, "%s\n", (fqrec1->seq.s)+strlen(curr->bc));

			fprintf (curr->bcfile1, "+%s", fqrec1->name.s);
			if (fqrec1->comment.l) fprintf (curr->bcfile1, " %s\n", fqrec1->comment.s);
			else fprintf (curr->bcfile1, "\n");

			fprintf (curr->bcfile1, "%s\n", (fqrec1->qual.s)+strlen(curr->bc));


			fprintf (curr->bcfile2, "@%s", fqrec2->name.s);
			if (fqrec2->comment.l) fprintf (curr->bcfile2, " %s\n", fqrec2->comment.s);
			else fprintf (curr->bcfile2, "\n");

			if (!both_have_barcodes) fprintf (curr->bcfile2, "%s\n", fqrec2->seq.s);
			else fprintf (curr->bcfile2, "%s\n", (fqrec2->seq.s)+strlen(curr->bc));

			fprintf (curr->bcfile2, "+%s", fqrec2->name.s);
			if (fqrec2->comment.l) fprintf (curr->bcfile2, " %s\n", fqrec2->comment.s);
			else fprintf (curr->bcfile2, "\n");

			if (!both_have_barcodes) fprintf (curr->bcfile2, "%s\n", fqrec2->qual.s);
			else fprintf (curr->bcfile2, "%s\n", (fqrec2->qual.s)+strlen(curr->bc));
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
		}
	}

	if (l1 < 0) {
		l2 = kseq_read (fqrec2);
		if (l2 >= 0) {
			fprintf (stderr, "Error: PE file 1 is shorter than PE file 2. Disregarding rest of PE file 2.\n");
		}
	}

	/*fprintf (stderr, "\nFastQ paired records kept: %d (%d pairs)\nFastQ single records kept: %d (from PE1: %d, from PE2: %d)\nFastQ paired records discarded: %d (%d pairs)\nFastQ single records discarded: %d (from PE1: %d, from PE2: %d)\n\n", kept_p, (kept_p/2), (kept_s1+kept_s2), kept_s1, kept_s2, discard_p, (discard_p/2), (discard_s1+discard_s2), discard_s1, discard_s2);*/

	kseq_destroy (fqrec1);
	kseq_destroy (fqrec2);
	gzclose (pe1);
	gzclose (pe2);
	fclose (unknownfile1);
	fclose (unknownfile2);
	fclose (barfile);

	curr = head;
	while (curr) {
		fclose (curr->bcfile1);
		fclose (curr->bcfile2);
		curr = curr->next;
	}

	return EXIT_SUCCESS;
}
