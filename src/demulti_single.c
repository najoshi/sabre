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


static struct option single_long_options[] = {
  {"fastq-file", required_argument, 0, 'f'},
  {"barcode-file", required_argument, 0, 'b'},
  {"unknown-output", required_argument, 0, 'u'},
  {"max-mismatch", optional_argument, 0, 'm'},
  {"quiet", optional_argument, 0, 'z'},
  {GETOPT_HELP_OPTION_DECL},
  {GETOPT_VERSION_OPTION_DECL},
  {NULL, 0, NULL, 0}
};

void single_usage (int status) {
  
  fprintf (stderr, "\nUsage: %s se -f <fastq sequence file> -b <barcode file> -u <unknown barcode output file>\n\
\n\
Options:\n\
-f, --fastq-file, Input fastq file (required)\n", PROGRAM_NAME);

  fprintf (stderr, "-b, --barcode-file, File with barcode and output file name per line (required)\n\
-u, --unknown-output, Output file that contains records with no barcodes found. (required)\n\
-m <n>, --max-mismatch <n>, Optional argument that is the maximum number of mismatches allowed in a barcode. Default 0.\n\
--quiet, don't output matching info\n\
--help, display this help and exit\n\
--version, output version information and exit\n\n");

  exit (status);
}

int single_main (int argc, char *argv[]) {

	gzFile se=NULL;
	kseq_t *fqrec;
	int debug=0;
	int optc;
	extern char *optarg;
	FILE* barfile = NULL;
	FILE* unknownfile=NULL;
	char *barfn=NULL;
	char *infn=NULL;
	char *unknownfn=NULL;
	barcode_data *curr, *head, *temp;
	char barcode [MAX_BARCODE_LENGTH];
	char baroutfn [MAX_FILENAME_LENGTH];
	int num_unknown=0;
	int total=0;
	int mismatch=0;
	int quiet=0;

	while (1) {
		int option_index = 0;
		optc = getopt_long (argc, argv, "df:b:u:m:z", single_long_options, &option_index);

		if (optc == -1) break;

		switch (optc) {
			if (single_long_options[option_index].flag != 0) break;

			case 'f':
				infn = (char*) malloc (strlen (optarg) + 1);
				strcpy (infn, optarg);
				break;

			case 'b':
				barfn = (char*) malloc (strlen (optarg) + 1);
				strcpy (barfn, optarg);
				break;

			case 'u':
				unknownfn = (char*) malloc (strlen (optarg) + 1);
				strcpy (unknownfn, optarg);
				break;

			case 'm':
				mismatch = atoi (optarg);
				break;

			case 'z':
				quiet=1;
				break;

			case 'd':
				debug = 1;
				break;

			case_GETOPT_HELP_CHAR(single_usage)
			case_GETOPT_VERSION_CHAR(PROGRAM_NAME, VERSION, AUTHORS);

			case '?':
				single_usage (EXIT_FAILURE);
				break;

			default:
				single_usage (EXIT_FAILURE);
				break;
		}
	}


	if (!infn || !barfn) {
		single_usage (EXIT_FAILURE);
	}

	if (!strcmp (infn, barfn)) {
		fprintf (stderr, "Error: Input file is same as barcode file.\n");
		return EXIT_FAILURE;
	}

	se = gzopen (infn, "r");
	if (!se) {
		fprintf (stderr, "Could not open input file '%s'.\n", infn);
		return EXIT_FAILURE;
	}

	barfile = fopen (barfn, "r");
	if (!barfile) {
		fprintf (stderr, "Could not open barcode file '%s'.\n", barfn);
		return EXIT_FAILURE;
	}

	unknownfile = fopen (unknownfn, "w");
	if (!unknownfile) {
		fprintf (stderr, "Could not open unknown output file '%s'.\n", unknownfn);
		return EXIT_FAILURE;
	}


	/* Set up a linked list of the barcode data */
	head = NULL;
	while (fscanf (barfile, "%s%s", barcode, baroutfn) != EOF) {
		curr = (barcode_data*) malloc (sizeof (barcode_data));
		curr->bc = (char*) malloc (strlen(barcode) + 1);
		strcpy (curr->bc, barcode);

		curr->bcfile = fopen (baroutfn, "w");
		curr->num_records = 0;

		curr->next = head;
		head = curr;
	}

	
	fqrec = kseq_init (se);

	while (kseq_read (fqrec) >= 0) {

		/* Go through linked list of barcode data and compare each barcode */
		/* with the sequence until a match is found or no match is found for any */
		curr = head;
		while (curr) {
			if (strncmp_with_mismatch (curr->bc, fqrec->seq.s, strlen (curr->bc), mismatch) == 0) {
				break;
			}

			curr = curr->next;
		}


		/* If barcode data is found, output to demultiplexed file, else output to unknown file */
		if (curr != NULL) {
			fprintf (curr->bcfile, "@%s:%s", fqrec->name.s, curr->bc);
			if (fqrec->comment.l) fprintf (curr->bcfile, " %s\n", fqrec->comment.s);
			else fprintf (curr->bcfile, "\n");

			fprintf (curr->bcfile, "%s\n", (fqrec->seq.s)+strlen(curr->bc));

			fprintf (curr->bcfile, "+%s", fqrec->name.s);
			if (fqrec->comment.l) fprintf (curr->bcfile, " %s\n", fqrec->comment.s);
			else fprintf (curr->bcfile, "\n");

			fprintf (curr->bcfile, "%s\n", (fqrec->qual.s)+strlen(curr->bc));

			curr->num_records++;
		}

		else {
			fprintf (unknownfile, "@%s", fqrec->name.s);
			if (fqrec->comment.l) fprintf (unknownfile, " %s\n", fqrec->comment.s);
			else fprintf (unknownfile, "\n");

			fprintf (unknownfile, "%s\n", fqrec->seq.s);

			fprintf (unknownfile, "+%s", fqrec->name.s);
			if (fqrec->comment.l) fprintf (unknownfile, " %s\n", fqrec->comment.s);
			else fprintf (unknownfile, "\n");

			fprintf (unknownfile, "%s\n", fqrec->qual.s);

			num_unknown++;
		}

		total++;
	}


	if (!quiet) {
		fprintf (stdout, "\nTotal FastQ records: %d\n\n", total);
		curr = head;
		while (curr) {
			fprintf (stdout, "FastQ records for barcode %s: %d\n", curr->bc, curr->num_records);
			curr = curr->next;
		}
		fprintf (stdout, "\nFastQ records with no barcode match: %d\n", num_unknown);
		fprintf (stdout, "\nNumber of mismatches allowed: %d\n\n", mismatch);
	}

	kseq_destroy (fqrec);
	gzclose (se);
	fclose (barfile);
	fclose (unknownfile);

	free (infn);
	free (unknownfn);
	free (barfn);

	curr = head;
	while (curr) {
		fclose (curr->bcfile);
		free (curr->bc);
		temp = curr;
		curr = curr->next;
		free (temp);
	}

	return 0;
}
