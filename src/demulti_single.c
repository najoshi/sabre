#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include "sabre.h"
#include "kseq.h"

__KS_GETC(gzread, BUFFER_SIZE)
__KS_GETUNTIL(gzread, BUFFER_SIZE)
__KSEQ_READ

static struct option single_long_options[] = {
  {"fastq-file", required_argument, 0, 'f'},
  {"barcode-file", required_argument, 0, 'b'},
  {GETOPT_HELP_OPTION_DECL},
  {GETOPT_VERSION_OPTION_DECL},
  {NULL, 0, NULL, 0}
};

void single_usage (int status) {
  
  fprintf (stderr, "\nUsage: %s se -f <fastq/a sequence file> -b <barcode file> -u <unknown barcode output file>\n\
\n\
Options:\n\
-f, --fastq-file, Input fastq/a file (required)\n", PROGRAM_NAME);

  fprintf (stderr, "-b, --barcode-file, File with barcodes and output file names, one per line (required)\n\
-u, --unknown-file, Output file that contains records with no barcodes found. (required)\n\
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
	barcode_data *curr, *head;
	char barcode [MAX_BARCODE_LENGTH];
	char baroutfn [MAX_FILENAME_LENGTH];

	while (1) {
		int option_index = 0;
		optc = getopt_long (argc, argv, "df:b:u:",single_long_options, &option_index);

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


	head = NULL;
	while (fscanf (barfile, "%s%s", barcode, baroutfn) != EOF) {
		curr = (barcode_data*) malloc (sizeof (barcode_data));
		curr->bc = (char*) malloc (strlen(barcode) + 1);
		strcpy (curr->bc, barcode);

		curr->bcfile = fopen (baroutfn, "w");

		curr->next = head;
		head = curr;
	}

	
	fqrec = kseq_init (se);

	while (kseq_read (fqrec) >= 0) {

		curr = head;
		while (curr) {
			if (strncmp (curr->bc, fqrec->seq.s, strlen (curr->bc)) == 0) {
				break;
			}

			curr = curr->next;
		}

		if (curr != NULL) {
			fprintf (curr->bcfile, "@%s", fqrec->name.s);
			if (fqrec->comment.l) fprintf (curr->bcfile, " %s\n", fqrec->comment.s);
			else fprintf (curr->bcfile, "\n");

			fprintf (curr->bcfile, "%s\n", (fqrec->seq.s)+strlen(curr->bc));

			fprintf (curr->bcfile, "+%s", fqrec->name.s);
			if (fqrec->comment.l) fprintf (curr->bcfile, " %s\n", fqrec->comment.s);
			else fprintf (curr->bcfile, "\n");

			fprintf (curr->bcfile, "%s\n", (fqrec->qual.s)+strlen(curr->bc));
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
		}
	}

	/*fprintf (stderr, "\nFastQ records kept: %d\nFastQ records discarded: %d\n\n", kept, discard);*/

	kseq_destroy (fqrec);
	gzclose (se);

	curr = head;
	while (curr) {
		fclose (curr->bcfile);
		curr = curr->next;
	}

	return 0;
}
