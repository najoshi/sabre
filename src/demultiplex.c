
/*
 * set softtab=4
 * set shiftwidth=4
 * set expandtab
 *
 */

/* 
 * sabre FASTQ files demultiplexing
 * demultiplex.c: FASTQ demultiplexing
 * 
 */

int demulti() {

    barcode_data_paired *curr, *head, *temp;
    char barcode [MAX_BARCODE_LENGTH];
    char s_name [MAX_FILENAME_LENGTH];

    /* Creating linked list of barcode data */
    // https://www.hackerearth.com/practice/data-structures/linked-list/singly-linked-list/tutorial/
    // where each node is represents one barcode from the barcode file
    // number of nodes should equal to number of barcodes (lines) in the file
    head = NULL;
    while (fscanf (barfile, "%s%s", barcode, s_name) != EOF) {
        curr = (barcode_data_paired*) malloc (sizeof (barcode_data_paired));
        curr->bc = (char*) malloc (strlen(barcode) + 1);
        strcpy(curr->bc, barcode);

        char *bcout_fn1 = get_bc_fn(barcode, s_name, 1);
        char *bcout_fn2 = get_bc_fn(barcode, s_name, 1);

        curr->bcfile1 = fopen (_mkdir(bcout_fn1), "w");
        curr->bcfile2 = fopen (_mkdir(bcout_fn2), "w");
        curr->num_records = 0;

        curr->next = head;
        head = curr;
    }

    fqrec1 = kseq_init (pe1);

    if(paired > 0) {
        fqrec2 = kseq_init (pe2);
    }

    /* Get reads, one at a time */
    while((l1 = kseq_read (fqrec1)) >= 0) {

	int n_crop_fq1;
	int n_crop_fq2;
	char *actl_bc_fq1 = [MAX_BARCODE_LENGTH];
	char *actl_bc_fq2 = [MAX_BARCODE_LENGTH];

        if(paired > 0) {
            l2 = kseq_read (fqrec2);
            if (l2 < 0) {
                fprintf (stderr, "ERROR: R2 file is shorter than R1 file. Disregarding rest of R1 file \n");
                break;
            }
        }

        /* Find matching barcode */
        curr = head;
        while (curr) {

            n_crop_fq1 = chk_bc_mtch(curr->bc, fqrec1->seq.s, mismatch, max_5prime_crop);
            if (n_crop_fq1 >= 0) {
                //found matching barcode
                break;
            }

	    if(paired > 0) {

		n_crop_fq2 = chk_bc_mtch(curr->bc, fqrec2->seq.s, mismatch, max_5prime_crop);

            	if (n_crop_fq2 < 0) {
		    // it is ok not to have matching barcode..
		    fprintf (stderr, "ERROR: R2 didn't have a matching barcode. \n");
            	}

            	if (n_crop_fq1 != n_crop_fq2) {
		    // this will go heand in heand with previous check
		    // but can be stand along thing as well, when one read has an overhand 
		    // and the other doesn't, shouldn't be the case though (I think)
		    fprintf (stderr, "ERROR: Number of cropped bases doesn't match between R1 and R2\n");
		}

		actl_bc_fq1 = strlen(fqrec1->seq.s)+strlen(curr->bc)+n_crop_fq1;
		actl_bc_fq2 = strlen(fqrec2->seq.s)+strlen(curr->bc)+n_crop_fq2;

                // didn't match if != zero
		if(strcmp(actl_bc_fq1, actl_bc_fq2) != 0) {
		    fprintf (stderr, "ERROR: Actual R1 and R2 barcodes didn't match, %s and %s. This is strange.. \n",
			              actl_bc_fq1,
				      actl_bc_fq2);
		}
		else {
		    //write read out to a matching barcode
		    break;
		}
	    }

            curr = curr->next;
        }

        /* Write read out into barcode specific file */
        if(curr != NULL) {
            // if UMI is shorter then 10, discard the reads
            //if(strlen((fqrec1->seq.s)+strlen(curr->bc)) >= min_umi_len) {

            fqrec1->name.s
            fqrec1->seq.s
            fqrec1->comment.l
            fqrec1->qual.s
            curr->bc

            char *trimed_fq1 = (fqrec1->seq.s)+strlen(curr->bc);

	    //}
    }
}
