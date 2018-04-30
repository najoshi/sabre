
/*
 * Not valid C code, just some snipets for later use
 *
 * TODO build on get_fqread function in utils. 
 * I think that function should also take char *umi
 * if umi == NULL, then none of 3? modes are true, return normal
 * if umi != NULL  then this must be one of 3? modes, merge R1 and R2 and return R1 only
 */

typedef struct listel_p {
	char* bc;
	int num_records;
	//FILE* bcfile1;
	//FILE* bcfile2;
	gzFile bcfile1;
	gzFile bcfile2;
	struct listel_p *next;
} barcode_data_paired;

    /* Creating linked list of barcode data */
    // https://www.hackerearth.com/practice/data-structures/linked-list/singly-linked-list/tutorial/
    // where each node is represents one barcode from the barcode file
    // number of nodes should equal to number of barcodes (lines) in the file
    head = NULL;
    char barcode [MAX_BARCODE_LENGTH];
    char s_name [MAX_SNAME_LENGTH];
    //while (fscanf (barfile, "%s%s%s", barcode, baroutfn1, baroutfn2) != EOF) {
    while (fscanf (barfile, "%s\t%s", barcode, s_name) != EOF) {
        char bcout_prefix [MAX_BARCODE_LENGTH+MAX_SNAME_LENGTH];
        char bcout_fn1 [MAX_FILENAME_LENGTH];
        char bcout_fn2 [MAX_FILENAME_LENGTH];

        curr = (barcode_data_paired*) malloc (sizeof (barcode_data_paired));
        curr->bc = (char*) malloc (strlen(barcode) + 1);
        strcpy (curr->bc, barcode);

        if(strlen(s_name) > MAX_FILENAME_LENGTH) {
            fprintf (stderr, "ERROR: Too many characters in your sample name; %s:%d \n", s_name, strlen(s_name));
        }
        //TODO make this into a function call later on.
        //want a function in utils.c get_bc_fn(s_name, barcode, 1|2) to return 
        //a string = bcout_fn to... maybe this isn't worth a function call..
        strcat(bcout_prefix, s_name);
        strcat(bcout_prefix, "_");
        strcat(bcout_prefix, barcode);

        strcpy(bcout_fn1, bcout_prefix);
        strcat(bcout_fn1, "_R1.fastq.gz");

        strcpy(bcout_fn2, bcout_prefix);
        strcat(bcout_fn2, "_R2.fastq.gz");

        curr->bcfile1 = gzopen(_mkdir(bcout_fn1), "wb");
        curr->bcfile2 = gzopen(_mkdir(bcout_fn2), "wb");
        curr->num_records = 0;

        curr->next = head;
        head = curr;
    }

    while (curr) {
        gzclose(curr->bcfile1);
        gzclose(curr->bcfile2);
        free (curr->bc);
        temp = curr;
        curr = curr->next;
        free (temp);
    }

