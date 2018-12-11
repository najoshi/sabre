
void* sanity_chk(void *arg)
    if (!fq1 || !fq2 || !unknownfn1 || !unknownfn2 || !barfn) {
        paired_usage (EXIT_FAILURE);
    }

    if (!strcmp (fq1, fq2) || !strcmp (fq1, unknownfn1) || !strcmp (fq1, unknownfn2) ||
            !strcmp (fq1, barfn) || !strcmp (fq2, unknownfn1) || !strcmp (fq2, unknownfn2) ||
            !strcmp (fq2, barfn) || !strcmp (unknownfn1, unknownfn2) || !strcmp (unknownfn1, barfn) ||
            !strcmp (unknownfn2, barfn)) {

        fprintf (stderr, "ERROR: Duplicate input and/or output file names.\n");
        return EXIT_FAILURE;
    }

    pe1 = gzopen (fq1, "r");
    if (!pe1) {
        fprintf (stderr, "ERROR: Could not open input file 1 '%s'.\n", fq1);
        return EXIT_FAILURE;
    }

    pe2 = gzopen (fq2, "r");
    if (!pe2) {
        fprintf (stderr, "ERROR: Could not open input file 2 '%s'.\n", fq2);
        return EXIT_FAILURE;
    }

    unknownfile1 = gzopen(unknownfn1, "wb");
    if (!unknownfile1) {
        fprintf (stderr, "ERROR: Could not open unknown output file 1 '%s'.\n", unknownfn1);
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

    if(threads < 0) {
        fprintf(stderr, "WARNING: Negative number of threads detected %d, setting threads to 1\n", threads);
        threads = 1;
    }
}
