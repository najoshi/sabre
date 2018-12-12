#include "sabre.h"
#include "utils.h"
#include "usage.h"
#include "demultiplex.h"

int main(int argc, char *argv[]) {

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
        {"threads", optional_argument, 0, 't'},
        {"version", optional_argument, NULL, 'v'},
        {"help", optional_argument, NULL, 'h'},
        {"story", optional_argument, NULL, 'o'},
        //{"quiet", no_argument, 0, 'z'},
        {NULL, 0, NULL, 0}
    };

    // this is on the stack so no need to malloc 
    // stack gets cleaned when function exits,
    // because of that no need to free either
    param_t params;
    set_default_params(&params);

    metrics_t metrics;

    //clock_t begin = clock();
    time_t start, end;
    start = time(NULL);

    FILE* barfile = NULL;

    char *unassigned1_fn=strdup("unassigned_R1.fq.gz");
    char *unassigned2_fn=strdup("unassigned_R2.fq.gz");
    char *umis_2_short_fn=strdup("umis_too_short.txt");

    FILE* log_file=NULL;
    int optc;
    extern char *optarg;

    char *fq1_fn=NULL;
    char *fq2_fn=NULL;
    char *log_fn=NULL;

    char *barfn=NULL;
    char s_name[MAX_FILENAME_LENGTH];
    barcode_data_t *curr, *head, *temp;
    char barcode [MAX_BARCODE_LENGTH];

    int threads=4;

    while (1) {
        int option_index = 0;
        //colon after a flag means should have arguments and no colon means just a flag i.e bool, no args after it
        optc = getopt_long (argc, argv, "dnucvof:r:b:z:w:m:s:l:z:a:t:", paired_long_options, &option_index);

        if (optc == -1) break;

        switch (optc) {
            if (paired_long_options[option_index].flag != 0) break;

            case 'f':
            fq1_fn = (char*) malloc (strlen (optarg) + 1);
            strcpy (fq1_fn, optarg);
            break;

            case 'r':
            fq2_fn = (char*) malloc (strlen (optarg) + 1);
            strcpy (fq2_fn, optarg);
            break;

            case 'b':
            barfn = (char*) malloc (strlen (optarg) + 1);
            strcpy (barfn, optarg);
            break;

            case 'z':
            if(unassigned1_fn) {
                free(unassigned1_fn);
            }
            unassigned1_fn = (char*) malloc (strlen (optarg) + 1);
            strcpy (unassigned1_fn, optarg);
            break;

            case 'w':
            if(unassigned2_fn) {
                free(unassigned2_fn);
            }
            unassigned2_fn = (char*) malloc (strlen (optarg) + 1);
            strcpy (unassigned2_fn, optarg);
            break;

            case 'c':
            params.combine=1;
            break;

            case 'u':
            params.umi=1;
            break;

            case 'm':
            params.mismatch = atoi (optarg);
            break;

            case 'l':
            params.min_umi_len = atoi (optarg);
            break;

            case 'a':
            params.max_5prime_crop = atoi (optarg);
            break;

            case 'n':
            params.no_comment = 1;
            break;

            case 't':
            threads = atoi (optarg);
            break;

            case 's':
            log_fn = (char*) malloc (strlen (optarg) + 1);
            strcpy (log_fn, optarg);
            break;

            case 'v':
            version(EXIT_SUCCESS);
            break;

	    //NOTE if user requrested the help menu i.e --help then
	    //return success for all other cases below while help menu
	    //is printed it wasn't intended by user (or at least we don't know that)
	    //and therefore exit code - fail
            case 'h':
            version(EXIT_SUCCESS);
            break;

            case 'o':
            little_story(EXIT_SUCCESS);
            break;

            case '?':
            usage(EXIT_FAILURE);
            break;

            default:
            usage(EXIT_FAILURE);
            break;
        }
    }

    params.fq1_fd = gzopen(fq1_fn, "r");
    params.fq2_fd = gzopen(fq2_fn, "r");
    params.unassigned1_fd = gzopen(unassigned1_fn, "wb");
    params.unassigned2_fd = gzopen(unassigned1_fn, "wb");
    params.umis_2_short_fd = fopen(umis_2_short_fn, "a");

    // ? where does this goes?
    fprintf(params.umis_2_short_fd, "name\tumi\tlen\tmin_len\n");

    //TODO plugin sanity_chk here
    
    if(params.fq2_fd) {
        params.paired = 1;
    }

    fprintf(stderr, "\n\
                     \n  Running: %s\
                     \n  Command line args:\
                     \n      --fq1 %s\
                     \n      --fq2 %s\
                     \n      --barcodes %s\
                     \n      --unassigned_R1 %s\
                     \n      --unassigned_R2 %s\
                     \n      --combine %d\
                     \n      --umi %d\
                     \n      --max-mismatch %d\
                     \n      --min-umi-len %d\
                     \n      --max-5prime-crop %d\
                     \n      --no-comment %d\
                     \n      --stats %s\
                     \n      --threads %d\
                     \n\
                     \n  In Progess...\
                     \n", PROGRAM_NAME,\
                     fq1_fn, fq2_fn,\
                     barfn,\
                     unassigned1_fn, unassigned2_fn,\
                     params.combine, params.umi,\
                     params.mismatch, params.min_umi_len,
                     params.max_5prime_crop, params.no_comment,\
                     log_fn, threads);

    char *bcout_fn1 = NULL;
    char *bcout_fn2 = NULL;
    /* Creating linked list of barcode data */
    // https://www.hackerearth.com/practice/data-structures/linked-list/singly-linked-list/tutorial/
    // where each node is represents one barcode from the barcode file
    head = NULL;
    while (fscanf (barfile, "%s%s", barcode, s_name) != EOF) {
        curr = (barcode_data_t*) malloc(sizeof(barcode_data_t));
        curr->bc = (char*) malloc(strlen(barcode) + 1);
        strcpy(curr->bc, barcode);

        bcout_fn1 = (char *) malloc(MAX_FILENAME_LENGTH*2);
        bcout_fn1[0] = '\0';
        get_bc_fn(&bcout_fn1, s_name, curr->bc, 1);
        //curr->bcfile1 = fopen (_mkdir(bcout_fn1), "w");
        curr->bcfile1 = gzopen(_mkdir(bcout_fn1), "wb");
        //curr->bcfile1 = popen(_mkdir(bcout_fn1), "wb");
        // popen returns file handler

        if(params.paired > 0 && params.combine < 0) {
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
    free(barfn);
    free(fq1_fn);
    free(fq2_fn);
    free(unassigned1_fn);
    free(unassigned2_fn);
    free(umis_2_short_fn);

    // Threading
    pthread_t tid[threads];
    pthread_mutex_t in_lock;
    pthread_mutex_t out_lock;
    pthread_cond_t cv;
    int line_num = 0;
    int out_line_num = 0;

    pthread_mutex_init(&in_lock, NULL);
    pthread_mutex_init(&out_lock, NULL);
    pthread_cond_init(&cv, NULL);

    thread_data_t thread_data[threads];

    for(int i=0; i < threads; i++) {

	thread_data->params = &params;
	thread_data->curr = &curr;
	thread_data->metrics = &metrics;
        thread_data->id = i;
        thread_data->in_lock = &in_lock;
        thread_data->out_lock = &out_lock;
        thread_data->line_num = &line_num;
        thread_data->out_line_num = &out_line_num;
        thread_data->cv = &cv;

        pthread_create(&(tid[i]), NULL, &demult_runner, thread_data);
    }

    for(int i=0; i < threads; i++) {
        pthread_join(tid[i], NULL);
    }
    printf("Threads all done\n");

    pthread_mutex_destroy(&in_lock);
    pthread_mutex_destroy(&out_lock);

    //if (!log_fn) { is this better?
    if (log_fn == NULL) {
        log_file = stdout;
    }
    else {
        log_file = fopen(log_fn, "w");
    }

    fprintf (log_file, "Barcode\tN_records\tN_pairs\tP_pairs\n");
    curr = head;
    int total_pairs = metrics.total/2;

    while (curr) {

        int n_pairs = curr->num_records/2;
        float percent_pairs = (float) n_pairs/total_pairs;

        fprintf (log_file,"%s\t%d\t%d\t%.2f\n", curr->bc, curr->num_records, n_pairs, percent_pairs);

        curr = curr->next;
    }

    int unknown_pairs = metrics.num_unknown/2;
    float percent_unknown = (float) unknown_pairs/total_pairs;
    float tot_chk = (float) total_pairs/total_pairs;

    fprintf(log_file, "unassigned\t%d\t%d\t%.2f\n", metrics.num_unknown, unknown_pairs, percent_unknown);
    fprintf(log_file, "total\t%d\t%d\t%.2f\n", metrics.total, total_pairs, tot_chk);

    end = time(NULL);
    fprintf(stderr, "\n All done :) \
                     \n It took %.2f minutes\n",
                     difftime(end, start)/60);

    // good read :)
    little_story(EXIT_SUCCESS);

    fclose(barfile);
    fclose(log_file);
    params_destroy(&params);

    free(log_fn);

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
