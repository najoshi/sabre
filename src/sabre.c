
#include "sabre.h"
#include "utils.h"
#include "usage.h"
#include "demultiplex.h"

FILE* my_fopen(const char* fname, int gz) {
    static char* compressor = NULL;

    if (gz) {

        if (!compressor) {
            // Guess whether to use pigz or gzip
            char tmp[100];
            FILE* fin = popen("pigz --version 2>&1", "r");
            char* str = fgets(tmp, 50, fin);
            int found = strncmp("pigz",tmp,4)==0;
            compressor = str && found ? "pigz -p 2" : "gzip";
        }

        char command[2048];
        sprintf(command, "%s > %s", compressor, fname);
        FILE* ret = popen(command, "w");
        return ret;
    } else {
        return fopen(fname, "w");
    }
}

int main(int argc, char *argv[]) {

    if(argc < 2) {
        usage(EXIT_SUCCESS);
    }

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
        {"gz-out", optional_argument, NULL, 'g'},
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
    metrics. num_unknown = 0;
    metrics.total = 0;

    //clock_t begin = clock();
    time_t start, end;
    start = time(NULL);

    char *fq1_fn=NULL;
    char *fq2_fn=NULL;

    FILE* bc_fd;
    char *bc_fn=NULL;

    char *unassigned1_fn=NULL;
    char *unassigned2_fn=NULL;
    char *umis_2_short_fn=strdup("umis_too_short.txt");

    FILE* log_file;
    char *log_fn=strdup("stats.txt");

    int optc;
    extern char *optarg;

    barcode_data_t *curr, *head, *temp;

    int threads=4;

    while (1) {
        int option_index = 0;
        //colon after a flag means should have arguments and no colon means just a flag i.e bool, no args after it
        optc = getopt_long (argc, argv, "dnucvogf:r:b:z:w:m:s:l:z:a:t:", paired_long_options, &option_index);

        if (optc == -1) break;

        switch (optc) {

            case 'f':
            fq1_fn = (char*) malloc (strlen (optarg) + 1);
            strcpy (fq1_fn, optarg);
            break;

            case 'r':
            fq2_fn = (char*) malloc (strlen (optarg) + 1);
            strcpy (fq2_fn, optarg);
            break;

            case 'b':
            bc_fn = (char*) malloc (strlen (optarg) + 1);
            strcpy (bc_fn, optarg);
            break;

            case 'z':
            unassigned1_fn = (char*) malloc (strlen (optarg) + 1);
            strcpy (unassigned1_fn, optarg);
            break;

            case 'w':
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

            case 'g':
            params.gz_out = 1;
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
            usage(EXIT_SUCCESS);
            break;

            case 'o':
            little_story(EXIT_SUCCESS);
            break;

            case '*':
            usage(EXIT_FAILURE);
            break;

            default:
            usage(EXIT_FAILURE);
            break;
        }
    }

    // TODO check that what's requre if not run usage

    params.fq1_fd = gzopen(fq1_fn, "r");
    if(!params.fq1_fd) {
        fprintf(stderr, "ERROR: Could not open input file R1 '%s'.\n", fq1_fn);
        exit(EXIT_FAILURE);
    }

    params.fq2_fd = gzopen(fq2_fn, "r");
    if(!params.fq2_fd) {
        fprintf(stderr, "ERROR: Could not open input file R2 '%s'.\n", fq2_fn);
        exit(EXIT_FAILURE);
    }

    if (!unassigned1_fn) {
        unassigned1_fn= params.gz_out ? strdup("unassigned_R1.fq.gz") : strdup("unassigned_R1.fq");
    }
    if (!unassigned2_fn) {
        unassigned2_fn= params.gz_out ? strdup("unassigned_R2.fq.gz") : strdup("unassigned_R2.fq");
    }
    params.unassigned1_fd = my_fopen(unassigned1_fn, params.gz_out);
    params.unassigned2_fd = my_fopen(unassigned2_fn, params.gz_out);
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
                     \n\
                     \n", PROGRAM_NAME,\
                     fq1_fn, fq2_fn,\
                     bc_fn,\
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
    bc_fd = fopen(bc_fn, "r");
    if (!bc_fd) {
        fprintf(stderr, "ERROR: Unable to barcode file\n");
        exit(EXIT_FAILURE);
    }

    head = NULL;
    curr = NULL;

    char line_buff[1024];
    int max_items = 6;
    while(fgets(line_buff, 1024, bc_fd)) {
        curr = (barcode_data_t*) calloc(1,sizeof(barcode_data_t));

        char *p = strtok(line_buff, "\t");
        char *s_name = strdup(p);

        p = strtok(NULL, "\t");
        curr->bc_grp = strdup(p);

        bcout_fn1 = (char *) malloc(MAX_FILENAME_LENGTH*2);
        bcout_fn1[0] = '\0';
        get_bc_fn(&bcout_fn1, s_name, curr->bc_grp, 1, params.gz_out);
        curr->bcfile1 = my_fopen(_mkdir(bcout_fn1), params.gz_out);

        if(params.paired > 0 && params.combine < 0) {
            bcout_fn2 = (char *) malloc(MAX_FILENAME_LENGTH*2);
            bcout_fn2[0] = '\0';
            get_bc_fn(&bcout_fn2, s_name, curr->bc_grp, 2, params.gz_out);
            curr->bcfile2 = my_fopen(_mkdir(bcout_fn2), params.gz_out);
        }

        //TODO for hardcode max limit of items in the barcodes file to 6
        curr->bc = calloc(max_items, sizeof(void*));

        int i=0;
        while(i <= max_items && (p = strtok(NULL, "\t\n"))) {
            // remove the token, new line char
            curr->bc[i] = strdup(p);
            fprintf(stdout, "  BC %s ", curr->bc[i]);
            i++;
        }
        fprintf(stdout, "\n");

        curr->num_records = 0;
        curr->next = head;
        head = curr;
    }

    free(bcout_fn1);
    free(bcout_fn2);
    free(bc_fn);
    free(fq1_fn);
    free(fq2_fn);
    free(unassigned1_fn);
    free(unassigned2_fn);
    //free(umis_2_short_fn);

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

        thread_data[i].params = &params;
        thread_data[i].curr = curr;
        thread_data[i].metrics = &metrics;
        thread_data[i].id = i;
        thread_data[i].in_lock = &in_lock;
        thread_data[i].out_lock = &out_lock;
        thread_data[i].line_num = &line_num;
        thread_data[i].out_line_num = &out_line_num;
        thread_data[i].cv = &cv;

        pthread_create(&(tid[i]), NULL, &demult_runner, thread_data+i);
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

    while(curr) {

        int n_pairs = curr->num_records/2;
        float percent_pairs = (float) n_pairs/total_pairs;

        fprintf(log_file,"%s\t%d\t%d\t%.2f\n", curr->bc_grp, curr->num_records, n_pairs, percent_pairs);

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


    fclose(bc_fd);
    fclose(log_file);
    params_destroy(&params);

    free(log_fn);

    curr = head;
    while (curr) {
        fclose(curr->bcfile1);
        if (curr->bcfile2)
            fclose(curr->bcfile2);

        free (curr->bc_grp);
        free (curr->bc);
        temp = curr;
        curr = curr->next;
        free (temp);
    }

    // good read :)
    little_story(EXIT_SUCCESS);
}
