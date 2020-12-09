
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

#include "demultiplex.h"

// Is m1 a better match than m2?  0 -> yes, use m1.  1 -> no, stick with m2
int better(match_ret_t m1, match_ret_t m2) {
    if (m1.cropped<0)
        return 0;
    if (m2.cropped<0)
        return 1;
    return (m1.cropped + m1.mismatches) < (m2.cropped+m2.mismatches);
}

// Is m the best possible match?
int best(match_ret_t m) {
    return m.cropped==0 && m.mismatches==0;
}

void write_out(const param_t *params, pthread_mutex_t *out_lock, metrics_t* metrics,
               match_ret_t best_match, barcode_data_t *best_bc,
               const char* actl_bc,
               fq_rec_t *fq_rec1,fq_rec_t *fq_rec2);

void* demult_runner(void *arg) {
    fq_rec_t *fq_rec1 = (fq_rec_t*) malloc(sizeof(fq_rec_t));
    fq_rec_t *fq_rec2 = (fq_rec_t*) malloc(sizeof(fq_rec_t));

    init_fq_rec(fq_rec1);
    init_fq_rec(fq_rec2);

    thread_data_t* thread_data = (thread_data_t*)arg;
    //int my_line_num;

    /* Get reads, one at a time */

    while(1) {

        // lock reading
        pthread_mutex_lock(thread_data->in_lock);

        //this is equivalent to if(false), which means this block
        //is always skipped, unless when there is an error/end of the file
        if(get_fq_rec(fq_rec1, thread_data->params->fq1_fd)) {
            // sanity check no more reads
            pthread_mutex_unlock(thread_data->in_lock);
            break;
        }

        if(thread_data->params->paired > 0) {
            if(get_fq_rec(fq_rec2, thread_data->params->fq2_fd)) {
                //error out there becuase if reached the end of the file
                //then we should hit first break, above, since the assumptions
                //that the files of equal length. If issues with R2 only this is an error
                fprintf (stderr, "\n\
                        \n ERROR: R2 file shorter than R1 file.\
                        \n Stopping here:\
                        \n %s\
                        \n",
                        fq_rec1->name);
                pthread_mutex_unlock(thread_data->in_lock);
                exit(1);
            }
        }

        // unlock reading
        // TODO this bit of code for ordered fastq files, implement later?
        //my_line_num = *(thread_data->line_num);
        //*thread_data->line_num += 1;
        pthread_mutex_unlock(thread_data->in_lock);

        // Store a copy of the barcode found in the read (including the mismatches)
        char *actl_bc = NULL;

        /* Step 1: Find matching barcode */
        match_ret_t best_match = {-1,-1};
        barcode_data_t *best_bc = NULL;
        for (barcode_data_t *curr = thread_data->curr; curr!=NULL; curr = curr->next) {
            for (int i=0; curr->bc[i]; i++) {
                match_ret_t mtch = chk_bc_mtch(curr->bc[i], fq_rec1->seq, thread_data->params->mismatch, thread_data->params->max_5prime_crop);
                if(better(mtch, best_match)) {
                    //found better match
                    best_match = mtch;
                    best_bc = curr;
                    if (actl_bc) free(actl_bc);
                    actl_bc = strndup( (fq_rec1->seq)+mtch.cropped, strlen(curr->bc[i]) );
                    if (best(best_match))
                        break;
                }
            }

            if (best(best_match))
                break;
        }

        /* Step 2: Write read out into barcode specific file */
        //TODO this bit of code to keep fastq files ordered as per original fastq files
        //which I don't think that needed? at least not at this stage
        // lock writing
        //while(*(thread_data->out_line_num) != my_line_num) {
        //    pthread_cond_wait(thread_data->cv, thread_data->out_lock);
        //}
        //*thread_data->out_line_num += 1;

        //pthread_cond_broadcast(thread_data->cv);  // Tell everyone it might be their turn!

        write_out(thread_data->params, thread_data->out_lock, thread_data->metrics, best_match, best_bc, actl_bc, fq_rec1, fq_rec2);
        if (actl_bc)
            free(actl_bc);

        thread_data->metrics->total += 2;
    }

    free(fq_rec1);
    free(fq_rec2);
    //free(thread_data); according to valgrind report this line isn't needed since no errors given out..
    return NULL;
}


void write_out(const param_t *params, pthread_mutex_t *out_lock, metrics_t* metrics,
               match_ret_t best_match, barcode_data_t *best_bc,
               const char* actl_bc,
               fq_rec_t *fq_rec1,fq_rec_t *fq_rec2) {

    char fqread1[MAX_READ_SIZE];
    char fqread2[MAX_READ_SIZE];

    fqread1[0] = '\0';
    fqread2[0] = '\0';

    char *umi_idx = NULL;

    if(best_bc != NULL) {
        //for now assume barcode and umi are in R1 read
        if(params->umi > 0) {

            const char *actl_umi_idx = (fq_rec1->seq)+strlen(actl_bc)+best_match.cropped;

            if(strlen(actl_umi_idx) < params->min_umi_len) {
                pthread_mutex_lock(out_lock);
                fprintf(params->umis_2_short_fd, "%s\t%s\t%zu\t%d\n", fq_rec1->name, actl_umi_idx, strlen(actl_umi_idx), params->min_umi_len);
                pthread_mutex_unlock(out_lock);
                return;
            }
            else {
                umi_idx = strdup(actl_umi_idx);
                umi_idx[params->min_umi_len] = '\0';
            }
        }

        if(params->combine > 0 && actl_bc != NULL) {
            get_merged_fqread(fqread1, fq_rec1, fq_rec2, actl_bc, umi_idx, params->no_comment, best_match.cropped);

            pthread_mutex_lock(out_lock);
            fputs(fqread1, best_bc->bcfile1);
            pthread_mutex_unlock(out_lock);

        }
        else {
            get_fqread(fqread1, fq_rec1, actl_bc, umi_idx, params->no_comment, best_match.cropped);

            pthread_mutex_lock(out_lock);
            fputs(fqread1, best_bc->bcfile1);

            if(params->paired > 0) {
                get_fqread(fqread2, fq_rec1, actl_bc, umi_idx, params->no_comment, best_match.cropped);

                fputs(fqread2, best_bc->bcfile2);

                //dont need to increment buff_cnt, assuming fq_read1 keeps the right count
                best_bc->num_records += 1;
            }
            pthread_mutex_unlock(out_lock);
        }
        best_bc->num_records += 1;
    }
    else {

        get_fqread(fqread1, fq_rec1, NULL, NULL, params->no_comment, 0);

        pthread_mutex_lock(out_lock);
        fputs(fqread1, params->unassigned1_fd);

        metrics->num_unknown += 1;

        if(params->paired > 0) {
            get_fqread(fqread2, fq_rec2, NULL, NULL, params->no_comment, 0);

            fputs(fqread2, params->unassigned2_fd);

            metrics->num_unknown += 1;
        }
        pthread_mutex_unlock(out_lock);
    }
}